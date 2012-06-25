/*
 * Machine reassignment solver for the ROADEF/EURO 2012 Challenge
 * Copyright (C) 2012 Michaël Gabay, Sofia Zaourar
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include "localIP.h"

#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinModel.hpp"
#include "CbcModel.hpp"
#include "CbcHeuristicFPump.hpp"
#include "ClpPresolve.hpp"

#include <math.h>
#include <set>

#include <fstream>

using namespace std;


void set_solve_param(CbcModel& model) {
    CbcHeuristicFPump fpump(model);
    model.addHeuristic(&fpump);
    //model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);
	//model.setDblParam(CbcModel::CbcMaximumSeconds, 60.0);
    model.solver()->messageHandler()->setLogLevel(0);
    model.setLogLevel(0);
}

/// Destructor
LocalIP::~LocalIP() {
    delete mod;
}

/// Résolution
double LocalIP::solve() {
    // Creates model
    CoinModel& model = *build_model();

    // Creates solver independant interface using Clp
    OsiClpSolverInterface solver1;
    OsiSolverInterface *solver = &solver1;
    
    // Loads model
    solver->loadFromCoinModel(model);
    
    ClpSolve solvectl;
    solvectl.setSolveType(ClpSolve::useDual);
    solvectl.setPresolveType(ClpSolve::presolveOn);
    solver1.setSolveOptions(solvectl);
    
    mod = new CbcModel(*solver);
    set_solve_param(*mod);
    mod->initialSolve(); // Solve LP
    mod->branchAndBound();

    delete &model;
    return 0.0;
}


CoinModel* LocalIP::build_model() {
    CoinModel* mod = new CoinModel();
    CoinModel& model = *mod;

    // Computes number of variables
    int process_variables = machines.size() * processes.size();
    int resources_variables = machines.size() * instance.num_resources;
    int load_cost_variables = resources_variables;
    int balance_cost_variables = instance.num_bcosts * machines.size();
    int service_move_costs_variables = 0;//instance.num_services + 1;
    int spread_variables = 0;
    // count spread variables :
    if (spread_list != NULL) {
        list<Spread_constraint>::iterator sit = spread_list->begin();
        for (;sit != spread_list->end(); ++sit) {
            spread_variables += (*sit).locations.size();
        }
    }

    int num_variables = process_variables + resources_variables +
        load_cost_variables + balance_cost_variables +
        service_move_costs_variables + spread_variables;

    // Creates all variables and sets bounds for the assignment variables
    model.setColumnBounds(num_variables - 1, 0.0, 1.0);
    
    // Define variables attributes
    int index = 0;
    int bound = process_variables;
    vector<const Process*>::const_iterator pit = processes.begin();
    vector<const Machine*>::const_iterator mit = machines.begin();
    // process variables p(process,machine)
    for (; pit != processes.end(); ++pit) {
        // a process is assigned to 1 machine
        // variable will be 1 if assigned, 0 otherwise

        // process move cost
        ul pmc = (*pit)->mv_cost;
        // machine move costs for this process
        ul * const mmc = (*pit)->init_assignment->mv_costs;
        for (mit=machines.begin(); mit != machines.end(); ++mit) {
            // Variable cost
            if (*mit == (*pit)->init_assignment) { // No move
                model.setObjective(index, 0.0);
            } else { // Process moved
                model.setObjective(index,pmc+mmc[(*mit)->id]);
            }
            model.setInteger(index); // integer variable
            ++index;
        }
    }

    // resources variables
    // U(machine,resource) = resource used on machine
    assert(index == bound);
    bound += resources_variables;
    mit = machines.begin();
    for (;mit != machines.end(); ++mit) {
        for (int j = 0; j < (int)instance.num_resources; ++j) {
            model.setColumnBounds(index, 0.0, (double)(*mit)->capacities[j]);
            model.setObjective(index, 0.0);         // do not infer in the objective
            ++index;
        }
    }

    // Load cost
    assert(index == bound);
    bound += load_cost_variables;
    // We add one machine and its resources before the next one...
    for(int r = 0; index < bound; ++index) {
        model.setColumnBounds(index, 0.0, COIN_DBL_MAX);
        model.setObjective(index, instance.load_cost[r]);
        r = (r + 1) % instance.num_resources;
    }


    // Balance cost
    assert(index == bound);
    bound += balance_cost_variables;
    // We add one machine and its balanced_costs before the next one...
    for(int b = 0; index < bound; ++index) {
        Balanced_cost& bc = *instance.bcosts[b];
        model.setColumnBounds(index, 0.0, COIN_DBL_MAX);
        model.setObjective(index, bc.weight);
        b = (b + 1) % instance.num_bcosts;
    }


    // Spread variables
    assert(index == bound);
    bound += spread_variables;
    if (spread_list != NULL) {
        list<Spread_constraint>::iterator sit = spread_list->begin();
        for (;sit != spread_list->end(); ++sit) {
            Spread_constraint& spc = (*sit);
            for (int j = 0; j < (int)spc.locations.size(); ++j) {
                model.setColumnBounds(index, 0.0, 1.0);
                model.setObjective(index, 0.0);         // do not infer in the objective
                model.setInteger(index); // integer variable
                ++index;
            }
        }
    }


    /*
    // Service move costs
    assert(index == bound);
    bound += service_move_costs_variables;
    // Number of process from a service moved
    for(; index < bound - 1; ++index) {
        model.setColumnBounds(index, 0.0, COIN_DBL_MAX);
        model.setObjective(index, 0.0);
    }
    model.setColumnBounds(index, 0.0, COIN_DBL_MAX);
    model.setObjective(index, instance.w_svc);
*/
    //============== Constraints ==============//
    int len = machines.size() * processes.size() + 1; // ok because the parameters are small !
    int * columns = new int[len];
    double * elements = new double[len];

    // One process is assigned to one machine
    index = 0;
    for (ul i = 0; i < processes.size(); ++i) {
        for (ul j = 0; j < machines.size(); ++j) {
            columns[j] = index++;
            elements[j] = 1.0;
        }
        model.addRow(machines.size(), columns, elements, 1.0, 1.0);
    }

    
    // Resource variables and transience constraints
    bound = process_variables - 1; // Before the first resource variable
    for (int i = 0; i < (int)machines.size(); ++i) {
        const Machine* const current_machine = machines[i];
        for (int r = 0; r < (int)instance.num_resources; ++r) {
            // Set usage and compute transience variables
            bool tr = instance.is_transient[r];
            ul capa = 0;
            vector<int> init;
            index = i;
            for (int j = 0; j < (int)processes.size(); ++j) {
                const Process* const p = processes[j];
                columns[j+1] = index;
                elements[j+1] = (double)(p->resources[r]);
                if (tr && p->init_assignment == current_machine) {
                    //capa += p->resources[r];
                    init.push_back(j);
                }
                index += machines.size();
            }
            if (tr) {
                for (int j = 0; j < (int)instance.num_processes; ++j) {
                    const Process* const p = instance.processes[j];
                    if (p->init_assignment == current_machine) {
                        capa += p->resources[r];
                    }
                }
            }
            ++bound;
            columns[0] = bound; elements[0] = -1.0;
            // Because of balance costs, the following constraint needs to be an equality!
            model.addRow(processes.size()+1, columns, elements, 0.0, 0.0);

            // Recall that resources constraints were added into variables bounds
            // Only transient resources need to be considered here 
            if (!tr) continue;
            elements[0] = 1.0;
            for (int j = 0; j < (int)init.size(); ++j) {
                columns[j+1] = init[j] * machines.size() + i;
                elements[j+1] = -(double)processes[init[j]]->resources[r];
            }
            model.addRow(init.size()+1, columns, elements, -COIN_DBL_MAX, current_machine->capacities[r] - capa);
        }
    }
    
    
    // Conflict constraints
    // Only added when necessary
    set<const Service*>* services_done = new set<const Service*>();
    index = 0;
    for (pit = processes.begin(); pit != processes.end(); ++index, ++pit) {
        if (services_done->find((*pit)->service) != services_done->end())
            continue;

        // Service not added yet
        int num = 0;
        const Service* serv = (*pit)->service;
        for (int j = index + 1; j < (int)processes.size(); ++j) {
            if (processes[j]->service == serv) {
                columns[++num] = j * machines.size();
                elements[num] = 1.0;
            }
        }
        if (num == 0) continue; // there is only one process in this service

        services_done->insert(serv);
        ++num;
        columns[0] = index * machines.size();
        elements[0] = 1.0;
        model.addRow(num, columns, elements, -COIN_DBL_MAX, 1.0);
        for (int j = 1; j < (int)machines.size(); ++j) {
            for (int k = 0; k < num; ++k) {
                columns[k] += 1; // next machine
            }
            model.addRow(num, columns, elements, -COIN_DBL_MAX, 1.0);
        }
    }
    delete services_done;


    // Load costs
    index = process_variables; // first resource variable
    bound = process_variables + resources_variables; // first load cost
    for (int i = 0; i < (int)machines.size(); ++i) {
        for (int r = 0; r < (int)instance.num_resources; ++r) {
            columns[0] = index++; columns[1] = bound++;
            elements[0] = -1.0; elements[1] = 1.0;
            model.addRow(2, columns, elements, -1.0 * ((double)machines[i]->sf_capacities[r]), COIN_DBL_MAX);
        }
    }


    // Balance costs
    bound = process_variables + resources_variables + load_cost_variables;
    index = process_variables;
    for (int i = 0; i < (int)instance.num_bcosts; ++i) {
        Balanced_cost& bc = *(instance.bcosts[i]);
        for (int j = 0; j < (int)machines.size(); ++j) {
            // balance cost variables
            columns[0] = bound + j * instance.num_bcosts + i; 
            elements[0] = 1.0;
            // + target * U(m,r1)
            columns[1] = index + j * instance.num_resources + bc.r0; 
            elements[1] = (double)bc.target;
            // - U(m,r2)
            columns[2] = index + j * instance.num_resources + bc.r1; 
            elements[2] = -1.0;

            double lower = (double)bc.target * (double)machines[j]->capacities[bc.r0] -
                (double)machines[j]->capacities[bc.r1];
            model.addRow(3, columns, elements, lower, COIN_DBL_MAX);
        }
    }
   

    // Spread constraints
    index = num_variables - spread_variables;
    if (spread_list != NULL) {
        list<Spread_constraint>::iterator sit = spread_list->begin();
        for (;sit != spread_list->end(); ++sit) {
            Spread_constraint& spc = (*sit);
            
            // value of the spread constraint on a location : Sum x_ij >= spread_loc
            // for the service and the location
            set<ul>::iterator current_loc = spc.locations.begin();
            for (int j = 0; current_loc != spc.locations.end(); ++j, ++current_loc) {
                // j is the location, sum all x_ij on the machines in this location
                columns[0] = index + j;
                elements[0] = -1.0;
                int num = 1;
                ul loc = *current_loc;
                for (int k = 0; k < (int)machines.size(); ++k) {
                    if (machines[k]->location != loc) continue;
                    for (int l = 0; l < (int)processes.size(); ++l) {
                        if (processes[l]->service != spc.svc) continue;
                        columns[num] = k + l * machines.size();
                        elements[num] = 1.0;
                        ++num;
                    }
                }
                model.addRow(num, columns, elements, 0, COIN_DBL_MAX);
            }

            // spread constraint : Sum spread(loc) >= spreadmin
            for (int j = 0; j < (int)spc.locations.size(); ++j) {
                columns[j] = index++;
                elements[j] = 1.0;
            }
            if (spc.locations.size() < spc.required) PRINT(spc.locations.size()<< " ===== " << spc.required);
            model.addRow(spc.locations.size(), columns, elements, spc.required, COIN_DBL_MAX);
        }
    }
    

    delete [] columns;
    delete [] elements;

    return mod;
}
