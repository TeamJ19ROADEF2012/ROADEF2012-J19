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
 
#include "instance.h"
#include "solution.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>

using namespace std;


Instance::~Instance() {
    delete [] is_transient;
    delete [] load_cost;

    for (ul m = 0; m < num_machines; ++m) {
        delete machines[m];
    }
    delete [] machines;

    for (ul p = 0; p < num_processes; ++p) {
        delete processes[p];
    }
    delete [] processes;

	for (ul s = 0; s < num_services; ++s) {
		delete services[s];
	}
	delete [] services;
    
    for (ul b = 0; b < num_bcosts; ++b) {
        delete bcosts[b];
    }
    delete [] bcosts;
}


Instance* parse_instance(char* file) {
    Instance * inst;
    ul tmp, tmp0;
    ul * rsc;
    ul num_res, num_mc, num_svc, num_proc, num_bc, pmv, svc, mmv;
    ifstream input(file);

    if (!input.is_open()) {
        cout << "Cannot open file " << file << endl;
        cout << "Aborting" << endl;
        exit(1);
    }

    // First pass - get const attributes
    input >> num_res;
    for (ul i = 0; i < num_res; ++i) input >> tmp >> tmp;

    input >> num_mc;
    // As much efficient as a single pass
    //for (ul i = 0; i <= num_mc; ++i) input.ignore(1000000,'\n');
    for (ul i = 0; i < num_mc; ++i) {
        // Ensures we are on a new line
        input >> tmp;
        // skip it
        //input.ignore(1000000,'\n');
        input >> tmp;
        // skip capacities and moving costs;
        for (ul j = 0; j < num_res; ++j) input >> tmp >> tmp;
        for (ul j = 0; j < num_mc; ++j) input >> tmp;
    }

    input >> num_svc;
    for (ul i = 0; i < num_svc; ++i) {
        // Skip spreadmin and dependencies
        input >> tmp >> tmp0;
        for (ul j = 0; j < tmp0; ++j) input >> tmp;
    }

    input >> num_proc;
    for (ul i = 0; i < num_proc; ++i) {
        // Skip service, requirements and pmc
        input >> tmp;
        for (ul j = 0; j < num_res; ++j) input >> tmp;
        input >> tmp;
    }

    input >> num_bc;
    for (ul i = 0; i < num_bc; ++i) {
        // Skip balanced costs
	    input >> tmp >> tmp >> tmp >> tmp;
    }

    input >> pmv >> svc >> mmv;

    inst = new Instance(num_res, num_mc, num_svc, num_proc,
            num_bc, pmv, svc, mmv);

    // Back to the beginning to parse data
    input.seekg(0, ios::beg);

    // Second pass - set attributes
    // Resources
    input >> tmp;
    for (ul i = 0; i < num_res; ++i) {
        input >> inst->is_transient[i] >> inst->load_cost[i];
    }

    // Machines
    input >> tmp;
    for (ul i = 0; i < num_mc; ++i) {
        input >> tmp; // neighborhood
        input >> tmp0; // location

		inst->machines[i] = new Machine(i, num_res, num_mc, tmp, tmp0);

        // Add the machine to its neighborhood (+resize if needed)
        if (tmp >= inst->neighborhoods.size())
            inst->neighborhoods.resize(tmp+1);
        inst->neighborhoods[tmp].push_back(inst->machines[i]);
        
        // Add the machine to its location (+resize if needed)
        if (tmp0 >= inst->locations.size())
            inst->locations.resize(tmp0+1);
        inst->locations[tmp0].push_back(inst->machines[i]);

        // Filling capacities
        for (ul j = 0; j < num_res; ++j)
            input >> inst->machines[i]->capacities[j];
        // Filling safety capacities
        for (ul j = 0; j < num_res; ++j)
            input >> inst->machines[i]->sf_capacities[j];
        // Filling moving costs in both places
        for (ul j = 0; j < num_mc; ++j) {
            input >> tmp;
            inst->machines[i]->mv_costs[j] = tmp * mmv;
        }
    }

    // Services
	for (ul s = 0; s < num_svc; ++s)
		inst->services[s] = new Service(s);

    input >> tmp;
    for (ul i = 0; i < num_svc; ++i) {
		input >> tmp0 >> tmp;
		inst->services[i]->spreadmin = tmp0;

        // set dependencies
        for (ul j = 0; j < tmp; ++j) {
            input >> tmp0;
			inst->services[i]->dependencies.push_back(inst->services[tmp0]);
			inst->services[tmp0]->rev_dependencies.push_back(inst->services[i]);
        }
    }

    // Processes
    input >> tmp;
    rsc = new ul[num_res];
    for (ul i = 0; i < num_proc; ++i) {
        input >> tmp;
		//inst->services[tmp]->processes.push_back(inst->processes[i]);
        for (ul j = 0; j < num_res; ++j) input >> rsc[j];
        input >> tmp0;
        inst->processes[i] = new Process(i, inst->services[tmp], num_res, num_mc, rsc, tmp0 * pmv);
		inst->services[tmp]->processes.push_back(inst->processes[i]);
    }
    delete[] rsc;

	// More computations to set no_conflict
	for (ul s = 0; s < num_svc; ++s) {
		assert(inst->services[s]->processes.size() != 0);
		inst->services[s]->no_conflict = ( inst->services[s]->processes.size() == 1 );
	}

    // balanced costs
    input >> tmp;
    for (ul i = 0; i < num_bc; ++i) {
        // Skip balanced costs
        input >> tmp >> tmp0 >> num_res >> num_mc;
        inst->bcosts[i] = new Balanced_cost(i, tmp, tmp0, num_res, num_mc);
    }

    input.close();
    return inst;
}


void Instance::set_init_assignment (Solution& init_sol) {
    for (ul i = 0; i < num_processes; ++i) {
        processes[i]->init_assignment = init_sol.assignment[i];
    }
}


#ifdef DEBUG
void Instance::print_model() {
	ul tmp;

	cout << num_resources << endl;
	for (ul i = 0; i < num_resources; ++i) {
		cout << is_transient[i] << " ";
		cout << load_cost[i] << endl;
	}

	cout << num_machines << endl;
	for (ul i = 0; i < num_machines; ++i) {
		tmp = machines[i]->neighborhood;
		cout << tmp << " ";
		vector<const Machine*>& v = neighborhoods[tmp];
		if (find(v.begin(), v.end(), machines[i]) == v.end()) {
			cout << "Error on neighboorhoods" << endl;
		} 
		tmp = machines[i]->location;
		cout << tmp << " ";
		vector<const Machine*>& v1 = locations[tmp];
		if (find(v1.begin(), v1.end(), machines[i]) == v1.end()) {
			cout << "Error on locations" << endl;
		}
		for(ul j = 0; j < num_resources; ++j) {
			cout << machines[i]->capacities[j] << " ";
		}
		for(ul j = 0; j < num_resources; ++j) {
			cout << machines[i]->sf_capacities[j] << " ";
		}
		for(ul j = 0; j < num_machines; ++j) {
			tmp = machines[i]->mv_costs[j];
			if (tmp != get_cost(i,j)) {
				cout << "Error on machine move costs" << endl;
			}
			cout << tmp / w_mmv << " ";
		} cout << endl;
	}

	cout << num_services << endl;
	for (ul i = 0; i < num_services; ++i) {
		cout << services[i]->spreadmin << " ";
		cout << services[i]->dependencies.size() << " ";
		for (ul j = 0; j < services[i]->dependencies.size(); ++j) {
			cout << services[i]->dependencies[j] << " ";
		} cout << endl;

	}

	cout << num_processes << endl;
	for (ul i = 0; i < num_processes; ++i) {
		cout << processes[i]->service->id << " ";
		for(ul j = 0; j < num_resources; ++j) {
			cout << processes[i]->resources[j] << " ";
		}
		cout << processes[i]->mv_cost / w_pmv << " " << endl;
	}

	cout << num_bcosts << endl;
	for (ul i = 0; i < num_bcosts; ++i) {
		cout << bcosts[i]->r0 << " " << bcosts[i]->r1 << " " <<
			bcosts[i]->target << endl << bcosts[i]-> weight << endl;
	}

	cout << w_pmv << " " << w_svc << " " << w_mmv << endl;
}

void Instance::print_properties() {
	// DEBUG
	cout << "Some instance caracteritics " 
		 << "--------------------------------" << endl;

	cout << "Num machines = " << num_machines << endl;
	cout << "Num processes = " << num_processes << endl;
	cout << "Num resources = " << num_resources;
	ul num_transient = 0;
	for (ul r = 0; r < num_resources; ++r)
		if (is_transient[r])
			num_transient++;
	cout << ", num transient = " << num_transient << endl;

	cout << "Num locations = " << locations.size() << endl;
	cout << "Num services = " << num_services << endl;

	ul num_dep = 0;
	for(ul s = 0; s < num_services; ++s)
		num_dep += services[s]->dependencies.size();
	cout << "Num dependencies = " << num_dep << endl;
	cout << "Num neighborhoods = " << neighborhoods.size() << endl;
	cout << "Num balance costs = " << num_bcosts << endl;

	cout << endl;
	for (ul n = 0; n < neighborhoods.size(); ++n) {
		cout << "Neighborhood " << n << endl;
		for (vector<const Machine*>::const_iterator it = neighborhoods[n].begin(); 
			 it != neighborhoods[n].end(); ++it)
			cout << (*it)->id << " " ;
		cout << endl;
	}
	
	cout << endl;
	for (ul l = 0; l < locations.size(); ++l) {
		cout << "Location " << l << endl;
		for (vector<const Machine*>::const_iterator it = locations[l].begin(); 
			 it != locations[l].end(); ++it)
			cout << (*it)->id << " " ;
		cout << endl;
	}
	cout << "--------------------------------------------------------" 
		 << endl;

}
#endif
