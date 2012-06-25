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

#include "localIPReassignment.h"
#include "localIP.h"
#include "nextMachine.h"

#include <unistd.h>
#include <iostream>
#include <assert.h>

using namespace std;


LocalIPReassignment::LocalIPReassignment(const Instance& inst, const Solution& init_sol):
	sol(new Solution(init_sol)), 
	inst(inst), 
	init_sol(init_sol),
	nb_neighborhoods(inst.neighborhoods.size()),
    services_locations(new map< ul, list<const Process*> >[inst.num_services]),
	proc_in_neighborhood(new vector<const Process*>[inst.neighborhoods.size()])
{
	// extract the processes assigned to each neighborhood from init_sol
	for (ul idp = 0; idp < inst.num_processes; ++idp) {
		ul n = init_sol.assignment[idp]->neighborhood;
		proc_in_neighborhood[n].push_back(inst.processes[idp]);
	}

    spread();
}


LocalIPReassignment::~LocalIPReassignment()
{
	delete sol;
	delete [] proc_in_neighborhood;
    delete [] services_locations;
}


void LocalIPReassignment::run() {
    // For each neighborhood
	for (ul n = 0; n < nb_neighborhoods; ++n) {
        run_ip(n);
    }
}


map< const Machine*, vector<const Process*> >*
LocalIPReassignment::divide(ul neighborhood) {
    const vector<const Process*>& processes = proc_in_neighborhood[neighborhood];

    map< const Machine*, vector<const Process*> >* proc_on_machines =
        new map< const Machine*, vector<const Process*> >();
    vector<const Process*>::const_iterator it = processes.begin();
    for (; it != processes.end(); ++it) {
        (*proc_on_machines)[sol->assignment[(*it)->id]].push_back(*it);
    }
    
    return proc_on_machines;
}


bool LocalIPReassignment::trivial_assignment(
        const std::vector<const Machine*>& mc, const std::vector<const Process*>& proc) {
    if (mc.size() == 0) return true;
    if (mc.size() == 1) { // No movement can be done in this neighborhood
        // updates solution and ends
        set_sol(mc[0], proc);
        return true;
    }
    return false;
}

void LocalIPReassignment::set_sol(const Machine* const mc, const vector<const Process*>& proc) {
    vector<const Process*>::const_iterator it = proc.begin();

    while (it != proc.end()) {
        set_sol(mc, *it);
        ++it;
    }
}

void LocalIPReassignment::set_sol(const Machine* const mc, const Process* const proc) {
        // set assignment
        sol->assignment[proc->id] = mc;

        // update spreads
        if (mc != init_sol.assignment[proc->id]) {
            list<const Process*>& lst = services_locations[proc->service->id][init_sol.assignment[proc->id]->location];
            
            lst.remove(proc);
            if (lst.empty()) {
                map< ul, std::list<const Process*> >::iterator it;
                it = services_locations[proc->service->id].find(init_sol.assignment[proc->id]->location);
                services_locations[proc->service->id].erase(it);
            }
            services_locations[proc->service->id][mc->location].push_back(proc);
        }
}


void LocalIPReassignment::run_ip(const Machine* const m1, const Machine* const m2) {
    vector<const Machine*> m;
    m.push_back(m1);
    m.push_back(m2);
    run_ip(m);
}


void LocalIPReassignment::run_ip(const vector<const Machine*>& mch) {
    const vector<const Process*>& proc = proc_in_neighborhood[mch.front()->neighborhood];
    vector<const Process*> processes;

    const vector<const Machine*>::const_iterator beg = mch.begin();
    const vector<const Machine*>::const_iterator end = mch.end();
    vector<const Process*>::const_iterator it = proc.begin();
    for (; it != proc.end(); ++it) {
        if (find(beg, end, sol->assignment[(*it)->id]) != end) {
            processes.push_back(*it);
        }
    }

    run_ip(processes, mch);
}


void LocalIPReassignment::run_ip(ul neighborhood) {
    const vector<const Process*>& processes = proc_in_neighborhood[neighborhood];
    const vector<const Machine*>& machines = inst.neighborhoods[neighborhood];

    if (trivial_assignment(machines, processes)) return;

    map< const Machine*, vector<const Process*> >& proc_on_machines = *(divide(neighborhood));
    map< const Machine*, vector<const Process*> >::iterator it = proc_on_machines.begin();
    // Remove big machines
    while (it != proc_on_machines.end()) {
        if (it->second.size() > NTHRESHOLD) {
            set_sol(it->first, it->second);
            proc_on_machines.erase(it++);
        } else {
            ++it;
        }
    }

    NextMachine nm;
    nm.sort_by_load_costs(inst, proc_on_machines);
    //nm.no_sort(proc_on_machines);
    //nm.random_sort(42, proc_on_machines);
    
    // Greedy decomposition
    while (!nm.empty()) {
        vector<const Machine*> mc;
        vector<const Process*> pc;
        int num_machines = 0;
        ul num_processes = NTHRESHOLD;
        while (num_machines < MTHRESHOLD) {
            pair < const Machine *, vector < const Process *>*> p = nm.next(num_processes);
            if (p.first == NULL) break;
            num_processes -= p.second->size();
            ++num_machines;
            mc.push_back(p.first);
            pc.insert(pc.end(),p.second->begin(),p.second->end());
        }
        run_ip(pc, mc);
    }

    delete &proc_on_machines;
}

/** Computes spreads for init_sol */
void LocalIPReassignment::spread() {
    const Machine ** const assignment = init_sol.assignment;

    for (ul i = 0; i < inst.num_processes; ++i) {
        Process* proc = inst.processes[i];
        const Machine* mc = assignment[i];
        services_locations[proc->service->id][mc->location].push_back(proc);
    }
}


/** Returns sets to be considered for spread constraints in the subproblem */
list<Spread_constraint>* LocalIPReassignment::spread_constraints(const vector<const Process*>& processes,
        const vector<const Machine*>& machines) {
    set<const Service*> ok; // Services done
    list<Spread_constraint>* spr = new list<Spread_constraint>;
    
    vector<const Process*>::const_iterator it = processes.begin();
    for (; it != processes.end(); ++it) {
        const Process* proc = *it;
        if (ok.find(proc->service) != ok.end()) continue;

        const Service* serv = proc->service;
        ok.insert(serv);
        
        set<ul> current_locations;
        vector<const Process*>::const_iterator it2 = serv->processes.begin();
        for (; it2 != serv->processes.end(); ++it2) {
            if (find(processes.begin(), processes.end(), *it2) != processes.end())
                continue;
            current_locations.insert(sol->assignment[(*it2)->id]->location);
        }

        if (current_locations.size() >= serv->spreadmin) continue;

        Spread_constraint s(serv, serv->spreadmin - current_locations.size());
        vector<const Machine*>::const_iterator itm = machines.begin();
        for (; itm != machines.end(); ++itm) {
            if (current_locations.find((*itm)->location) != current_locations.end())
                continue;
            s.locations.insert((*itm)->location);
        }
        spr->push_back(s);
    }

    return spr;
}


void LocalIPReassignment::run_ip(const vector<const Process*>& processes,
        const vector<const Machine*>& machines) {

    if (trivial_assignment(machines, processes)) return;
    
    list<Spread_constraint>* spr = spread_constraints(processes, machines);
    LocalIP ip(inst, machines, processes, spr);
    
    // Solve
    ip.solve();
    
    // Updates solution
    const double * columns = ip.get_assignment();
    for (ul j = 0; j < processes.size(); ++j) {
        int index = j * machines.size();
        for (ul k = 0; k < machines.size(); ++k) {
            if (columns[index] > 0.5) { // == 1
                set_sol(machines[k], processes[j]);
                break;
            }
            ++index;
        }
    }

    delete spr;
}


#ifdef DEBUG_FUNC
/// Prints the list of processes in each neighborhood
void LocalIPReassignment::print_processes() {
    for (ul n = 0; n < nb_neighborhoods; ++n) {
	 	cout << endl << "Process in neighborhood " << n << endl;
	 	for (vector<const Process*>::iterator p = proc_in_neighborhood[n].begin(); 
	 		 p != proc_in_neighborhood[n].end(); ++p)
	 		cout << (*p)->id << " " ;
	}
}
#endif
