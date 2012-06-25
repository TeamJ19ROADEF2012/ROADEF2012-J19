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

#include "solution.h"

#include <fstream>
#include <iostream>
#include <assert.h>

using namespace std;

Solution::Solution(char* file, const Instance& inst):
	tmp_sol(NULL),
	inst(inst),
    assignment(new Machine const *[inst.num_processes]),
	services_of_machine(new set<const Service *>[inst.num_machines]),
	locations_of_service(new map<ul, ul>[inst.num_services]),
	services_of_neighborhood(new map<ul, ul>[inst.neighborhoods.size()]),
    moved_services(new ul[inst.num_services]), max_moved_services(pair<ul, ul>(0,0))
{
	pthread_mutex_init(&plock, NULL);

    ifstream input(file);
	ul m;

    if (!input.is_open()) {
        cout << "Cannot open file " << file << endl;
        cout << "Aborting" << endl;
        exit(1);
    }

	// Initialization of services_of_neighborhood
	for (ul n = 0; n < inst.neighborhoods.size(); ++n)
		for (ul s = 0; s < inst.num_services; ++s)
			services_of_neighborhood[n][s] = 0;


	pair<map<ul, ul>::iterator, bool> loc_info;

	// Update for each assignment
    for (ul p = 0; p < inst.num_processes; ++p) {
		const Process & proc = *inst.processes[p];
        input >> m;
		assert(m < inst.num_machines);
		assignment[p] = inst.machines[m];

		services_of_machine[m].insert(proc.service);

		// Ajouter une location au service, si la location est déjà présente, incrémenter le bnre de processus
		loc_info = locations_of_service[proc.service->id].insert( pair<ul, ul>(inst.machines[m]->location, 1) );
		if (loc_info.second == false) loc_info.first->second += 1;

        services_of_neighborhood[inst.machines[m]->neighborhood][proc.service->id] += 1;
    }

    input.close();
}


void Solution::init(const Solution* init_sol) {
    for (ul i = 0; i < inst.num_services; ++i) {
        locations_of_service[i].clear();
    }

    for (ul i = 0; i < inst.num_machines; ++i) {
		services_of_machine[i].clear();
    }

	// Initialization of services_of_neighborhood
	for (ul n = 0; n < inst.neighborhoods.size(); ++n)
		for (ul s = 0; s < inst.num_services; ++s)
			services_of_neighborhood[n][s] = 0;

	pair<map<ul, ul>::iterator, bool> loc_info;

	// Update for each assignment
    for (ul p = 0; p < inst.num_processes; ++p) {
		const Process & proc = *inst.processes[p];
        ul m = assignment[p]->id;
		assert(m < inst.num_machines);

		services_of_machine[m].insert(proc.service);

		// Ajouter une location au service, si la location est déjà présente, incrémenter le bnre de processus
		loc_info = locations_of_service[proc.service->id].insert( pair<ul, ul>(inst.machines[m]->location, 1) );
		if (loc_info.second == false) loc_info.first->second += 1;

        services_of_neighborhood[inst.machines[m]->neighborhood][proc.service->id] += 1;
    }

    compute_cost(init_sol);
}


Solution::~Solution() {
    delete [] assignment;
    delete [] services_of_machine;
	delete [] locations_of_service;
    delete [] services_of_neighborhood;
    delete [] moved_services;
}


ul Solution::compute_cost(const Solution* init_sol, ul** store_balanced_resources) {
	
    cost = 0;
    ul tmp, t0, t1;
    ul resources[inst.num_machines][inst.num_resources];
    ul * res;
    const ul * const lc = inst.load_cost;
    Machine * mc;
    Process * pr, ** processes = inst.processes;
    Balanced_cost * bc;


    // Load cost and set resources
    for (ul i = 0; i < inst.num_machines; ++i) {
        mc = inst.machines[i];
        res = resources[i];
        memset(res, 0, sizeof(ul)*inst.num_resources);
        for (ul j = 0; j < inst.num_processes; ++j) {
            if (assignment[j] != mc) continue;
            pr = processes[j];
            for (ul k = 0; k < inst.num_resources; ++k) {
                res[k] += pr->resources[k];
            }
        }
        for (ul k = 0; k < inst.num_resources; ++k) {
            if (res[k] > mc->sf_capacities[k]) {
				cost += lc[k] * (res[k] - mc->sf_capacities[k]);
			}
        }
    }
	
    // Balance cost
    bool store = (store_balanced_resources != NULL);
    for (ul j = 0; j < inst.num_bcosts; ++j) {
        bc = inst.bcosts[j];
        tmp = 0;
        for (ul i = 0; i < inst.num_machines; ++i) {
            res = resources[i];
            mc = inst.machines[i];
            t0 = mc->capacities[bc->r0] - res[bc->r0];
            t1 = mc->capacities[bc->r1] - res[bc->r1];
            if (bc->target * t0 > t1) {
                t0 = bc->target * t0 - t1;
                if (store)
                    store_balanced_resources[i][j] = t0;
                tmp += t0;
			} else if (store) {
                store_balanced_resources[i][j] = 0;
            }
        }
        cost += tmp * bc->weight;
    }

    if (init_sol == NULL) return cost;

    // Process and machine move costs
    memset(moved_services, 0, inst.num_services * sizeof(ul));
    for (ul i = 0; i < inst.num_processes; ++i) {
        const Machine * const m = init_sol->assignment[i];
        if (assignment[i] != m) {
            pr = processes[i];
            moved_services[pr->service->id] += 1;
            cost += pr->mv_cost; 
            cost += inst.get_cost(m, assignment[i]);
        }
    }

    // Service move cost
    max_moved_services.first = 0;
	max_moved_services.second = 0;
    for (ul i = 0; i < inst.num_services; ++i) {
        if (moved_services[i] > max_moved_services.first) {
            max_moved_services.first = moved_services[i];
			max_moved_services.second = 1;
		}
		else if (moved_services[i] == max_moved_services.first)
			max_moved_services.second += 1;
    }

    cost += max_moved_services.first * inst.w_svc;

    return cost;
}


void Solution::print(const char* file) {
    ostream* fp = &cout;
    ofstream fout(file);
    const Solution* sol = this;

    pthread_mutex_lock(&plock);
    // verrou obtenu

    // sinon verrou non obtenu => operator= bloqué dans son mutex
    // on travaille donc sur tmp_sol et on relâchera plock car le flot
    // dans le mutex a été interrompu par le signal

    if (!fout.is_open()) {
        cout << "Cannot write output file " << file << endl;
        cout << "Printing solution to the standard output" << endl;
    } else fp = &fout;

    ostream& out = *fp;

    for (ul i = 0; i < inst.num_processes; ++i) {
        out << sol->assignment[i]->id << " ";
    }

    if (fp == &fout) fout.close();

    pthread_mutex_unlock(&plock);
}

void Solution::apply_move_assignment(const Process& proc, const Machine* const m) {
    // We are modifying an assignment
    const Machine * const oldm = assignment[proc.id];
    //Machine& oldm = *(inst.machines[old_assign]);
    map<ul,ul>::iterator it2;

    // Updating moved_services
    if (oldm == proc.init_assignment) {
            
		moved_services[proc.service->id] += 1;
            
		if (moved_services[proc.service->id] > max_moved_services.first) {
            max_moved_services.first += 1; // the gap cannot be greater than 1
			max_moved_services.second = 1; // only this service reaches te max
		}
		else if (moved_services[proc.service->id] == max_moved_services.first)
			max_moved_services.second += 1;

    } else if (m == proc.init_assignment) {
		assert(moved_services[proc.service->id] >= 1);
        moved_services[proc.service->id] -= 1;
		if (moved_services[proc.service->id] == max_moved_services.first - 1) {
			if (max_moved_services.second == 1) {
				max_moved_services.first -= 1;
				max_moved_services.second = 0;
				for (ul s = 0; s < inst.num_services; ++s) {
					if (moved_services[s] == max_moved_services.first)
						max_moved_services.second += 1;
				}
			}else
				max_moved_services.second -= 1;
		}
    }

	assert(locations_of_service[proc.service->id][oldm->location] >= 1);
	if (locations_of_service[proc.service->id][oldm->location] == 1)
		locations_of_service[proc.service->id].erase(oldm->location);
	else 
		locations_of_service[proc.service->id][oldm->location] -= 1;

    // Updates neighborhoods
    // it2 = services_of_neighborhood[oldm->neighborhood].find(proc.service);
    // it2->second -= 1;
    // if (it2->second == 0)
    //     services_of_neighborhood[oldm->neighborhood].erase(it2);

	assert(services_of_neighborhood[oldm->neighborhood][proc.service->id] >= 1);
	services_of_neighborhood[oldm->neighborhood][proc.service->id] -= 1;
}


void Solution::apply_swap(const Process& proc1, const Machine* const m1,
            const Process& proc2, const Machine* const m2) {
    pthread_mutex_lock(&plock);

    if (proc1.service != proc2.service) {
        services_of_machine[assignment[proc1.id]->id].erase(proc1.service);
        services_of_machine[assignment[proc2.id]->id].erase(proc2.service);
    }

    apply_move_assignment(proc1, m1);
    apply_assignment(proc1, m1, false, false);
    apply_move_assignment(proc2, m2);
    apply_assignment(proc2, m2, false, false);

    pthread_mutex_unlock(&plock);
}

void Solution::apply_assignment(const Process& proc, const Machine* const m, bool move, bool lock) {
    if (lock)
        pthread_mutex_lock(&plock);

    if (move) {
        // The following line is ok iff the previous solution respects conflict constraints
	    assert(services_of_machine[assignment[proc.id]->id].find(proc.service) !=
                services_of_machine[assignment[proc.id]->id].end());
        services_of_machine[assignment[proc.id]->id].erase(proc.service);

        apply_move_assignment(proc, m);
    }

	assignment[proc.id] = m;
	services_of_machine[m->id].insert(proc.service);
    services_of_neighborhood[m->neighborhood][proc.service->id] += 1;
	pair<map<ul, ul>::iterator, bool> loc_info;
	loc_info = locations_of_service[proc.service->id].insert(pair<ul, ul>(m->location, 1));
	if (loc_info.second == false) loc_info.first->second += 1;
    
    if (lock)
        pthread_mutex_unlock(&plock);
}


/** Returns all processes assigned to the given neighborhood */
vector<const Process*>* Solution::neighborhood_processes(ul neighborhood) {
    vector<const Process*>* processes = new std::vector<const Process*>;
    for (int i = 0; i < (int)inst.num_processes; ++i) {
        if (assignment[i]->neighborhood == neighborhood) {
            processes->push_back(inst.processes[i]);
        }
    }
    return processes;
}

