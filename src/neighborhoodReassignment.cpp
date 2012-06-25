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

#include "neighborhoodReassignment.h"
#include "vbp_heuristic.h"

#include <unistd.h>
#include <iostream>
#include <assert.h>

using namespace std;


NeighborhoodReassignment::NeighborhoodReassignment(const Instance& inst, const Solution& init_sol):
	sol(new Solution(init_sol)), 
	inst(inst), 
	init_sol(init_sol),
	nb_neighborhoods(inst.neighborhoods.size()),
    services_locations(new map< ul, list<const Process*> >[inst.num_services]),
	real_cap_in_neighborhhod(new map<const Machine *, ul *>[inst.neighborhoods.size()]),
	proc_in_neighborhood(new vector<const Process*>[inst.neighborhoods.size()])
{
	compute_real_cap();
	// extract the processes assigned to each neighborhood from init_sol
	for (ul idp = 0; idp < inst.num_processes; ++idp) {
		ul n = init_sol.assignment[idp]->neighborhood;
		proc_in_neighborhood[n].push_back(inst.processes[idp]);
	}
}


NeighborhoodReassignment::~NeighborhoodReassignment()
{
	delete sol;
	delete [] proc_in_neighborhood;
    delete [] services_locations;
	for (ul n = 0; n < inst.neighborhoods.size(); ++n) {
		map<const Machine *, ul *>::iterator itm = real_cap_in_neighborhhod[n].begin();
		for (; itm != real_cap_in_neighborhhod[n].end(); ++itm) {
			delete [] itm->second;
		}
	}
	delete [] real_cap_in_neighborhhod;
}


void NeighborhoodReassignment::run() {
    // print_processes();
	
    // For each neighborhood ( = machines)
	for (ul mn = 0; mn < nb_neighborhoods; ++mn) {
	 	cout << "Neighborhood " << mn << ": " << endl;
	 	// For each set of processes (associated with the neighborhoods)
		for (ul pn = 0; pn < nb_neighborhoods; ++pn) {
	 		cout << "Procs of " << pn << ": " << endl;
	 		run_vbp(mn, pn);
	 	}
	}
}

void NeighborhoodReassignment::run_vbp(ul machines_neighborhood, ul processes_neighborhood)
{
	// Call a vbp heuristic
	Vbp_heuristic * vbp_h = new Vbp_heuristic(inst, inst.neighborhoods[machines_neighborhood], 
											  proc_in_neighborhood[processes_neighborhood],
											  real_cap_in_neighborhhod[machines_neighborhood]);
	vbp_h->run();
	
	// Add the new assignments to sol
	// for (vector<const Process*>::iterator p = proc_in_neighborhood[processes_neighborhood].begin(); 
	// 	 p != proc_in_neighborhood[processes_neighborhood].end(); ++p) {
	// 	sol->assignment[(*p)->id] = vbp_h->assignment[*p];
	// }

	delete vbp_h;
}

void NeighborhoodReassignment::compute_real_cap()
{
	// Initialize real_cap_in_neighborhhod with macines capacities
	for (ul m = 0; m < inst.num_machines; ++m) {
		const Machine * mach = inst.machines[m];
		ul n = mach->neighborhood;
		real_cap_in_neighborhhod[n][mach] = new ul[inst.num_resources];
		memcpy(real_cap_in_neighborhhod[n][mach], mach->capacities, inst.num_resources * sizeof(ul));
	}

	// Substract initial usage for transient resources
	for (ul p = 0; p < inst.num_processes; ++p) {
		const Process * proc = inst.processes[p];
		const Machine * mach = proc->init_assignment;
		const ul n = mach->neighborhood;

		for (ul r = 0; r < inst.num_resources; ++r) {
			if (inst.is_transient[r])	
				real_cap_in_neighborhhod[n][mach][r] -= proc->resources[r];
		}
	}
}

