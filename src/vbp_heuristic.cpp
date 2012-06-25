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

#include "vbp_heuristic.h"
#include "util.h"

#include <unistd.h>
#include <assert.h>
#include <iostream>
#include <iterator>

using namespace std;

Vbp_heuristic::Vbp_heuristic(const Instance& inst, 
							 const vector<const Machine*> & neighborhood, 
							 const vector<const Process*> & processes,
							 const map<const Machine *, ul *> & real_remaining_cap):
	inst(inst), 
	neighborhood(neighborhood),
	processes(processes),
	nb_proc(processes.size()), nb_mach(neighborhood.size()), 
	real_remaining_cap(real_remaining_cap),
	total_cap(new ul[inst.num_resources]),
	total_req(new ul[inst.num_resources])
{ }


Vbp_heuristic::~Vbp_heuristic() {

	for (RemainingMachines::iterator it = rem_mach.begin(); it != rem_mach.end(); ++it)
		delete [] it->second.rem_cap;

	delete [] total_req;
	delete [] total_cap;
}


void Vbp_heuristic::run() {

	compute_resources_priorities();
	compute_proc_volumes(); 
	compute_mach_volumes();

	// First run the 2 variants on the machine ordred in increasing volume order
	if (run_bin_centering()) {
			return;
	}
	// Reinitialize rem_mach and rem_proc
	reinit_rem_proc(); reinit_rem_mach();

	if (run_bin_balancing()) {
		return;
	}
	
	// Random order on the machines
	ul * shuffled_id_machines = new ul[nb_mach];
	generator& gen = create_rand_generator(0);

	for (ul t = 0; t < 10; ++t) {

		// Reinitialize rem_mach (random order), rem_proc
		reinit_rem_proc();
		generate_shuffled(shuffled_id_machines, nb_mach, gen);
		reinit_rem_mach(shuffled_id_machines);

		// First run bin centering
		if (run_bin_centering()) {
			break;
		}

		// Reinitialize rem_mach (random order), rem_proc
		reinit_rem_proc();
		generate_shuffled(shuffled_id_machines, nb_mach, gen);
		reinit_rem_mach(shuffled_id_machines);

		// Bin balancing
		if (run_bin_balancing())
			break;
	}

	delete [] shuffled_id_machines;
	delete &gen;
}

bool Vbp_heuristic::run_bin_centering() {
	
	RemainingMachines::iterator itm;
	RemainingProcesses::reverse_iterator itp;
	ul * updated_rem_cap = new ul[inst.num_resources];

	// For each machine of the neighborhood (in volume increasing order)
	for (itm = rem_mach.begin(); itm != rem_mach.end(); ++itm) {

		// Consider the machine mach
		Mach * mach = &itm->second;
		//cout << "Machine " << mach->m->id << endl;
		memcpy(updated_rem_cap, mach->rem_cap, inst.num_resources * sizeof(ul));

		// For each process not assignened yet (in volume decreasing order)
		itp = rem_proc.rbegin();
		while (itp != rem_proc.rend()) {
			
			const Process * p = itp->second;

			// Check conflict constraint
			if (!p->service->no_conflict && 
				mach->services.find(p->service) != mach->services.end()) {
				++itp;
				continue;
			}

			// Check capacity constraints 
			bool feasible = true;
			for (ul r = 0; r < inst.num_resources; ++r) {
				if (inst.is_transient[r] && p->init_assignment == mach->m) continue;
				if (p->resources[r] > mach->rem_cap[r]) {
					feasible = false;
					break;
				} 
				updated_rem_cap[r] -= p->resources[r];		
			}

			if (feasible) {
				memcpy(mach->rem_cap, updated_rem_cap, inst.num_resources * sizeof(ul));
				if (!p->service->no_conflict) mach->services.insert(p->service);
				rem_proc.erase(--(itp.base())); // erase itp and pass to the next proc
				assignment[p] = mach->m;
			} else {
				memcpy(updated_rem_cap, mach->rem_cap, inst.num_resources * sizeof(ul));
				++itp;
			}
		}
	}
	
	if (rem_proc.size() != 0) 
		cout << "Infeasible. Nb remaining proc: "  << rem_proc.size() 
			 << "/" << nb_proc << endl;
	else 
		cout << "Feasible !!" << endl;

	delete [] updated_rem_cap;

	return (rem_proc.size() == 0);
}

bool Vbp_heuristic::run_bin_balancing() {

	RemainingProcesses::reverse_iterator itp = rem_proc.rbegin();
	RemainingMachines::iterator itm = rem_mach.begin();
	Mach * mach = &itm->second;
	ul * updated_rem_cap = new ul[inst.num_resources];
	memcpy(updated_rem_cap, mach->rem_cap, inst.num_resources * sizeof(ul));

	// For each process
	while(itp != rem_proc.rend()){
		
		const Process * p = itp->second;
		bool feasible = true;
		ul nb_try = 1;

		// Search the first machine where p fits starting with the current machine
		do {

			feasible = true;
		
			// Conflict
			if (!p->service->no_conflict && 
				mach->services.find(p->service) != mach->services.end()) {
				// We pass to the next machine
				if (++itm == rem_mach.end()) itm = rem_mach.begin();
				mach = &itm->second;					
				memcpy(updated_rem_cap, mach->rem_cap, inst.num_resources * sizeof(ul));
				nb_try++;
				feasible = false;
				continue;
			}

			// Capacity
			for (ul r = 0; r < inst.num_resources; ++r) {
				if (inst.is_transient[r] && p->init_assignment == mach->m) 
					continue;
				if (p->resources[r] > updated_rem_cap[r]) {
					// We pass to the next machine
					if (++itm == rem_mach.end()) itm = rem_mach.begin(); mach = &itm->second;
					memcpy(updated_rem_cap, mach->rem_cap, inst.num_resources * sizeof(ul));
					feasible = false;
					nb_try++;
					break;
				}
				updated_rem_cap[r] -= p->resources[r];
			}

		} while (!feasible && nb_try <= nb_mach);
		
		if (!feasible) {
			break;
		}

		// Apply the assignment
		memcpy(mach->rem_cap, updated_rem_cap, inst.num_resources * sizeof(ul));
		if (!p->service->no_conflict) mach->services.insert(p->service);
		assignment[p] = mach->m;
		rem_proc.erase(--itp.base()); // erase itp contents, itp points on the next proc
		if (++itm == rem_mach.end()) itm = rem_mach.begin(); mach = &itm->second;
		memcpy(updated_rem_cap, mach->rem_cap, inst.num_resources * sizeof(ul));
	}
	delete [] updated_rem_cap;

	// Have we assigned all the processes ?
	if (rem_proc.size() != 0)
		cout << "Infeasible. Nb remaining proc: "  << rem_proc.size() << "/" << nb_proc << endl;
	else 
		cout << "Feasible !!" << endl;
	return (rem_proc.size() == 0);
}




void Vbp_heuristic::compute_resources_priorities()
{
	// Compute total requirement of each resource
	memset(total_req, 0, inst.num_resources * sizeof(ul));
	vector<const Process*>::const_iterator p;
	for (p = processes.begin(); p != processes.end(); ++p) {
		for (ul r = 0; r < inst.num_resources; ++r) {
			total_req[r] += (*p)->resources[r];
		}
	}

	// Compute total capacity for each resource
	memset(total_cap, 0, inst.num_resources * sizeof(ul));
	vector<const Machine*>::const_iterator m;
	for (m = neighborhood.begin(); m != neighborhood.end(); ++m) {
		for (ul r = 0; r < inst.num_resources; ++r) {
			total_cap[r] += (*m)->capacities[r];
		}
	}
}

void Vbp_heuristic::compute_proc_volumes()
{
	// Compute the volume of each process and insert it to rem_proc
	vector<const Process*>::const_iterator p;
	for (p = processes.begin(); p != processes.end(); ++p) {
		ul volume = 0;
		for (ul r = 0; r < inst.num_resources; ++r) {
			volume += (10000 * nb_proc * (*p)->resources[r])/total_req[r];
		}
		rem_proc.insert( pair<ul, const Process*>(volume, *p) );
	}
}

void Vbp_heuristic::compute_mach_volumes()
{
	// Compute the volume of each machine and insert it to rem_mach
	vector<const Machine*>::const_iterator m;
	for (m = neighborhood.begin(); m != neighborhood.end(); ++m) {
		ul volume = 0;
		for (ul r = 0; r < inst.num_resources; ++r) {
			volume += (10000 * nb_mach * (*m)->capacities[r]) / total_cap[r];
		}
		Mach mach(*m, new ul[inst.num_resources]);
		memcpy(mach.rem_cap, real_remaining_cap.find(*m)->second, inst.num_resources * sizeof(ul));
		rem_mach.insert( pair<ul, Mach>(volume, mach) );
	}
}


void Vbp_heuristic::reinit_rem_mach() {
	RemainingMachines::iterator itm;
	for (itm = rem_mach.begin(); itm != rem_mach.end(); ++itm) {
		itm->second.services.clear();
		memcpy(itm->second.rem_cap, real_remaining_cap.find(itm->second.m)->second, inst.num_resources * sizeof(ul));
	}

}

void Vbp_heuristic::reinit_rem_mach(ul * shuffled_id_machines) {

	// Clear and free memory
	RemainingMachines::iterator itm;
	for (itm = rem_mach.begin(); itm != rem_mach.end(); ++itm) {
		delete [] itm->second.rem_cap;
	}
	rem_mach.clear();

	// Create rem_mach with a random order
	for (ul id = 0; id < nb_mach; ++id) {
		const Machine * m = neighborhood[shuffled_id_machines[id]];
		Mach mach(m, new ul[inst.num_resources]);
		memcpy(mach.rem_cap, real_remaining_cap.find(m)->second, inst.num_resources * sizeof(ul));
		rem_mach.insert( pair<ul, Mach>(0, mach) );
	}
}

void Vbp_heuristic::reinit_rem_proc() {
	rem_proc.clear(); 
	compute_proc_volumes();
}
