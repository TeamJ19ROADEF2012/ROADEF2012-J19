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

#ifndef __VBP_HEURISTIC_H
#define __VBP_HEURISTIC_H

#include "instance.h"
#include <set>

/** Implementation of two Vector Bin Packing heuristics : First Fit
 * Decreasing (Bin Centering) and Next Fit Decreasing (Bin Balancing).
 */

class Vbp_heuristic {

 public:

	/// Constructor with a vector of machines and a vector of processes to be assigned to
	Vbp_heuristic(const Instance& inst, const std::vector<const Machine*> & neighborhood, 
				  const std::vector<const Process*> & processes,
				  const std::map<const Machine *, ul *> & real_remaining_cap);
	
	/// Destructor
	~Vbp_heuristic();

	/// Run different vbp heuristics
	void run();

	/// For p in "processes", assignment[p] = a machine of "neighborhood"
	std::map<const Process*, const Machine*> assignment;

 private:
	
	const Instance& inst;
	const std::vector<const Machine*> & neighborhood;
	const std::vector<const Process*> & processes;
	const ul nb_proc;
	const ul nb_mach;

	/// For m machine pointer, r resource, real_remaining_cap[m][r] = m->capacity[r] - transient_usage(m, r)"
	const std::map<const Machine *, ul *> & real_remaining_cap;

	/// For r resource, resource_priority[r] = sum_m C(m, r)
	ul * total_cap;

	/// For r resource, resource_priority[r] = sum_p R(p, r)
	ul * total_req;

	/// Structure containing the current machine data
	typedef struct _Mach {
		const Machine * m;                 // Machine pointer
		ul * rem_cap;                      // remaining capacity for each resource
		std::set<const Service*> services; // Services currently assigned to this machine

	    _Mach(const Machine * m, ul * rem_cap):m(m), rem_cap(rem_cap) {};
	} Mach;

	/// Sorted machines in incresing order regarding their initial volume
	typedef std::multimap<ul, Mach> RemainingMachines;
	RemainingMachines rem_mach;

	/// Sorted processes in incresing order regarding their volume
	typedef std::multimap<ul, const Process*> RemainingProcesses;
	RemainingProcesses rem_proc;

	/// Computes the vectors total_cap and total_req
	void compute_resources_priorities();

	/// Computes the rem_proc 
	void compute_proc_volumes();

	/// Computes the rem_mach
	void compute_mach_volumes();

	/// Run the bin centering heuristic, returns true iff it computed a feasible assignment
	bool run_bin_centering();

	/// Run the bin balancing heuristic, returns true iff it computed a feasible assignment
	bool run_bin_balancing();

	/// Reinitialize rem_proc, with the order of the volumes
	void reinit_rem_proc();
	
	/// Reinitialize rem_mach, with the order of the volumes
	void reinit_rem_mach();

	/// Reinitialize rem_mach, with random order
	void reinit_rem_mach(ul * shuffled_id_machines);

};

#endif
