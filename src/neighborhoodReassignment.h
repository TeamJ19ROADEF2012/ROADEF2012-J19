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

#ifndef __NEIGHBORHOODREASSIGNMENT_H
#define __NEIGHBORHOODREASSIGNMENT_H

#include "instance.h"
#include "solution.h"

#include <list>
#include <map>


/**
 * Decomposes the original problem into subproblems associated to the
 * neighborhoods (one subproblem per neighborhood) and uses
 * Vbp_heuristic to solve them
 */
class NeighborhoodReassignment {

 public:

	/// Constructor 
	NeighborhoodReassignment(const Instance& inst, const Solution& init_sol);

	/// Destructor
	~NeighborhoodReassignment();

	/// Assign each process to a neighborhood (as in the init sol)
	/// and call vbp_heuristic on each neighborhood
	void run();

	/// The solution obtained by concatenating machine assignment
	/// in each neighborhood
	Solution * sol;

 private:

	const Instance& inst;
	const Solution& init_sol;
	const ul nb_neighborhoods;

    /// To manage spread constraints : location of all processes from a service
    std::map< ul, std::list<const Process*> > * services_locations;

	/// for m in neighborhood n, real_capacity[n][m][r] = m->capacity[r] - "m->transient_usage[r]"
	std::map<const Machine *, ul *> * real_cap_in_neighborhhod;

	/// proc_in_neighborhood[n] = vector of the processes 
	/// that we have to assign to the machines of n
	std::vector<const Process*> * proc_in_neighborhood;

	/// Compute the vector real_cap_in_neighborhhod
	void compute_real_cap();

    /// Runs vector bin packing
    void run_vbp(ul machines_neighborhood, ul processes_neighborhood);
};

#endif
