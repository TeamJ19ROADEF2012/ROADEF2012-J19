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

#ifndef __LOCALIPREASSIGNMENT_H
#define __LOCALIPREASSIGNMENT_H

#include "instance.h"
#include "solution.h"

#include <list>
#include <map>

// Maximum number of processes to be considered at once
#define NTHRESHOLD 60
// Maximum number of machines to be considered at once
#define MTHRESHOLD 2

/** 
 * Structure to help manage spread constraints locally :
 * For a service, contains the required spread on a given subset of machines 
 */
typedef struct _spread_constraint {
    /// Service
    const Service* svc;
    
    /// Required spread (local spreadmin)
    const ul required;

    /// Machines
    std::set<ul> locations;

    _spread_constraint(const Service* s, ul requirement):
        svc(s), required(requirement) { };

    _spread_constraint(const _spread_constraint& s):
        svc(s.svc), required(s.required), locations(s.locations) { };

} Spread_constraint;


/**
 * Given an instance and a current solution, it decomposes the problem
 * into subproblems associated to subsets of machines, and solve them
 * using LocalIP heuristic. The subsets sizes are defined by two
 * constants.
 */
class LocalIPReassignment {

 public:

	/// Constructor 
	LocalIPReassignment(const Instance& inst, const Solution& init_sol);

	/// Destructor
	~LocalIPReassignment();

	/// Assign each process to a neighborhood (as in the init sol)
	/// and call ip heuristic on each neighborhood
	void run();

	/// The solution obtained by concatenating machine assignment
	/// in each neighborhood
	Solution * sol;

    /// Optimizes the assignment on two machines - does not check problem size 
    void run_ip(const Machine* const m1, const Machine* const m2);

    /// Optimizes the assignment on a series of machines - does not check problem size
    void run_ip(const std::vector<const Machine*>& mch);

 private:

    /// The instance
	const Instance& inst;

    /** The solution before applying localIP.
     * Does not have to be the initial solution ! */
	const Solution& init_sol;

    /// Number of neighborhoods
	const ul nb_neighborhoods;

    // To manage spread constraints : location of all processes from a service
    std::map< ul, std::list<const Process*> > * services_locations;

	// for m in neighborhood n, real_capacity[n][m][r] = m->capacity[r] - "m->transient_usage[r]"
	std::map<const Machine *, ul *> * real_cap_in_neighborhhod;

	/// proc_in_neighborhood[n] = vector of the processes 
	/// that we have to assign to the machines of n
	std::vector<const Process*> * proc_in_neighborhood;

	// Compute the vector real_cap_in_neighborhhod
	void compute_real_cap();

    /// Runs integer programming on a neighborhood
    void run_ip(ul neighborhood);

    /// Runs integer programming on the given subneighborhood
    void run_ip(const std::vector<const Process*>& processes,
            const std::vector<const Machine*>& Machine);

    /// Divides processes on a neighborhood depending on their machines
    std::map< const Machine*, std::vector<const Process*> >* divide(ul neighborhood);

    /** Computes spreads for init_sol */
    void spread();

    /** Returns sets to be considered for spread constraints in the subproblem */
    std::list<Spread_constraint>* spread_constraints(const std::vector<const Process*>& processes,
        const std::vector<const Machine*>& machines);

    /// Assign the given processes to the given machine
    void set_sol(const Machine* const mc, const std::vector<const Process*>& proc);

    /// Assign the given process to the given machine
    void set_sol(const Machine* const mc, const Process* const proc);

    /// Tests if the instance is trivial. If so sets the solution and returns true. Otherwise, returns false.
    bool trivial_assignment(const std::vector<const Machine*>& mc, const std::vector<const Process*>& proc);

    /// Sort machines by their load costs
    std::map< ul, const Machine*>* sort_by_load_costs(
        const std::map< const Machine*,std::vector<const Process*> >& proc_on_machines);


#ifdef DEBUG_FUNC
    // DEBUG procedures
    void print_processes();
#endif

};

#endif
