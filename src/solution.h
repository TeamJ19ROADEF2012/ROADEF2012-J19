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

#ifndef __SOLUTION_H
#define __SOLUTION_H

#include <cstdlib>
#include <set>
#include <map>
#include <pthread.h>

#include "instance.h"

/**
 * Stores a solution of the problem - including the assignment and
 * some structures with solution details, along with methods to
 * compute the cost.
 */

class Solution {
private:

    /// Buffer for interruptions and thread safety
    const Solution* tmp_sol;

    /// Mutex to change solution
    pthread_mutex_t plock;
	

public:
    /// Solution cost - may not be updated
    ul cost;
	
	/// Instance
	const Instance& inst;

    /// assignment[id_p] = m <=> process p runs on machine m
	/// Array of pointers that point to constant Machines
    Machine const ** assignment;

	/// services_of_machine[m] = set of the services assigned to machine m
	std::set<const Service *> * services_of_machine;
		
    /** Contains the number of process from services on each location
	 * where the service runs: 
	 * services_of_location[s][l] = nb of proc of s in l (>0)
     * updated in apply_assignment */
	std::map<ul, ul> * locations_of_service;
	
    /** Contains the number of process from services of each neighborhood
	 * services_of_neighborhood[n][s] = nb of proc of s in n
	 * updated in apply_assignment */
    std::map<ul, ul> * services_of_neighborhood;

    /// Moved Services - computed using compute_cost, not copied when =
    ul * moved_services;

    /** Maximum number of moved services - computed using
     * compute_cost, not copied when =
	 * max_moved_services.first = maximal number of processes moved in a service
	 * max_moved_services.sercond = number of services that reach this max
     */
	std::pair<ul, ul> max_moved_services;
	
    /// Constructor
    Solution(const Instance& inst):tmp_sol(NULL), cost(0), inst(inst), 
		assignment(new Machine const *[inst.num_processes]),
		services_of_machine(new std::set<const Service*> [inst.num_machines]),
		locations_of_service(new std::map<ul, ul>[inst.num_services]),
	    services_of_neighborhood(new std::map<ul, ul>[inst.neighborhoods.size()]),
        moved_services(new ul[inst.num_services]) {
		
        pthread_mutex_init(&plock, NULL);
	}
			
	/// Copy constructor
    Solution(const Solution& sol): tmp_sol(NULL), cost(sol.cost), inst(sol.inst),
		assignment(new Machine const *[inst.num_processes]),
		services_of_machine(new std::set<const Service*> [inst.num_machines]),
		locations_of_service(new std::map<ul, ul>[inst.num_services]),
	    services_of_neighborhood(new std::map<ul, ul>[inst.neighborhoods.size()]),
        moved_services(new ul[inst.num_services]) {
				
        pthread_mutex_init(&plock, NULL);
		memcpy(assignment, sol.assignment, inst.num_processes * sizeof(Machine*));

		for (ul m = 0; m < inst.num_machines; ++m) {
			services_of_machine[m] = sol.services_of_machine[m];
		}

		for (ul s = 0; s < inst.num_services; ++s)
			locations_of_service[s] = sol.locations_of_service[s];

		for (ul n = 0; n < inst.neighborhoods.size(); ++n)
			services_of_neighborhood[n] = sol.services_of_neighborhood[n];

		/* for (ul s = 0; s < inst.num_services; ++s) */
		/* 	moved_services[s] = sol.moved_services[s]; */
		/* max_moved_services = sol.max_moved_services; */
    }

    /** Initializes solution structures :
     * services_of_machine, services_of_neighborhood, locations_of_service */
    void init(const Solution* init_sol);

    /// Destructor
    virtual ~Solution();

    /// Operator =
    Solution& operator=(const Solution& sol) {
        tmp_sol = &sol;

        pthread_mutex_lock(&plock);

        // pour se protéger des modifs // ininteressantes.
        // if (cost >= sol.cost) return *this;
        cost = sol.cost;
        memcpy(assignment, sol.assignment, inst.num_processes * sizeof(Machine*));
		
        for (ul m = 0; m < inst.num_machines; ++m) {
            services_of_machine[m] = sol.services_of_machine[m];
		}

		for (ul s = 0; s < inst.num_services; ++s)
			locations_of_service[s] = sol.locations_of_service[s];

		for (ul n = 0; n < inst.neighborhoods.size(); ++n)
			services_of_neighborhood[n] = sol.services_of_neighborhood[n];

		for (ul s = 0; s < inst.num_services; ++s)
			moved_services[s] = sol.moved_services[s];
		max_moved_services = sol.max_moved_services;
		
        pthread_mutex_unlock(&plock);
        tmp_sol = NULL;

		return *this;
    }

    /// Parser: constructs a solution based on a file
    Solution(char* file, const Instance& inst);

	/// To keep a copy of the best solution
	void partial_copy(const Solution & sol) {
        tmp_sol = &sol;
        pthread_mutex_lock(&plock);
		cost = sol.cost;
		memcpy(assignment, sol.assignment, inst.num_processes * sizeof(Machine*));
        pthread_mutex_unlock(&plock);
	}
	
    ul get_cost() {
        ul tmp;

        pthread_mutex_lock(&plock);
        tmp = cost;
        pthread_mutex_unlock(&plock);

        return tmp;
	}


    /** Computes the solution cost
     * Pass NULL as init_sol parameter if it is the initial solution
     * If store_balance_costs != NULL then the balance resources (not weighted)
     * are stored into it (store_balance_costs double index array
     * should have already been allocated)
     */
    ul compute_cost(const Solution * init_sol, ul ** store_balance_resources = NULL);

    /// Output solution to the given file
    void print(const char* file); 

    /** Apply a new assignmet for a given process,
     * it is a new assignment, move has to be false, otherwise - if
     * we are modifying an assignment where p has already been assigned -
     * move should be true.
     * lock should be true iff we need to lock the mutex */
    void apply_assignment(const Process& proc, const Machine* const m,
            bool move = false, bool lock = true);

    /** Apply a swap */
    void apply_swap(const Process& proc1, const Machine* const m1,
                    const Process& proc2, const Machine* const m2);

    /** Returns all processes assigned to the given neighborhood */
    std::vector<const Process*>* neighborhood_processes(ul neighborhood);

private:
    /// Apply common operations needed when moving p on m
    void apply_move_assignment(const Process& proc, const Machine* const m);
};

#endif
