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

#ifndef __LOCALOPT_H
#define __LOCALOPT_H

#include "instance.h"
#include "solution.h"
#include "util.h"

/**
 * A local search heuristic. Two main moves are implemented:
 * moving processes and swaping processes. Some routines
 * that make smart use of those mouvements are also provided.
 */

class LocalOpt {
public:
	int seed;

    /** Initializes the local search heuristic
     * REQUIRES a feasible initial solution
     */
    LocalOpt(const Instance& instance, const Solution& init_sol,
			 const Solution& sol, int seed);
	virtual ~LocalOpt();

    /// Run the moves until the solution cannot be improved
	void run_move_proc_downward(bool NO_ALEA = false);
	void run_swap_downward();

	/// Current solution
	Solution solution;

	/// Initialize machine_assignment using solution - require machine_assignment "empty"
	void init_machine_assignment();

	/** Compute or recompute ordered_lc_machines and machine_lc using solution
     * Can be used to initialize the structure or to update it totally */
	void compute_ordered_lc_machines();

    /** "Smart moves" : moves processes in order to reduce load costs
     * Tries to move processes from machines with large load cost to
     * machines with lower load cost. Accepts the first feasible move.
     * Version 1 assigns a process only to machines with negative load costs.
     * Version 2 assigns a process to all machines with lower load cost than its current machine
     * \param version true for version 1, false for version 2 */
	void reduce_lc(const bool version = true);

    /// Computes lower bound on load costs
    ul compute_lc_lower_bound();

private:
	const Instance& inst;
	const Solution& init_sol;
	
	/// available_resources[m][r] = < Cap(m,r) - U(m,r),  max(0, Sf_Cap(m, r) - U(m, r)) >
	std::pair<ul, ul> ** available_resources;
	ul ** transient_resources;
    ul ** balance_costs;

	/// machine_assignment[idm] = set of processes assigned to idm in solution
	std::set<const Process*> * machine_assignment;

	/// Sturcture to store machine pointers ordered by a key value
	typedef std::multimap<long, const Machine*> OrderedMachines;

	/// Machines ordered in increasing signed load cost in solution order
	/// here lc = U(m, r) - Sf_Cap(m, r)
	OrderedMachines ordered_lc_machines;

	// machine_lc[idm] = signed load cost of the machine in solution
	long * machine_lc;


	// Move processes

    // If gen == NULL, then no alea
	/* void run_opt(generator* gen); */
	void run_move_proc(generator* gen);
    
    /// Moves the given process on the first interesting machine
    void optimize_proc(Process& proc, generator* gen);

    /** Returns increase and decrease in costs if we remove
     * the specified process from assigned_machine. Also computes
     * balanced costs if assign_balance_costs != NULL
     * \return a pair (increase, decrease) */
    std::pair<ul, ul> move_proc_benefit(const Process& proc,
            const Machine& assigned_machine, ul * assign_balance_costs);

    /** Tries to move process proc from assigned_machine to current_machine,
     * if it is feasible and reduces cost, process is moved and true is returned
     * o.w. nothing is done and false is returned
     * REQUIRES : assigned_machine != current_machine */
    bool move_proc(const Process& proc, const Machine& assigned_machine,
            const Machine& current_machine, const std::pair<ul,ul> process_move_change,
            const ul * const new_assign_balance_costs);
     

	// Swap processes

	/// Try to swap diffrent random pairs of processes in solution
	void run_swap(generator& gen);

	/// Try to Swap p1 and p2 in solution, returns true if it succeeds
	bool swap(const Process* p1, const Process* p2, 
			  const Machine* m1, const Machine* m2);

	/// Is the swap of p1 and p2 feasible in solution ?
	bool is_feasible_swap(const Process* p1, const Process* p2, 
						  const Machine* m1, const Machine* m2);

	/// Compute delta cost of the swap of p1 and p2 in solution
	long delta_cost_swap(const Process* p1, const Process* p2,
						 const Machine* m1, const Machine* m2);

	/// Update solution by swapping p1 and p2 -> TODO: utiliser le apply_assignent de solution
	void apply_swap(const Process* p1, const Process* p2,
					const Machine* m1, const Machine* m2,
					ul delta_cost);
};

#endif
