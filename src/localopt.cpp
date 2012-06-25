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

#include "localopt.h"
#include "check_constraints.h"

#include <unistd.h>
#include <assert.h>

using namespace std;

LocalOpt::LocalOpt(const Instance &inst, const Solution& init_sol,
				   const Solution& sol, int seed): 
	seed(seed),
	solution(sol), inst(inst), init_sol(init_sol),
	available_resources(new pair<ul, ul>*[inst.num_machines]),
	transient_resources(new ul*[inst.num_machines]),
	balance_costs(new ul*[inst.num_machines]),
	machine_assignment(new set<const Process*>[inst.num_machines]),
	machine_lc(new long[inst.num_machines])
{
    ul res; // temporary variable to store resources

	// Initialization of available_resources, transient_resources,
    // and creation of balance costs array
	Machine ** mc = inst.machines;
	for (ul m = 0; m < inst.num_machines; ++m) {
        balance_costs[m] = new ul[inst.num_bcosts];
		available_resources[m] = new pair<ul, ul>[inst.num_resources];
		transient_resources[m] = new ul[inst.num_resources];
        memset(transient_resources[m], 0, inst.num_resources * sizeof(ul));
		for (ul r = 0; r < inst.num_resources; ++r) {
			available_resources[m][r].first = mc[m]->capacities[r]; 
            available_resources[m][r].second = mc[m]->sf_capacities[r]; 
		}
	}

	// initializing resources
    for (ul p = 0; p < inst.num_processes; ++p) {
        ul m_init = init_sol.assignment[p]->id;
        ul m_curr = sol.assignment[p]->id;
        for (ul r = 0; r < inst.num_resources; ++r) {
			res = inst.processes[p]->resources[r];
			assert(available_resources[m_curr][r].first >= res);
			available_resources[m_curr][r].first -= res; 
			if (available_resources[m_curr][r].second < res)
                available_resources[m_curr][r].second = 0; 
            else
			    available_resources[m_curr][r].second -= res; 
            if (m_init != m_curr && inst.is_transient[r])
                transient_resources[m_init][r] += res;
         }
    }
}


LocalOpt::~LocalOpt(){
    for (ul m = 0; m < inst.num_machines; ++m) {
        delete [] available_resources[m];
        delete [] transient_resources[m];
        delete [] balance_costs[m];
    }
    delete [] available_resources;
    delete [] transient_resources;
    delete [] balance_costs;
	delete [] machine_assignment;
	delete [] machine_lc;
}


void LocalOpt::run_move_proc_downward(bool NO_ALEA){
    ul prev_cost;
    generator* gen;
    if (NO_ALEA)
        gen = NULL;
    else
        gen = &(create_rand_generator(seed));

    // Compute costs, initializes balance_costs and service move costs
    solution.compute_cost(&init_sol, balance_costs);

    do {
        prev_cost = solution.cost;
		assert(solution.cost == solution.compute_cost(&init_sol));
        run_move_proc(gen);
    } while (prev_cost > solution.cost);

    if (gen != NULL)
        delete gen;
}

void LocalOpt::run_swap_downward(){
    ul prev_cost;
	generator& gen = create_rand_generator(seed);
	do {
		prev_cost = solution.cost;
		assert(solution.cost == solution.compute_cost(&init_sol));
        run_swap(gen);
    } while (prev_cost > solution.cost);
    delete &gen;
}


void LocalOpt::run_move_proc(generator* gen){
    ul pid;
    ul * shuffled_proc = NULL;
    
    // Using an randomly sorted array of processes
    if (gen != NULL) {
        shuffled_proc = new ul[inst.num_processes];
        generate_shuffled(shuffled_proc, inst.num_processes, *gen);
    }

	// For each process, try to change its machine with benefit
	for (ul p = 0; p < inst.num_processes; ++p) {
        if (gen == NULL) pid = p;
        else pid = shuffled_proc[p];
        // alternative : pid = (*gen)() * (inst.num_processes);
        optimize_proc(*inst.processes[pid], gen);
    }
}


pair<ul, ul> LocalOpt::move_proc_benefit(const Process& proc,
        const Machine& assigned_machine, ul * assign_balance_costs) {
    ul dlc = 0; // decrease
    ul ibc = 0; // increase
    ul tmp, dec_cost, inc_cost; // temp variables
    const ul assign = assigned_machine.id;
    pair<ul, ul> * const assign_resources = available_resources[assign];
	
    // Computes decrease in load cost if machine is changed (in dlc)
    dlc = 0;
	for (ul j = 0; j < inst.num_resources; ++j) {
        if (assign_resources[j].second != 0) {
            // no load cost for this resource previously 
            // hence no decrease after move
            continue;
        }

        tmp = assigned_machine.capacities[j] - assign_resources[j].first; // Used
        // tmp > assigned_machine.sf_capacities[j] because of the 1st if
        // current load cost (not weighted) for this resource
        dec_cost = tmp - assigned_machine.sf_capacities[j];
        // new resource usage
        tmp -= proc.resources[j];
        // new (not weighted) load cost
        if (tmp > assigned_machine.sf_capacities[j])
            tmp -= assigned_machine.sf_capacities[j];
        else
            tmp = 0;
        dlc += inst.load_cost[j] * (dec_cost - tmp);
	}

    // adds decrease in move cost
    if (!(&assigned_machine == proc.init_assignment)) {
        dlc += inst.process_mvcost(proc,assigned_machine);
    }

    // Computes (in/de)crease in balance costs
    ibc = 0;
    for (ul bb = 0; bb < inst.num_bcosts; ++bb) {
        const Balanced_cost& bc = *(inst.bcosts[bb]);
        inc_cost = bc.target * (assign_resources[bc.r0].first + proc.resources[bc.r0]);
        dec_cost = assign_resources[bc.r1].first + proc.resources[bc.r1];
        if (inc_cost > dec_cost)
            inc_cost -= dec_cost;
        else inc_cost = 0;
        // inc_cost is the new - not weighted - balance cost
        if (assign_balance_costs != NULL)
            assign_balance_costs[bb] = inc_cost;
        if (inc_cost > balance_costs[assign][bb]) {
            // increased balance cost
            ibc += (inc_cost - balance_costs[assign][bb]) * bc.weight;
        } else {
            // decreased balance cost
            dlc += (balance_costs[assign][bb] - inc_cost) * bc.weight;
        }
    }

    return make_pair(ibc, dlc);
}


bool LocalOpt::move_proc(const Process& proc, const Machine& assigned_machine,
        const Machine& current_machine, const std::pair<ul,ul> process_move_change,
        const ul * const new_assign_balance_costs) {
    const ul ibc = process_move_change.first;
    const ul dlc = process_move_change.second;
    // cost = old_cost + inc_cost - dec_cost ... beware ... unsigned
    ul tmp, dec_cost = 0, inc_cost = 0;
    pair<ul, ul> * mresources; // current machine resources
    pair<ul, ul> * const assign_resources = available_resources[assigned_machine.id];
    ul current_balance_costs[inst.num_bcosts];

    mresources = available_resources[current_machine.id];
    // adds process and machine move costs
    inc_cost = (&current_machine == proc.init_assignment) ?
        0 : inst.process_mvcost(proc,current_machine);
    inc_cost += ibc;

    // Eligible ? + compute increase in load cost 
	for (ul j = 0; j < inst.num_resources; ++j) {
        // tmp = resources needed on current machine to assign process
        // taking transience into account
        if (inst.is_transient[j]) {
            tmp = transient_resources[current_machine.id][j];
            if (&current_machine != proc.init_assignment)
                tmp += proc.resources[j];
        } else {
            tmp = proc.resources[j];
        }
        // enough resources ?
		if (mresources[j].first < tmp) {
            // not feasible
            return false;
		}
			
        if (mresources[j].second < proc.resources[j]) {
			inc_cost += (proc.resources[j] - mresources[j].second) *
				inst.load_cost[j];
		}
    }
        
    // Yes => Better cost ?
    dec_cost = dlc;
    // change in service move cost
    // decrease ?
    if (&current_machine == proc.init_assignment &&
		solution.moved_services[proc.service->id] == solution.max_moved_services.first &&
		solution.max_moved_services.second == 1) {
		dec_cost += inst.w_svc;
	}

	// increase ?
    else if ((&assigned_machine == proc.init_assignment) &&
            solution.moved_services[proc.service->id] == solution.max_moved_services.first) {
        // current_machine.id == proc.init_assignment is impossible 
        // because it would be skipped
        inc_cost += inst.w_svc;
    }
        
    // change in balance cost
    for (ul bb = 0; bb < inst.num_bcosts; ++bb) {
        const Balanced_cost& bc = *(inst.bcosts[bb]);
        tmp = bc.target * (mresources[bc.r0].first - proc.resources[bc.r0]);
        const ul tmp1 = mresources[bc.r1].first - proc.resources[bc.r1];
        if (tmp > tmp1)
            tmp -= tmp1;
        else tmp = 0;
        // tmp contains the new - unweighted - balance cost
        current_balance_costs[bb] = tmp;
        if (tmp > balance_costs[current_machine.id][bb]) {
            // increase balance cost
            inc_cost += (tmp - balance_costs[current_machine.id][bb]) * bc.weight;
        } else {
            // decrease balance cost
            dec_cost += (balance_costs[current_machine.id][bb] - tmp) * bc.weight;
        }
    }

    if (inc_cost >= dec_cost) return false; // no benefit
        
    // Yes => check constraints
    // Resources have already been checked (load costs)
    if (!check_conflict_constraints(proc, current_machine, solution) ||
            !check_spread_constraints(inst, proc, current_machine, solution) ||
            !check_dependency_constraints(inst, proc, current_machine, solution)) {
        // some constraint is violated
        return false;
    }
        
    // Constraints verified, assign
    solution.apply_assignment(proc, &current_machine, true);
    solution.cost += inc_cost;
    solution.cost -= dec_cost;
        

	// Update machine_lc (before updating capacities)
	long delta_lc_mas = 0;
	long delta_lc_mcur = 0;
	for (ul r = 0; r < inst.num_resources; ++r) {
		
		// For assigned_machine
		assert(assigned_machine.capacities[r] >= assign_resources[r].first);
		ul old_usage = assigned_machine.capacities[r] - assign_resources[r].first; 
		assert(old_usage >= proc.resources[r]);
		ul new_usage = old_usage - proc.resources[r];
		
		// Si les lc sont des ul
		// if (old_usage > assigned_machine.sf_capacities[r]) {
		// 	delta_lc_mas -= inst.load_cost[r] * (old_usage - assigned_machine.sf_capacities[r]);
		// }
		// if (new_usage > assigned_machine.sf_capacities[r]) {
		// 	delta_lc_mas += inst.load_cost[r] * (new_usage - assigned_machine.sf_capacities[r]);
		// }

		// Si les lc sont des long
		delta_lc_mas -= inst.load_cost[r] * (old_usage - assigned_machine.sf_capacities[r]);
		delta_lc_mas += inst.load_cost[r] * (new_usage - assigned_machine.sf_capacities[r]);

		// For current_machine
		assert(current_machine.capacities[r] >= mresources[r].first);
		old_usage = current_machine.capacities[r] - mresources[r].first;
		new_usage = old_usage + proc.resources[r];
		
		// Si les lc sont des ul
		// if (old_usage > current_machine.sf_capacities[r]) {
		// 	delta_lc_mcur -= inst.load_cost[r] * (old_usage - current_machine.sf_capacities[r]);
		// }
		// if (new_usage > current_machine.sf_capacities[r]) {
		// 	delta_lc_mcur += inst.load_cost[r] * (new_usage - current_machine.sf_capacities[r]);
		// }

		// Si les lc sont des long
		delta_lc_mcur -= inst.load_cost[r] * (old_usage - current_machine.sf_capacities[r]);
		delta_lc_mcur += inst.load_cost[r] * (new_usage - current_machine.sf_capacities[r]);

	}
	machine_lc[assigned_machine.id] += delta_lc_mas; 
	machine_lc[current_machine.id] += delta_lc_mcur;


    // update capacities
	for (ul j = 0; j < inst.num_resources; ++j) {
		mresources[j].first -= proc.resources[j];
        assign_resources[j].first += proc.resources[j];
            
        if (inst.is_transient[j]) {
            if (&current_machine == proc.init_assignment) {
                // assign != proc.init_assignment
                transient_resources[proc.init_assignment->id][j] -= proc.resources[j];
            } else if (&assigned_machine == proc.init_assignment) {
                // current_machine != proc.init_assignment
                transient_resources[proc.init_assignment->id][j] += proc.resources[j];
            }
        }
			
        if (mresources[j].second > proc.resources[j])
            mresources[j].second -= proc.resources[j];
        else
            mresources[j].second = 0;

        tmp = assigned_machine.capacities[j] - assign_resources[j].first; // used
        if (tmp < assigned_machine.sf_capacities[j])
            assign_resources[j].second = assigned_machine.sf_capacities[j] - tmp;
        // if tmp >= sf_cap then it was already = 0 and is still = 0
    }

	// Update balance costs
    memcpy(balance_costs[current_machine.id], current_balance_costs, inst.num_bcosts * sizeof(ul));
    memcpy(balance_costs[assigned_machine.id], new_assign_balance_costs, inst.num_bcosts * sizeof(ul));



    return true;
}


void LocalOpt::optimize_proc(Process& proc, generator* gen) {
    // cost = old_cost + inc_cost - dec_cost ... beware ... unsigned
    const Machine* assigned_machine = solution.assignment[proc.id];
    ul assign_balance_costs[inst.num_bcosts];
    ul shuffled_machines[inst.num_machines];

    if (gen == NULL) {
        for (ul i = 0; i < inst.num_machines; ++i) {
            shuffled_machines[i] = i;
        }
    } else {
        generate_shuffled(shuffled_machines, inst.num_machines, *gen);
    }

    
    pair<ul, ul> p = move_proc_benefit(proc, *assigned_machine, assign_balance_costs);

	for (ul ii = 0; ii < inst.num_machines; ++ii) {
        const Machine* current_machine = inst.machines[shuffled_machines[ii]];
        if (current_machine == assigned_machine) continue;
		
        if (move_proc(proc, *assigned_machine, *current_machine, p, assign_balance_costs)) break;
    }
}


// try to swap processes in the initial solution
void LocalOpt::run_swap(generator& gen){

    solution.compute_cost(&init_sol, balance_costs); // necessaire ???

	// for (ul sw = 0; sw < 5000; sw++) { // VARIANTE 0
	// 	// Choose randomly 2 processes to swap
	// 	ul idp2, idp1 = inst.num_processes * gen();
	// 	do
	// 		idp2 = inst.num_processes * gen();
	// 	while(solution.assignment[idp1] == solution.assignment[idp2]);

	// Choose a random order of the processes
	ul shuffled_proc[inst.num_processes];
    generate_shuffled(shuffled_proc, inst.num_processes, gen);	

	// For each process idp1
	for (ul i1 = 0; i1 < inst.num_processes; ++i1) { // VARIANTE 1
		// for (ul i1 = 0; i1 < inst.num_processes-1; ++i1) { // VARIANTE 2 
		ul idp1 = shuffled_proc[i1];

		//Try to swap it with other processes (assigned to other machines)
		for (ul i2 = 0; i2 < inst.num_processes; ++i2) {  // Variante 1
			// for (ul i2 = i1+1; i2 < inst.num_processes; ++i2) { // Variante 2
			ul idp2 = shuffled_proc[i2];
			if (solution.assignment[idp1] == solution.assignment[idp2]) continue;
			
			if (swap(inst.processes[idp1], inst.processes[idp2], 
					 solution.assignment[idp1], solution.assignment[idp2])) {
				break; // surtout pour variante 1
			}
		}
	}
}


bool LocalOpt::swap(const Process* p1, const Process* p2, const Machine* m1, const Machine* m2){
	// If the swap if feasible
	if (is_feasible_swap(p1, p2, m1, m2)) {
		
		// and interseting
		long delta_cost = delta_cost_swap(p1, p2, m1, m2);
		if ( delta_cost > 0) {
			
			// apply_it to solution
			apply_swap(p1, p2, m1, m2, delta_cost);
			return true;
		}
	}
	return false;
}

// is the swap of p1 and p2 feasible in solution
bool LocalOpt::is_feasible_swap(const Process* p1, const Process* p2, const Machine* m1, const Machine* m2){

	// Pointers to important data
	pair<ul, ul> * resources1 = available_resources[m1->id];
	pair<ul, ul> * resources2 = available_resources[m2->id];

	ul * tresources1 = transient_resources[m1->id];
	ul * tresources2 = transient_resources[m2->id];

	const Service * s1 = p1->service;
	const Service * s2 = p2->service;

	// Conflict constraints
	assert((s1 == s2) == (s1->id == s2->id));
	if ((s1 != s2) && 
		((solution.services_of_machine[m1->id].find(s2) != 
		  solution.services_of_machine[m1->id].end()) ||
		 (solution.services_of_machine[m2->id].find(s1) != 
		  solution.services_of_machine[m2->id].end()))) {
		return false;
	}

	// Dependency constraints
    solution.services_of_neighborhood[m1->neighborhood][p1->service->id] -= 1;
    solution.services_of_neighborhood[m2->neighborhood][p2->service->id] -= 1;
	if (m1->neighborhood != m2->neighborhood && s1 != s2) {
		if (!check_dependency_constraints(inst, *p1, *m2, solution, true) ||
			!check_dependency_constraints(inst, *p2, *m1, solution, true)) {
            solution.services_of_neighborhood[m1->neighborhood][p1->service->id] += 1;
            solution.services_of_neighborhood[m2->neighborhood][p2->service->id] += 1;
			return false;
        }
	}
    solution.services_of_neighborhood[m1->neighborhood][p1->service->id] += 1;
    solution.services_of_neighborhood[m2->neighborhood][p2->service->id] += 1;

	// Capacity constraints
	for (ul r = 0; r < inst.num_resources; ++r) {
		if ((resources1[r].first + p1->resources[r] < p2->resources[r]) ||
			(resources2[r].first + p2->resources[r] < p1->resources[r]))
			return false;
	}

	// Spread constraints 
	if (s1 != s2 && m1->location != m2->location) {
		if (!check_spread_constraints(inst, *p1, *m2, solution) ||
			!check_spread_constraints(inst, *p2, *m1, solution))
			return false;
	}

	// Transient usage constraints
	ul c11 = ul((init_sol.assignment[p1->id] == m1));
	ul c12 = ul((init_sol.assignment[p1->id] == m2));
	ul c21 = ul((init_sol.assignment[p2->id] == m1));
	ul c22 = ul((init_sol.assignment[p2->id] == m2));		

	for (ul r = 0; r < inst.num_resources; ++r) {
		if (!inst.is_transient[r])
			continue;

		if ((resources1[r].first + (1 - c11) * p1->resources[r] < (1 - c21) * p2->resources[r] + tresources1[r]) ||
			(resources2[r].first + (1 - c22) * p2->resources[r] < (1 - c12) * p1->resources[r] + tresources2[r]))
			return false;
	}

	return true;
}

// compute delta cost of the swap of p1 and p2 in solution 
// > 0 if the swap is interesting
long LocalOpt::delta_cost_swap(const Process* p1, const Process* p2, const Machine* m1, const Machine* m2){

	long delta_cost = 0;
		
	const pair<ul,ul> * available_res_m1 = available_resources[m1->id];
	const pair<ul,ul> * available_res_m2 = available_resources[m2->id];

	// Load cost
	for (ul r = 0; r < inst.num_resources; ++r) {
		
		// m1
		assert(m1->capacities[r] >= available_res_m1[r].first);
		ul old_usage = m1->capacities[r] - available_res_m1[r].first;
		assert(old_usage + p2->resources[r] >= p1->resources[r]);
		ul new_usage = old_usage - p1->resources[r] + p2->resources[r];
		
		if (old_usage > m1->sf_capacities[r]) {
			delta_cost -= inst.load_cost[r] * (old_usage - m1->sf_capacities[r]);
		}
		if (new_usage > m1->sf_capacities[r]) {
			delta_cost += inst.load_cost[r] * (new_usage - m1->sf_capacities[r]);
		}

		// m2
		assert(m2->capacities[r] >= available_res_m2[r].first);
		old_usage = m2->capacities[r] - available_res_m2[r].first;
		assert(old_usage + p1->resources[r] >= p2->resources[r]);
		new_usage = old_usage - p2->resources[r] + p1->resources[r];
		
		if (old_usage > m2->sf_capacities[r]) {
			delta_cost -= inst.load_cost[r] * (old_usage - m2->sf_capacities[r]);
		}
		if (new_usage > m2->sf_capacities[r]) {
			delta_cost += inst.load_cost[r] * (new_usage - m2->sf_capacities[r]);
		}
	}

	// Balance cost
	for (ul b = 0; b < inst.num_bcosts; ++b) {

		const  Balanced_cost & bc = *inst.bcosts[b];
		
		// For m1
		long old_bc = bc.target * available_res_m1[bc.r0].first - available_res_m1[bc.r1].first;
		long new_bc = old_bc
			+ bc.target * (p1->resources[bc.r0] - p2->resources[bc.r0])
			- p1->resources[bc.r1] + p2->resources[bc.r1];

		if (old_bc > 0) {
			if (new_bc > 0) {
				delta_cost += bc.weight * (new_bc - old_bc);
			} else {
				delta_cost -= bc.weight * old_bc;
			}
		} else {
			if (new_bc > 0) {
				delta_cost += bc.weight * new_bc;
			}
		}
		
		// For m2
		old_bc = bc.target * available_res_m2[bc.r0].first - available_res_m2[bc.r1].first;
		new_bc = old_bc
			+ bc.target * (p2->resources[bc.r0] - p1->resources[bc.r0])
			- p2->resources[bc.r1] + p1->resources[bc.r1];

		if (old_bc > 0) {
			if (new_bc > 0) {
				delta_cost += bc.weight * (new_bc - old_bc);
			} else {
				delta_cost -= bc.weight * old_bc;
			}
		} else {
			if (new_bc > 0) {
				delta_cost += bc.weight * new_bc;
			}
		}
	}

	// Process and machine move cost
	// For p1
	if (m1 == p1->init_assignment) {
		delta_cost += inst.process_mvcost(*p1, *m2);
	} 
	else if (m2 == p1->init_assignment) {
		delta_cost -= inst.process_mvcost(*p1, *m1);
	}
	else {
		delta_cost += inst.get_cost(p1->init_assignment, m2)
			- inst.get_cost(p1->init_assignment, m1);
	}

	// For p2
	if (m2 == p2->init_assignment) {
		delta_cost += inst.process_mvcost(*p2, *m1);
	} 
	else if (m1 == p2->init_assignment) {
		delta_cost -= inst.process_mvcost(*p2, *m2);
	}
	else {
		delta_cost += inst.get_cost(p2->init_assignment, m1)
			- inst.get_cost(p2->init_assignment, m2);
	}

	// Service move cost
	ul max_mvs = solution.max_moved_services.first;
	ul nb_max_mvs = solution.max_moved_services.second;

	if (p1->service != p2->service) {
		// cout << "Differents services" << endl;
		if ((solution.moved_services[p1->service->id] == max_mvs && m1 == p1->init_assignment) ||
			(solution.moved_services[p2->service->id] == max_mvs && m2 == p2->init_assignment)) {
			// cout << "Increase" << endl;
			delta_cost += inst.w_svc;
		} 
		else if ((nb_max_mvs == 1 && 
				  ((((solution.moved_services[p1->service->id]  == max_mvs && m2 == p1->init_assignment) &&
					(solution.moved_services[p2->service->id]+1 != max_mvs || m2 != p2->init_assignment))) ||
				  ((solution.moved_services[p2->service->id]    == max_mvs && m1 == p2->init_assignment) &&
				   (solution.moved_services[p1->service->id] +1 != max_mvs || m1 != p1->init_assignment)))) ||
				 (nb_max_mvs == 2 &&
				  (solution.moved_services[p1->service->id] == max_mvs && m2 == p1->init_assignment &&
				   solution.moved_services[p2->service->id] == max_mvs && m1 == p2->init_assignment))) {
			
			// cout << "Decrease nb_max_mvs = " << nb_max_mvs << endl;
			delta_cost -= inst.w_svc;
		}
	}
	// p1->service == p2->service
	else if (solution.moved_services[p1->service->id] + 1 >= max_mvs) {
		// cout << "Meme service" << endl;
		int nb_moved_proc = 0;
		if (p1->init_assignment == m1)
			nb_moved_proc++;
		else if (p1->init_assignment == m2)
			nb_moved_proc--;
		if (p2->init_assignment == m2)
			nb_moved_proc++;
		else if (p2->init_assignment == m1)
			nb_moved_proc--;
		
		if ((nb_moved_proc >= 0 || (nb_moved_proc == -1 && nb_max_mvs == 1)) &&
			solution.moved_services[p1->service->id] == max_mvs) {
			delta_cost += nb_moved_proc * inst.w_svc;
		}
		else if (nb_moved_proc == 2 && (solution.moved_services[p1->service->id] == max_mvs -1)) {
			delta_cost += inst.w_svc;
		}
		else if (nb_moved_proc == -2 && nb_max_mvs == 1 && solution.moved_services[p1->service->id] == max_mvs) {

			// FIXME : en stockant le nb_moved_service[i] 
			// = nb de services dont i proc sont déplacés, pour tout i
			bool decrease_by_1 = false;
			for (ul s = 0; s < inst.num_services; ++s) {
				if (solution.moved_services[s] == max_mvs -1) {
					decrease_by_1 = true;
					break;
				}
			}
			if (decrease_by_1) {
				delta_cost -= inst.w_svc;
			} else {
				delta_cost -= 2 * inst.w_svc;
			}
		}
	}

	return - delta_cost;
}

// Update solution by swapping p1 and p2
void LocalOpt::apply_swap(const Process* p1, const Process* p2, const Machine* m1, const Machine* m2, ul delta_cost) {

	// Cost
	solution.cost -= delta_cost;

	solution.apply_swap(*p1, m2, *p2, m1);

	// available_resources
	for (ul r = 0; r < inst.num_resources; ++r) {
		// On m1
		assert(available_resources[m1->id][r].first + p1->resources[r] >= p2->resources[r]);
		available_resources[m1->id][r].first += p1->resources[r] - p2->resources[r];
		assert(m1->capacities[r]  >= available_resources[m1->id][r].first);
		ul new_usage = m1->capacities[r] - available_resources[m1->id][r].first;
		if (new_usage >= m1->sf_capacities[r])
			available_resources[m1->id][r].second = 0;
		else
			available_resources[m1->id][r].second = m1->sf_capacities[r] - new_usage;
		// On m2
		assert(available_resources[m2->id][r].first + p2->resources[r] >= p1->resources[r]);
		available_resources[m2->id][r].first += p2->resources[r] - p1->resources[r];
		assert(m2->capacities[r] >= available_resources[m2->id][r].first);
		new_usage = m2->capacities[r] - available_resources[m2->id][r].first;
		if (new_usage >= m2->sf_capacities[r])
			available_resources[m2->id][r].second = 0;
		else
			available_resources[m2->id][r].second = m2->sf_capacities[r] - new_usage;
	}
	// transient_resources
	for (ul r = 0; r < inst.num_resources; ++r) {
		if (!inst.is_transient[r]) continue;
		// On m1
		if (m1 == p1->init_assignment)
			transient_resources[m1->id][r] += p1->resources[r];
		if (m1 == p2->init_assignment) {
			assert(transient_resources[m1->id][r] >= p2->resources[r]);
			transient_resources[m1->id][r] -= p2->resources[r];
		}
		// On m2
		if (m2 == p2->init_assignment)
			transient_resources[m2->id][r] += p2->resources[r];
		if (m2 == p1->init_assignment) {
			assert(transient_resources[m2->id][r] >= p1->resources[r]);
			transient_resources[m2->id][r] -= p1->resources[r];
		}
	} 
}


void LocalOpt::compute_ordered_lc_machines() {

	ordered_lc_machines.clear();

	// For each machine m
	for (ul idm = 0; idm < inst.num_machines; ++idm) {
		const Machine * m = inst.machines[idm];

		// Compute its load cost m_lc
		ul m_lc = 0;
		for (ul r = 0; r < inst.num_resources; ++r) {
			// cap - available_re = usage => condition is if usage > sfcap
			//if (m->capacities[r] > available_resources[idm][r].first + m->sf_capacities[r])
			m_lc += (long)inst.load_cost[r] * 
				((long)m->capacities[r] - (long)available_resources[idm][r].first - (long)m->sf_capacities[r]);
		}
		// Insert the pair in ordered_lc_machines
		ordered_lc_machines.insert( pair<ul, const Machine*>(m_lc, m) );
		machine_lc[idm] = m_lc;
	}
}


ul LocalOpt::compute_lc_lower_bound() {
    ul used_resources[inst.num_resources];  // overall used resources 
    ul sf_res[inst.num_resources];          // overall safety capacities 
    ul load_cost_lb = 0;                    // lower bound on load costs


    memset(used_resources, 0, inst.num_resources * sizeof(ul));
    memset(sf_res, 0, inst.num_resources * sizeof(ul));

	// For each machine m
	for (ul idm = 0; idm < inst.num_machines; ++idm) {
		const Machine * m = inst.machines[idm];

		for (ul r = 0; r < inst.num_resources; ++r) {
            if (m->capacities[r] < available_resources[idm][r].first) PRINT("WTF !!");
            used_resources[r] += m->capacities[r] - available_resources[idm][r].first;
            sf_res[r] += m->sf_capacities[r];
		}
	}

    // Computes the lower bound
	for (ul r = 0; r < inst.num_resources; ++r) {
        if (used_resources[r] < sf_res[r]) continue;
        load_cost_lb += (used_resources[r] - sf_res[r]) * inst.load_cost[r];
    }

    return load_cost_lb;
}


void LocalOpt::init_machine_assignment() {

	for (ul idp = 0; idp < inst.num_processes; ++idp) {
		const Machine * m  = solution.assignment[idp];
		machine_assignment[m->id].insert(inst.processes[idp]);
	}
}

void LocalOpt::reduce_lc(const bool version) {

	solution.compute_cost(&init_sol, balance_costs);

	// initializing machine_assignment
	init_machine_assignment();
	
	// Init machine load costs
	compute_ordered_lc_machines();

	OrderedMachines::reverse_iterator itm_expve;
	OrderedMachines::iterator itm_cheap;
	set<const Process *>::iterator itp_expve;

	const Machine * m_cheap;
	const Machine * m_expve;
	const Process * p_expve;
	
 startAgain:

	// For each machine m_expve with positive load cost
	for (itm_expve = ordered_lc_machines.rbegin(); itm_expve != ordered_lc_machines.rend() && itm_expve->first > 0 ; ++itm_expve) {
		m_expve = itm_expve->second;

		// For each process p_expve assigned to m_expve
		for (itp_expve = machine_assignment[m_expve->id].begin(); itp_expve != machine_assignment[m_expve->id].end(); ) {
			p_expve = *itp_expve;
			ul bc[inst.num_bcosts];
			pair<ul,ul> inc_dec = move_proc_benefit(*p_expve, *m_expve, bc);

			// For each machine m_cheap with negative load cost
			//for (itm_cheap = ordered_lc_machines.begin(); itm_cheap != ordered_lc_machines.end() && itm_cheap->first <= 0; ++itm_cheap) {       // version1
			//for (itm_cheap = ordered_lc_machines.begin(); itm_cheap != ordered_lc_machines.end() && itm_cheap->second != m_expve; ++itm_cheap) {  // version2
            bool keep_going = true;
			for (itm_cheap = ordered_lc_machines.begin(); itm_cheap != ordered_lc_machines.end() && keep_going; ++itm_cheap) {  // version2
				m_cheap = itm_cheap->second;

				// If the move is interesting, apply it
				if (move_proc(*p_expve, *m_expve, *m_cheap, inc_dec, bc)) {

					// Update machine_assignment
					machine_assignment[m_cheap->id].insert(p_expve);
					machine_assignment[m_expve->id].erase(itp_expve++);

					// Update ordered_lc_machine
					ordered_lc_machines.erase(itm_expve->first);
					ordered_lc_machines.insert(pair<long, const Machine*>(machine_lc[m_expve->id], m_expve));
					ordered_lc_machines.erase(itm_cheap);
					ordered_lc_machines.insert(pair<long, const Machine*>(machine_lc[m_cheap->id], m_cheap));

					goto startAgain; // because the most expensive/cheap machine may have changed
				}

                if (version) {
                    // Version 1
			        keep_going = (itm_cheap->first <= 0);
                } else {
                    // Version 2
			        keep_going = (itm_cheap->second != m_expve);
                }
			
			}
			// We failed to assign p_expve to any cheap machine, we try with the next process
            itp_expve++;
		}
	}
}
