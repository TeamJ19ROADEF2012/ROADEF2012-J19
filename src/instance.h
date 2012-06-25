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
 
#ifndef __INSTANCE_H
#define __INSTANCE_H


#include <cstring>
#include <vector>
#include <map>

// A useful debug macro
#ifdef DEBUG_FUNC
#define PRINT(a) cout << a << endl
#else
#define PRINT(a) 
#endif

typedef unsigned long ul;

class Solution;
class Instance;

/// Factory / Parsing function
Instance* parse_instance(char* file);

typedef struct _service Service;
typedef struct _machine Machine;
typedef struct _process Process;
typedef struct _balanced_cost Balanced_cost;


/** Machine characteristics */
struct _machine {
    /// Machine's id
    const ul id;

    /// Machine's neighborhood
    const ul neighborhood;

    /// Machine's location
    const ul location;

    /// Machine's capacities for each resource
    ul * const capacities;

    /// Machine's safety capacities for each resource
    ul * const sf_capacities;

    /// Weighted machine move costs
    ul * const mv_costs;

    _machine(ul id, ul num_res, ul num_mc,
            ul ngb, ul loc):
        id(id), neighborhood(ngb), location(loc),
        capacities(new ul[num_res]), sf_capacities(new ul[num_res]),
        mv_costs(new ul[num_mc]) {}

    ~_machine() {
        delete[] capacities;
        delete[] sf_capacities;
        delete[] mv_costs;
    }
};


/** Process characteristics */
struct _process {
    /// Process id
    const ul id;
	
    /// Initial assignment
    Machine const* init_assignment;

    /// Belongs to service... 
    const Service* service;
    
    /// Resources requirements
    ul * const resources;
    
    /// Weighted process move cost
    const ul mv_cost;

    /// Constructor
    _process(ul id, Service* svc, ul num_res, ul num_mc,
			 ul * rsc, ul cost):
	id(id), service(svc), resources(new ul[num_res]),
        mv_cost(cost)
    {
        memcpy(resources, rsc, num_res * sizeof(ul));
    }

    /// Destructor
    ~_process() {
        delete[] resources;
    }
};


/** Service characteristics */
struct _service {

	/// Service id
	const ul id;

	/// True if the service contains only one process
	bool no_conflict;

	/// Spread min (for spread constraints)
	ul spreadmin;

	/// Processes ids that belong to service
	std::vector<const Process*> processes;

	/// Les services dont dépend ce service
	/// Services which the service depends on 
	std::vector<const Service*> dependencies;

	/// Les services qui dépendent de ce service
	/// Services that depend on the given service
	std::vector<const Service*> rev_dependencies;

    _service(ul id):id(id){}
};


/** Balanced cost characteristics */
struct _balanced_cost {
    /// Triple <resource0, resource1, target>
    const ul id;

    /// resource 0
    const ul r0;

    /// resource 1
    const ul r1;

    /// Target
    const ul target;

    /// Balance cost weight
    const ul weight;

    _balanced_cost(ul id, ul rs0, ul rs1,
				   ul target, ul weight) : id(id), r0(rs0), r1(rs1),
        target(target), weight(weight) { }
};


/**
 * Stores all instance data of the considered problem
 */
class Instance {
public:
    /// Number of resources
    const ul num_resources;

    /// Number of machines
    const ul num_machines;

    /// Number of services
    const ul num_services;

    /// Number of processes
    const ul num_processes;

    /// Number of balanced costs
    const ul num_bcosts;
	
    /// Weight of process move cost
    const ul w_pmv;

    /// Weight of service move cost
    const ul w_svc;

    /// Weight of machine move cost
    const ul w_mmv;

    /// is_transient[i] = 1 <=> resource i is transient
    bool * const is_transient;

    /// Weight_loadCost
    ul * const load_cost;

    /// Machines
    Machine ** machines;

	/// Services 
	Service ** services;

    /// Processes 
    Process ** processes;

    /// Balanced costs
    Balanced_cost ** bcosts;

    /// Neighborhoods : Neighborhoods[n] = all machines in neighborhood n
    std::vector< std::vector<const Machine*> > neighborhoods;

    /// Locations locations[l] = all the machines in location l
    std::vector< std::vector<const Machine*> > locations;


    /// Constructor
    Instance(ul num_res, ul num_mc, ul num_svc,
            ul num_proc, ul num_bc, ul pmv,
            ul svc, ul mmv):
        num_resources(num_res), num_machines(num_mc), num_services(num_svc),
        num_processes(num_proc), num_bcosts(num_bc),
        w_pmv(pmv), w_svc(svc), w_mmv(mmv),
        is_transient(new bool[num_res]), load_cost(new ul[num_res]),
        machines(new Machine*[num_mc]),
        services(new Service*[num_svc]),
        processes(new Process*[num_proc]),
        bcosts(new Balanced_cost*[num_bc])
    {
    }

    /// Destructor
    virtual ~Instance();

    /// Returns move cost from machine i to j
    inline ul get_cost(ul i, ul j) const {
        return machines[i]->mv_costs[j];
    }

    /// Returns move cost from machine m1 to m2
    inline ul get_cost(const Machine* const m1, const Machine* const m2) const {
        //return get_cost(m1->id, m2->id);
        return m1->mv_costs[m2->id];
    }

    /** Returns move cost for process p to machine m
     * (which is 0 if m is the initial machine) */
    inline ul process_mvcost(const Process& p, const Machine& m) const {
        if (m.id == p.init_assignment->id) return 0;
        return p.mv_cost + get_cost(p.init_assignment->id, m.id);
    }

    /** Returns move cost for process p to machine m
     * (which is 0 if m is the initial machine) */
    inline ul process_mvcost(const Process& p, ul m) const {
        if (m == p.init_assignment->id) return 0;
        return p.mv_cost + get_cost(p.init_assignment, machines[m]);
    }

    /// Affect the original assignment to processes
    void set_init_assignment (Solution& init_sol);

    /// Parse a file and returns the corresponding Instance
    friend Instance* parse_instance(char* file);

	#ifdef DEBUG
    /// Prints the model and checks whether it is coherent
    void print_model();

	/// Prints some properties of the instance
	void print_properties();
	#endif
};

#endif
