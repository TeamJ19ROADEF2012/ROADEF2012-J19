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

#include <cassert>
#include "check_constraints.h"

#include <algorithm>

using namespace std;

// does not work for partial solutions
// Requires : the current solution is feasible at this point
bool check_spread_constraints(const Instance& inst, const Process& p,
        const Machine& m, const Solution& sol)
{
	ul spread_min = p.service->spreadmin;
	
    // defensive
    // because of the requirement, if it is <= then it is ==
	if (sol.locations_of_service[p.service->id].size() > spread_min ||
		sol.locations_of_service[p.service->id].find(m.location) == 
		sol.locations_of_service[p.service->id].end())
		return true;

	const Machine& m_prev = *sol.assignment[p.id];
	assert(sol.locations_of_service[p.service->id].find(m_prev.location) != 
		   sol.locations_of_service[p.service->id].end());
	return (sol.locations_of_service[p.service->id][m_prev.location] > 1);
}


bool check_dependency_constraints(const Instance& inst, const Process& p,
								  const Machine& m, const Solution& sol, bool swap) {

	const Service & s = *p.service;
	const Machine & m_prev = *sol.assignment[p.id];

    if (m.neighborhood == m_prev.neighborhood)
        return true;

    // Are the dependencies related to the process satisfied ?
    map<ul, ul>& services_new_neighborhood = sol.services_of_neighborhood[m.neighborhood];
    vector<const Service*>::const_iterator it, end = s.dependencies.end();

    for (it = s.dependencies.begin(); it != end; ++it) {
        // Yes - is it satisfied ? - in the new neighborhood
        //if (services_new_neighborhood.find(*it) == services_new_neighborhood.end())
		if (services_new_neighborhood[(*it)->id] == 0)
            return false;
    }

    // Is some other service depending on this one ?
    map<ul, ul>& services_old_neighborhood =  sol.services_of_neighborhood[m_prev.neighborhood];
    
	// if (services_old_neighborhood.find(p.service)->second >= 2) return true;
    if (swap) {
        if (services_old_neighborhood[s.id] >= 1) return true;
    } else {
        if (services_old_neighborhood[s.id] >= 2) return true;
    }

    end = s.rev_dependencies.end();
    for (it = s.rev_dependencies.begin(); it != end; ++it) {
        // Yes - are there such processes in the previous neighboorhood ?
        // if (services_old_neighborhood.find(*it) != services_old_neighborhood.end()) {
		if (services_old_neighborhood[(*it)->id] != 0){
            return false;
        }
    }

    return true;
}
