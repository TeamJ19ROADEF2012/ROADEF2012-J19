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

/** Functions used to check problem contraints */

#ifndef __CHECK_CONSTRAINTS_H
#define __CHECK_CONSTRAINTS_H

#include "instance.h"
#include "solution.h"

inline bool check_conflict_constraints(const Process& p,
        const Machine& m, const Solution& sol)
{
	return (sol.services_of_machine[m.id].find(p.service) ==
			sol.services_of_machine[m.id].end());
}

// does not work for partial solutions
// Requires : the current solution is feasible at this point
bool check_spread_constraints(const Instance& inst, const Process& p,
        const Machine& m, const Solution& sol);

// Checks that the dependencies are still satisfied
// Requires : the dependencies were previously satisfied
bool check_dependency_constraints(const Instance& inst, const Process& p,
        const Machine& m, const Solution& sol, bool swap = false);

#endif
