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

#ifndef __LIP_H
#define __LIP_H

#include "instance.h"
#include "solution.h"
#include "localIPReassignment.h"

#include <vector>
#include "CoinModel.hpp"
#include "CbcModel.hpp"
#include "OsiSolverInterface.hpp"

Solution* divided_solve(Solution& sol);


/**
 * Generates and solves an ILP model for a restriction of the original
 * problem to a subset of processes and a subset of machines. The
 * model does not include dependency constraints (implicitely checked)
 * and service move costs (because it cannot be computed locally).
 */

class LocalIP {
public:
    /** Constructor
     * Requires : all machines belong to the same neighborhood
     */
    LocalIP(const Instance& inst, const std::vector<const Machine*>& mch,
            const std::vector<const Process*>& proc,
            std::list<Spread_constraint>* spr):
        instance(inst), machines(mch), processes(proc), spread_list(spr), mod(NULL) { }

    /// Destructor
    virtual ~LocalIP();

    /// Solves and returns an assignment
    double solve();

    /// Prints solution attributes
    void print(OsiSolverInterface& si);

    /// Returns the computed assignment
    inline const double* get_assignment() {   
        return mod->solver()->getColSolution();
    }

private:
    // Attributes
    /// The instance
    const Instance& instance;

    /// Machines to consider
    const std::vector<const Machine*>& machines;

    /// Processes to assign
    const std::vector<const Process*>& processes;

    std::list<Spread_constraint>* spread_list;

    /// The model
    CbcModel* mod;

    /// Builds model
    CoinModel* build_model();
};
#endif
