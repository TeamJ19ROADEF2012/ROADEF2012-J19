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

#ifndef __NEXTMACHINE_H
#define __NEXTMACHINE_H

#include "instance.h"

#include <map>

/**
 * Returns the next machine to be considered for integration into
 * the subset of machines for the ILP according to the order
 * defined by a sort.
 */

class NextMachine {
private:
    std::multimap<ul, std::pair<const Machine*, std::vector<const Process*>* > > mch; 

    /** true if the next value to return is the max,
     * false if it is the min */
    bool up;
public:
    NextMachine():up(true) {}

    /** Returns next machine to use and its processes, with at
     * most max_proc processes. The returned pair is (NULL, NULL) if no
     * acceptable machine was found */
    std::pair<const Machine*, std::vector<const Process*>* > next(ul max_proc);

    /// Sort machines by their load costs
    void sort_by_load_costs(const Instance& inst,
            std::map< const Machine*,std::vector<const Process*> >& proc_on_machines);
    
    /// Do not sort machines
    void no_sort(std::map< const Machine*,std::vector<const Process*> >& proc_on_machines);

    /// Randomly sort machines
    void random_sort(long seed,
            std::map< const Machine*,std::vector<const Process*> >& proc_on_machines);

    /** Returns true if it is empty */
    bool empty() {
        return (mch.empty());
    }
};

#endif
