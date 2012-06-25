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

#include "nextMachine.h"
#include "solution.h"
#include "util.h"

using namespace std;

void NextMachine::sort_by_load_costs(const Instance& inst,
        map< const Machine*,vector<const Process*> >& proc_on_machines) {
    ul res[inst.num_resources];
    long cost;

    // Compute load costs
    map< const Machine*,vector<const Process*> >::iterator mit = proc_on_machines.begin();
    for (; mit != proc_on_machines.end(); ++mit) {
        const Machine * const mc = (*mit).first;
        memset(res, 0, sizeof(ul)*inst.num_resources);
        vector<const Process*>::const_iterator it = (*mit).second.begin();
        for (; it != (*mit).second.end(); ++it) {
            for (ul k = 0; k < inst.num_resources; ++k) {
                res[k] += (*it)->resources[k];
            }
        }
        cost = 0;
        for (ul k = 0; k < inst.num_resources; ++k) {
            if (res[k] > mc->sf_capacities[k]) {
				cost += inst.load_cost[k] * (res[k] - mc->sf_capacities[k]);
            }
        }
        mch.insert(make_pair(cost, make_pair(mc,&((*mit).second))));
    }
}


void NextMachine::no_sort(map< const Machine*,vector<const Process*> >& proc_on_machines) {
    map< const Machine*,vector<const Process*> >::iterator mit = proc_on_machines.begin();

    for (ul i = 0; mit != proc_on_machines.end(); ++i, ++mit) {
        const Machine * const mc = (*mit).first;
        mch.insert(make_pair(i, make_pair(mc,&((*mit).second))));
    }
}


void NextMachine::random_sort(long seed,
        map< const Machine*,vector<const Process*> >& proc_on_machines) {
    map< const Machine*,vector<const Process*> >::iterator mit = proc_on_machines.begin();
    ulgenerator& gen = create_rand_ulgenerator(seed);

    for (ul i = 0; mit != proc_on_machines.end(); ++i, ++mit) {
        const Machine * const mc = (*mit).first;
        mch.insert(make_pair(gen(), make_pair(mc,&((*mit).second))));
    }

    delete &gen;
}


pair < const Machine *, vector < const Process *>*> NextMachine::next(ul max_proc)
{
    pair < const Machine *, vector < const Process *>*>ret(NULL, NULL);

    if (up) {
	    multimap < ul, pair < const Machine *,
	        vector <const Process*>* > >::iterator it, end;
	    it = mch.begin();
	    end = mch.end();
	    for (; it != end; ++it) {
	        if ((*it).second.second->size() > max_proc) continue;
		    up = false;
		    ret = (*it).second;
		    mch.erase(it);
		    return ret;
	    }
    } else {
	    multimap < ul, pair < const Machine *,
	        vector <const Process*> *> >::reverse_iterator it, end;
	    it = mch.rbegin();
	    end = mch.rend();
	    for (; it != end; ++it) {
	        if ((*it).second.second->size() > max_proc) continue;
		    up = true;
		    ret = (*it).second;
		    mch.erase(--(it.base()));
		    return ret;
	    }
    }

    return ret;
}
