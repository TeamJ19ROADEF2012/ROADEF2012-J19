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

#include <cstdlib>
#include <iostream>

#include "instance.h"
#include "solution.h"
#include "neighborhoodReassignment.h"

using namespace std;

void usage() {
    cout << "Usage : ./test_neighborhoodReassignment [[-t time_limit] -p instance_path ";
    cout << "-i original_solution_path -o output_file -s seed] [-name]" << endl;
}


int main(int argc, char** argv) {
    int time_limit = 300; // time limit in seconds
    char* instance_filename = NULL;
    char* original_solution = NULL;
    char* output = NULL;
    int seed = time((time_t*)NULL);

    // Juste pour éviter le warning
    --time_limit; 

    if (argc == 1) {
        usage();
        exit(0);
    }

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-t")) time_limit = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-p")) instance_filename = argv[++i];
        else if (!strcmp(argv[i], "-i")) original_solution = argv[++i];
        else if (!strcmp(argv[i], "-o")) output = argv[++i];
        else if (!strcmp(argv[i], "-s")) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-name")) {
            cout << "Team J19: Michaël Gabay - Sofia Zaourar" << endl;
            if (argc == 2) exit (0);
        }
    }

    if (instance_filename == NULL) {
        cout << "Please specify an instance" << endl;
        usage();
        exit(1);
    } else if (original_solution == NULL) {
        cout << "Please specify an initial solution" << endl;
        usage();
        exit(1);
    } if (output == NULL) {
        cout << "Please specify the output path" << endl;
        usage();
        exit(1);
    }

    Instance * inst = parse_instance(instance_filename);
    Solution * init_sol = new Solution(original_solution, *inst);
    inst->set_init_assignment(*init_sol);
    init_sol->compute_cost(NULL);
	
	// inst->print_properties();	

	NeighborhoodReassignment * nr_solver = new NeighborhoodReassignment(*inst, *init_sol);
	nr_solver->run();
	
	//nr_solver->sol->print(output); //A decomm si la sol est realisable

	delete nr_solver;
	delete init_sol;
	delete inst;

    return 0;
}
