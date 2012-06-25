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
#include <csignal>
#include <sys/time.h>

#include "instance.h"
#include "solution.h"
#include "localopt.h"

#define NUM_THREADS 2


/** \mainpage ROADEF/EURO Challenge 2011, team J19
 *
 * MachineReassignment is the main program, it reads the problem,
 * runs the solving process in parallel on two cores and writes the
 * output when time limit is reached.
 *
 */

using namespace std;

// Global variables for time management, clean end and thread access

/* Program parameters */
int time_limit = 300;           // time limit in seconds, default = 5 mins
char* instance_filename = NULL; // Input file
char* output = NULL;            // Output file
char* original_solution = NULL; // Initial assignment input file
long seed = 0;                  // Seed, default is deterministic
pthread_t tid[NUM_THREADS];     // Threads tids

/* Problem variables */
Instance * inst = NULL;         // The instance
Solution * init_sol = NULL;     // Initial solution
LocalOpt * swap_solver = NULL;
LocalOpt * swap_downward_solver = NULL;


/// Help message
void usage() {
    cout << "Usage : ./machineReassignment [[-t time_limit] -p instance_path ";
    cout << "-i original_solution_path -o output_file -s seed] [-name]" << endl;
}


/// Sets the timer
void set_timer(int time_limit);

/// Parse parameters
void parse_params(int argc, char** argv);

/// Runs the solving process at program level
void solve();

/// Runs the swap solving proces at thread level
void* wrap_swap(void *arg);

/// Runs the idownward search + swap solving proces at thread level
void* wrap_downward_swap(void *arg);

/// Gets the best solution and prints it to the output file
void end_program(int sig);


/// Main program
int main(int argc, char** argv) {

    // Parse parmeters
    parse_params(argc, argv);

    // Sets timer
    set_timer(time_limit);

    // Parse files and creates objects
	inst = parse_instance(instance_filename);
	init_sol = new Solution(original_solution, *inst);

    inst->set_init_assignment(*init_sol);
    init_sol->compute_cost(NULL);

    // Solve
    solve();

    // Ends the program if it finishes before time limit
    // Actually unreachable 
    end_program(0);
}


void solve() {
    // Creates threads
    // This thread uses only swap (preceded by reduce load costs v1)
    swap_solver = new LocalOpt(*inst, *init_sol, *init_sol, seed);
    pthread_create(tid,NULL,wrap_swap, (void*)swap_solver);

    // This thread reduces load costs (v2), then runs downward search and swap
    swap_downward_solver = new LocalOpt(*inst, *init_sol, *init_sol, seed);
    pthread_create(tid+1,NULL,wrap_downward_swap, (void*)swap_downward_solver);

    pthread_join(tid[0],NULL);
	pthread_join(tid[1],NULL);
}


void* wrap_swap(void *arg)
{
    LocalOpt* solver = (LocalOpt*)arg;

    // reduce load_cost version 1
    solver->reduce_lc(true);

    while (1) {
        // Runs swap
		solver->run_swap_downward();
		
        // Runs a downward search to unstuck
	    solver->run_move_proc_downward(false);
    }

    return NULL;
}


void* wrap_downward_swap(void *arg)
{
    LocalOpt* solver = (LocalOpt*)arg;

    // reduce load_cost version 2
    solver->reduce_lc(false);

    // Runs the first downward search
    // Uses random moves
	solver->run_move_proc_downward(false);

    while (1) {
        // Runs swap
		solver->run_swap_downward();

        // Runs the next downward search
        // NO_ALEA is false
	    solver->run_move_proc_downward(false);
    }

    return NULL;
}


void end_program(int sig) {
    Solution * best = NULL;
    ul cost = -1, tmp;

    if (swap_solver != NULL) {
        best = &(swap_solver->solution);
        cost = swap_solver->solution.get_cost();
    }

    if (swap_downward_solver != NULL) {
        tmp = swap_downward_solver->solution.get_cost();
        if (best == NULL || cost > tmp) {
            best = &(swap_downward_solver->solution);
            cost = tmp;
        }
    }

    if (best == NULL) {
        cout << "Error: The program did not finish" << endl;
        exit(2);
    }

    best->print(output);

    // Clean exit, for debug purposes
    // Requires to kill threads before !
    // delete localopt_solver;
    // delete swap_solver;
    // delete swap_downward_solver;

    exit(0);
}


void set_timer(int time_limit) {
    struct itimerval itimer;
    struct timeval tilim;

    if (time_limit < 2) time_limit = 300;

    signal(SIGALRM, end_program);
    gettimeofday(&tilim, NULL);

    // Time limit
    // One second security in the computation,
    // 1.500 seconds security
    time_limit -= 2;
    tilim.tv_sec += time_limit;
    itimer.it_interval.tv_sec = time_limit;
    itimer.it_interval.tv_usec = 500000;
    itimer.it_value.tv_sec = time_limit;
    itimer.it_value.tv_usec = 500000;
    setitimer(ITIMER_REAL, &itimer, NULL);
}


void parse_params(int argc, char** argv) {
    // Missing arguments
    if (argc == 1) {
        usage();
        exit(1);
    }

    // Parse parameters
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-t")) time_limit = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-p")) instance_filename = argv[++i];
        else if (!strcmp(argv[i], "-i")) original_solution = argv[++i];
        else if (!strcmp(argv[i], "-o")) output = argv[++i];
        else if (!strcmp(argv[i], "-s")) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-name")) {
            cout << "J19" << endl;
            if (argc == 2) exit (0);
        }
    }

    // Check input parameters
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
}


