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

#ifndef __UTIL_H
#define __UTIL_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

/**
 * Some useful definitions for random sort and random number generation
 */

typedef boost::variate_generator< boost::mt19937, boost::uniform_real<double> > generator;
typedef boost::variate_generator< boost::mt19937, boost::uniform_int<ul> > ulgenerator;

// Shuffle an array of n elements (indexes 0..n-1)
// Knuth shuffle modified to generate directly the array
inline void generate_shuffled(ul * array, ul n, generator& gen) {
    ul i, j;

    array[0] = 0;
    for (i = 1; i < n; ++i) {
        j = gen() * (i + 1);
        array[i] = array[j];
        array[j] = i;
    }
}

inline generator& create_rand_generator(boost::mt19937::result_type seed) {
    boost::mt19937 gener(seed);
    boost::uniform_real<double> uni_dist(0,1);
    generator* uni = new generator(gener, uni_dist);
	return *uni;
}


inline ulgenerator& create_rand_ulgenerator(boost::mt19937::result_type seed) {
    boost::mt19937 gener(seed);
    boost::uniform_int<ul> uni_dist(0,((ul)-1));
    ulgenerator* uni = new ulgenerator(gener, uni_dist);
	return *uni;
}

#endif
