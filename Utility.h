/*
 * Utility.h
 *
 *  Created on: Mar 1, 2013
 *      Author: kballard
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <iostream>
#include <string>

namespace ThermalSolverNS
{

enum Verbosity {
	MIN_V,
	MAX_V,
	DEBUG_V
};

void Status(std::string message, int verbosity, int verbosity_required) {
	if (verbosity >= verbosity_required)
		std::cout << message << std::endl;
}

}

#endif /* UTILITY_H_ */
