/*
 * TimeContainer.h
 *
 *  Created on: Apr 25, 2013
 *      Author: kballard
 */

#ifndef TIMECONTAINER_H_
#define TIMECONTAINER_H_

#include <deal.II/base/exceptions.h>
#include "Utility.h"
#include <string>

namespace FEASolverNS
{

using namespace dealii;

class TimeContainer
{
public:
	TimeContainer(Verbosity verb) : verbosity(verb), curr_index(0) {};

	// Process the transient commands
	bool process_cmd(std::vector<std::string> tokens);

	// Increments time and returns false if the there are no more increments to go to
	// In other words, if the increment method returns false, the analysis is over
	bool increment_time();
	// Returns the current time value
	double get_current_time();
	// Returns the previous time value
	double get_prev_time();
	// Returns the time step index (0 based) as a string for use in file names for output
	std::string index_str();

private:
	Verbosity verbosity;
	std::vector<double> time_steps;
	int curr_index;
};

bool TimeContainer::process_cmd(std::vector<std::string> tokens)
{
	// TODO have transient analysis in thermalelasticityproblem and have options block marked by end
	if (tokens[0] == "TransientAnalysis") {
		Status("Analysis has been set to transient.", verbosity, MIN_V);
		Assert(tokens.size() == 6, ExcMessage("TransientAnalysis command in the input script did not meet the expected parameters of start_time end_time n_time_steps."))
		// Expected tokens are:
		// start_time end_time n_time_steps

		double start_time = atof(tokens[1].c_str());
		double end_time = atof(tokens[2].c_str());
		int n_time_steps = atoi(tokens[3].c_str());
		time_steps.resize(n_time_steps, 0.0);

		// Create the time steps vector
		double delta_time = (end_time - start_time)/((double) n_time_steps);
		for (int i = 0; i < n_time_steps; i++)
			time_steps[i] = start_time + delta_time*i;

		// Move to problem is_transient = true;
	}
	else {
		return false;
	}

	return true;
}

}


#endif /* TIMECONTAINER_H_ */
