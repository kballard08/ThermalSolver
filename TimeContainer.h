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
#include <sstream>

namespace FEASolverNS
{

using namespace dealii;

/**
 * This class contains information about the time for transient problems. It was discovered after the creation of this class
 * that deal II has other methods for handling time dependant problems, but since this was already implemented and working
 * it was kept.  It also allows the opportunity for further extension if needed, such as variable time steps dependant on
 * error, etc.
 */
class TimeContainer
{
public:
	TimeContainer() : curr_index(0) {};

	/// Intializes the TimeContainer
	/// For now there is just a instialize method for a linear spacing of time steps
	void initialize(double start_t, double end_t, int n);

	/// Increments time and returns false if the there are no more increments to go to
	/// In other words, if the increment method returns false, the analysis is over
	bool increment_time();
	/// Returns the current time value
	double get_current_time() { return time_steps[curr_index]; };
	/// Returns the previous time value
	double get_prev_time() { Assert(curr_index > 0, ExcMessage("TimeContainer method get_prev_time cannot be called when the current time step is 0")); return time_steps[curr_index-1]; };
	/// Returns the time step index (0 based)
	int index() const { return curr_index; };
	/// Returns the time step index (0 based) as a string for use in file names for output
	std::string index_str() const;

private:
	/// Vector of arbitrary times, steps do not have to be uniform
	std::vector<double> time_steps;
	/// current time step index
	unsigned int curr_index;
};

void TimeContainer::initialize(double start_t, double end_t, int n)
{
	Assert(start_t <= end_t, ExcMessage("The start time must be less than or equal to the end time"))
	Assert(n >=0 , ExcMessage("The start time must be less than the end time"))
	// Create the time steps vector
	time_steps.resize(n);
	double delta_time = (end_t - start_t)/((double) n);
	for (int i = 0; i < n; i++) {
		time_steps[i] = start_t + delta_time*i;
	}
}

bool TimeContainer::increment_time() {
	if (curr_index+1 >= time_steps.size())
		return false;
	// If the index increment would not exceed the time steps
	curr_index++;
	return true;
}

std::string TimeContainer::index_str() const {
	std::stringstream ss;
	ss << curr_index;
	return ss.str();
}

}


#endif /* TIMECONTAINER_H_ */
