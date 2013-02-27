/*
 * main.cpp
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#include "ThermalProblem.h"

int main ()
{
	try
	{
		using namespace dealii;

		deallog.depth_console (0);

		ThermalProblem<2> thermal_problem;
		thermal_problem.run();
	}
	catch (std::exception &exc)
	{
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;

		std::cerr << "Exception on processing: " << std::endl
				<< exc.what() << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;

		return 1;
	}
	catch (...)
	{
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
	}

	return 0;
}


