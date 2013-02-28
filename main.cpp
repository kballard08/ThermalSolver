/*
 * main.cpp
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#include "ThermalProblem.h"
#include "ScriptReader.h"

int main (int argc, char* argv[])
{
	// TODO: delete this comment block.  This was simply for testing the ScriptReader
	/*
	ThermalSolverNS::ScriptReader sr("script_test.txt");
	std::vector<std::string> tokens;
	int lineNum = 0;
	while(sr.get_next_line(tokens)) {
		std::cout << "Line: " << lineNum++ << std::endl;
		for (unsigned int i = 0; i < tokens.size(); i++)
			std::cout << tokens[i] << std::endl;
	}
	*/

	try
	{
		using namespace dealii;

		deallog.depth_console (0);

		ThermalSolverNS::ThermalProblem<2> thermal_problem;

		if (argc ==1)
			thermal_problem.run();
		else
			thermal_problem.run(); // change to use script
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


