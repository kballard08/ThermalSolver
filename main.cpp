/*
 * main.cpp
 *
 *  Created on: Feb 24, 2013
 *      Author: kballard
 */

#include <deal.II/base/exceptions.h>

#include "Executive.h"
#include "ScriptReader.h"

using namespace dealii;

int main (int argc, char** argv)
{
	try
	{
		using namespace dealii;

		// Initialize the MPI world
		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
		{
			deallog.depth_console (0);

			if (argc == 1) {
				// TODO: run a test instead?
				std::cout << "A path to an input script must be taken in as an argument." << std::endl;
				return 0;
			}
			else if (argc == 2) {
				FEASolverNS::ScriptReader* sr = new FEASolverNS::ScriptReader(argv[1]);

				std::vector<std::string> tokens;
				sr->get_next_line(tokens);
				if (tokens[0] == "SetDim") {
					Assert(tokens.size() == 2, dealii::ExcMessage("The command SetDim in the input script has the wrong ammound of arguments."))

					if (tokens[1] == "2") {
						FEASolverNS::Executive<2> ex;
						ex.run(sr); // pass ScriptReader to executive class to finish processing input script
					}
					else if (tokens[1] == "3") {
						FEASolverNS::Executive<3> ex;
						ex.run(sr); // pass ScriptReader to executive class to finish processing input script
					}
					else {
						Assert(false, ExcNotImplemented())
					}

				}
				else {
					std::cout << "The command SetDim must be the first command in the input script." << std::endl;
					return 0;
				}
			}
			else {
				std::cout << "Usage: " << argv[0] << " [optional] script_file_path" << std::endl;
				std::cout << "\tIf no script file path is provided, then the program will run a simple test." << std::endl;
				return 0;
			}
		}
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


