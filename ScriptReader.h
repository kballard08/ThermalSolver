/*
 * ScriptReader.h
 *
 *  Created on: Feb 27, 2013
 *      Author: kballard
 */

#ifndef SCRIPTREADER_H_
#define SCRIPTREADER_H_

#include <deal.II/base/exceptions.h>

#include <iterator>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace FEASolverNS
{

using namespace dealii;

/**
 * This class it meant to read a ascii based script and return each line as a vector of strings.  The strings are split
 * by whitespace in the script.
 */
class ScriptReader
{
public:
	ScriptReader() {};
	ScriptReader(std::string filePath);
	~ScriptReader();

	/// Use if the default constructor was called
	void open(std::string filePath);

	/// Return file indicates if the read was valid,
	/// false will be returned if the eof was reached
	bool get_next_line(std::vector<std::string> &tokens);

private:
	std::ifstream scriptFile;
}; // ScriptReader

ScriptReader::ScriptReader(std::string filePath)
{
	open(filePath);
}

ScriptReader::~ScriptReader()
{
	scriptFile.close();
}

void ScriptReader::open(std::string filePath)
{
	// Close scriptFile if was open
	if (scriptFile)
		scriptFile.close();

	// Assumes file exists, later add check to see if it exists
	scriptFile.open(filePath.c_str());
	if (!scriptFile) {
		Assert(false, ExcMessage("File path passed to ScriptReader does not exist."))
	}
	std::cout << "Successfully open script: " << filePath << std::endl;
}

bool ScriptReader::get_next_line(std::vector<std::string> &tokens)
{
	if (!scriptFile) {
		Assert(false, ExcMessage("ScriptReader has tried to read the stream before a file has been opened."))
	}

	bool foundTokens = false;
	while (!foundTokens) {
		// Copy the line from the file to a string
		std::string line;
		getline(scriptFile, line);

		// Test if there was an error
		if (scriptFile.bad()) {
			Assert(false, ExcMessage("Error when reading stream."))
		}

		// Check end of file
		if (scriptFile.eof())
			return false;

		// Convert the line into tokens
		std::istringstream ss(line);
		std::istream_iterator<std::string> current(ss), end;
		tokens.clear();

		// Insert tokens but watch for comments (indicated by //)
		while (current != end) {
			// Test if the token contains the comment characters "//"
			if ((*current).find("//") != std::string::npos)
				break;

			// Push back the token
			tokens.push_back(*current);

			// Increment iterator
			current++;
		}

		// If there are no tokens then the line was just "/n" so get the next one
		if (tokens.size() > 0)
			foundTokens = true;
	}

	// get_next_line succeeded
	return true;
}

} // namespace

#endif /* SCRIPTREADER_H_ */
