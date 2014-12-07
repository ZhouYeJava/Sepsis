
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "DenseVector.h"
#include "SparseVector.h"
#include "Sparm.h"
#include "data_type_api.h"


//! convert string to variable of type T. Used to reading floats, int etc from files
template<typename T>
inline T Scan(const std::string &input)
{
	std::stringstream stream(input);
	T ret;
	stream >> ret;
	return ret;
}

//! Specialisation for performance
template<>
inline int Scan<int>(const std::string &input)
{
	return atoi(input.c_str());
}

//! Specialisation for performance
template<>
inline float Scan<float>(const std::string &input)
{
	return (float) atof(input.c_str());
}

//! just return input
template<>
inline std::string Scan<std::string>(const std::string &input)
{
	return input;
}

inline std::vector<std::string> Tokenize(const std::string& str, const std::string& delimiters = " \t")
{
	std::vector<std::string> tokens;
	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (std::string::npos != pos || std::string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;
}


void construct_intervals(string input_filename, Sparm &sparm); 
vector<Example> read_input_examples(string input_filename, Sparm &sparm);
DenseVector* LoadWeights(std::string filename, Sparm &sparm);






