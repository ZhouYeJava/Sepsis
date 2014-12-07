
#include <string>
#include <sstream>

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

//! Specialisation to understand yes/no y/n true/false 0/1
template<>
bool Scan<bool>(const std::string &input);

//! convert vectors of string to vectors of type T variables
template<typename T>
inline std::vector<T> Scan(const std::vector< std::string > &input)
{
	std::vector<T> output(input.size());
	for (size_t i = 0 ; i < input.size() ; i++)
	{
		output[i] = Scan<T>( input[i] );
	}
	return output;
}


/** tokenise input string to vector of string. each element has been separated by a character in the delimiters argument. 
		The separator can only be 1 character long. The default delimiters are space or tab
*/
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

//! tokenise input string to vector of type T
template<typename T>
inline std::vector<T> Tokenize( const std::string &input, const std::string& delimiters = " \t") 
{
	std::vector<std::string> stringVector = Tokenize(input, delimiters);
	return Scan<T>( stringVector );
}
