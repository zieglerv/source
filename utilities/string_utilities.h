/// \file string_utilities.h
/// Set of string manipulation functions: \n
/// \author \n Maurizio Ungaro
/// \author mail: ungaro@jlab.org\n\n\n

#ifndef string_utilities_H
#define string_utilities_H

// Qt4 headers
#include <QString>
#include <QVariant>

// C++ headers
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

// mlibrary
#include "gstring.h"
using namespace gstring;

inline string stringify(double x)
{
	ostringstream o;
	o << x;
	return o.str();
}

inline string stringify(int x)
{
	ostringstream o;
	o << x;
	return o.str();
}





vector< vector<string> > dimensionstype(string);    ///< Returns dimensions nomenclature for different solid type
vector<string> get_strings_except(string, string);  ///< returns a vector of strings from a stringstream, space is delimiter, ignore string with second argument
void print_vstring(vector<string>);                 ///< prints each element of a string vector
vector<string> get_info(string);                    ///< get information from strings such as "5*GeV, 2*deg, 10*deg", parses out parenthesis, commas, quotes
vector<string> get_info(string, string);            ///< get information from strings such as "5*GeV, 2*deg, 10*deg", parses out strings in second argument
string get_variation(string);                       ///< parse variation name from string
bool is_main_variation(string);                     ///< returns 1 if the string "main:" is found on the input

ostream &operator<<(ostream &stream, map<string, string>);  ///< overload << for map<string, string>


// returns a double from a QVariant
inline double get_number(QVariant input)
{
	return get_number(trimSpacesFromString(qv_tostring(input)));
}








