// gemc headers
#include "string_utilities.h"

// mlibrary
#include "gstring.h"
using namespace gstring;

// C++ headers
#include <algorithm>

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

vector<string> get_info(string input, string chars)
{
	string stripped = replaceCharInStringWithChars(input, chars, " ");
	
	return getStringVectorFromString(stripped);
}

vector<string> get_info(string input)
{
	
	string chars("(),\"");
	
	string stripped = replaceCharInStringWithChars(input, chars, " ");
	
	return getStringVectorFromString(stripped);
}




// returns a vector of strings from a stringstream, space is delimeter
// ignore instances of second string
vector<string> get_strings_except(string input, string ignore)
{
	vector<string> pvalues;
	stringstream plist(input);
	while(!plist.eof())
	{
		string tmp;
		plist >> tmp;
		
		if(tmp.find(ignore) == string::npos)
			pvalues.push_back(trimSpacesFromString(tmp));
	}
	
	return pvalues;
}


double scan_number(const char *str)
{
	// Scan the c_string str for numbers only, then return the value as a float.
	// The str is not allowed to have spaces or alphanumerics, only 0-9 and .
	int i=0;
	while(char c=str[i++]) if(isalpha(c) && !(c=='-' || c=='+' || c=='e' || c=='E') )
	{
		cout << "WARNING: Unexpected Alphanumberic character found in number string:" << str << endl;
	}
	
	return( stringToDouble(str));
}






void print_vstring(vector<string> s)
{
	for(unsigned int i=0; i<s.size(); i++)
		cout << "string element: " << i << "  content: " << s[i] << endl;
}



string get_variation(string var)
{
	string variation = var;
	vector<string> vars = getStringVectorFromString(var);
	if(vars[0] == "main" && vars.size() > 1)
	{
		variation = vars[1];
	}
	return variation;
}

bool is_main_variation(string var)
{
	if(var.find("main:") == 0)
		return 1;
	
	return 0;
}



// overloading << for map<string, string>
ostream &operator<<(ostream &stream, map<string, string> smap)
{
	cout << endl;
	for(map<string, string>::iterator it = smap.begin(); it != smap.end(); it++)
		cout << it->first << " " << it->second << endl;
	
	cout << endl;
	
	return stream;
}





