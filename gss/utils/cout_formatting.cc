#include <gss/homomorphism.hh>

#include <iostream>

using std::cout;
using std::endl;
using std::string;

void format_cout_with_string_value(string key, string value, bool json_output)
{
    if (json_output)
        cout << "\"" + key + "\"" << ": \"" << value << "\"," << endl;
    else
        cout << key << " = " << value << endl;
};

void format_cout_with_int_value(string key, int value, bool json_output)
{
    if (json_output)
        cout << "\"" + key + "\"" << ": " << value << "," << endl;
    else
        cout << key << " = " << value << endl;
};

string format_extra_results(string key, string value, bool json_output)
{
    if (json_output)
        return "\"" + key + "\"" + ": " + value + ",";
    return key + " = " + value;
}