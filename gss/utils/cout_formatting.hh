#ifndef GLASGOW_SUBGRAPH_SOLVER_COUT_FORMATTING_HH
#define GLASGOW_SUBGRAPH_SOLVER_COUT_FORMATTING_HH
#include <__fwd/string.h>

using std::string;

void format_cout_with_string_value(string key, string value, bool json_output);

void format_cout_with_int_value(string key, int value, bool json_output);

string format_extra_results(string key, string value, bool json_output);


#endif // GLASGOW_SUBGRAPH_SOLVER_COUT_FORMATTING_HH