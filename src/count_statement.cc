#include <gss/counting_statement.hh>

#include <cstdio>
#include <fstream>
#include <iostream>

#include <cxxopts.hpp>

using namespace gss;

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::ofstream;
using std::string;

// Reads a counting statement (dev_docs/counting-statements.md), checks that it is
// well-formed, prints its count, and optionally writes its counting formula and its
// side-condition instance. The count is only right if the side-condition instance is
// unsatisfiable.
auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Count a counting statement", "Get started by using option --help");

        options.add_options("Program options") //
            ("help", "Display help information") //
            ("opb", "Write the counting formula to this file", cxxopts::value<string>()) //
            ("side-conditions", "Write the side-condition instance to this file (not written if there is nothing to prove)", cxxopts::value<string>());

        options.add_options() //
            ("statement-file", "Specify the statement file", cxxopts::value<string>());

        options.parse_positional({"statement-file"});

        auto options_vars = options.parse(argc, argv);

        if (options_vars.count("help")) {
            cout << options.help() << endl;
            return EXIT_SUCCESS;
        }

        if (! options_vars.count("statement-file")) {
            cout << "Usage: " << argv[0] << " [options] statement-file" << endl;
            return EXIT_FAILURE;
        }

        ifstream in{options_vars["statement-file"].as<string>()};
        if (! in)
            throw CountingStatementError{"cannot open " + options_vars["statement-file"].as<string>()};
        auto statement = read_counting_statement(in);
        auto count = count_counting_statement(statement);
        auto side_conditions = counting_side_conditions(statement);

        if (options_vars.count("opb")) {
            ofstream out{options_vars["opb"].as<string>()};
            write_counting_opb(out, statement);
            if (! out)
                throw CountingStatementError{"error writing " + options_vars["opb"].as<string>()};
        }

        if (options_vars.count("side-conditions")) {
            auto name = options_vars["side-conditions"].as<string>();
            bool wrote;
            {
                ofstream out{name};
                wrote = write_side_condition_opb(out, statement);
                if (! out)
                    throw CountingStatementError{"error writing " + name};
            }
            if (! wrote)
                std::remove(name.c_str());
        }

        cout << "variables = " << statement.number_of_variables << endl;
        cout << "values = " << statement.number_of_values << endl;
        cout << "nodes = " << statement.nodes.size() << endl;
        cout << "side_conditions = " << side_conditions.size() << endl;
        cout << "count = " << count << endl;

        return EXIT_SUCCESS;
    }
    catch (const cxxopts::exceptions::exception & e) {
        cerr << "Error: " << e.what() << endl;
        cerr << "Try " << argv[0] << " --help" << endl;
        return EXIT_FAILURE;
    }
    catch (const exception & e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}
