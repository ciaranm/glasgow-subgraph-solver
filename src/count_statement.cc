#include <gss/counting_statement.hh>

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
// well-formed, prints its count, and optionally writes its tagged formula, which has
// exactly that many solutions.
auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Count a counting statement", "Get started by using option --help");

        options.add_options("Program options") //
            ("help", "Display help information") //
            ("opb", "Write the tagged formula to this file", cxxopts::value<string>());

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

        if (options_vars.count("opb")) {
            ofstream out{options_vars["opb"].as<string>()};
            write_tagged_opb(out, statement);
            if (! out)
                throw CountingStatementError{"error writing " + options_vars["opb"].as<string>()};
        }

        cout << "variables = " << statement.domain_sizes.size() << endl;
        cout << "nodes = " << statement.nodes.size() << endl;
        cout << "count = " << statement.nodes.at(statement.root).count << endl;

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
