#include <gss/counting_statement.hh>

#include <cxxopts.hpp>

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace gss;

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::ofstream;
using std::string;

// The standalone half of the verification model in dev_docs/decomposition-counting.md:
// reads a counting-language statement, checks what it can, prints the count, and
// optionally writes the PB translation that a proof is checked against. It knows
// nothing about graphs or inclusion-exclusion.
auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Count a counting-language statement", "Get started by using option --help");

        options.add_options("Program options") //
            ("help", "Display help information") //
            ("opb", "Write the PB translation to this file", cxxopts::value<string>());

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

        auto file_name = options_vars["statement-file"].as<string>();
        ifstream in{file_name};
        if (! in)
            throw std::runtime_error{"cannot open " + file_name};
        auto statement = read_counting_statement(in);
        auto result = count_counting_statement(statement);

        if (options_vars.count("opb")) {
            ofstream opb{options_vars["opb"].as<string>()};
            write_counting_statement_opb(opb, statement);
            if (! opb)
                throw std::runtime_error{"error writing " + options_vars["opb"].as<string>()};
        }

        long disjoint = 0, implies = 0;
        for (auto & o : result.obligations)
            ++(o.kind == CountingObligation::Kind::Disjoint ? disjoint : implies);

        cout << "count = " << result.count << endl;
        cout << "nodes = " << statement.nodes.size() << endl;
        // the count is only right if these hold, and nothing here has checked them
        cout << "unchecked_disjoint_sums = " << disjoint << endl;
        cout << "unchecked_differences = " << implies << endl;

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
