#include <gss/counting_statement.hh>
#include <gss/counting_statement_views.hh>
#include <gss/formats/read_file_format.hh>

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
using std::vector;

// Reads a counting statement (dev_docs/counting-statements.md), checks that it is
// well-formed, prints its count, and optionally writes its tagged formula, which has
// exactly that many solutions. It can also read the statement as a compressed form of the
// solutions it counts: render it, list them, or sample them uniformly
// (gss/counting_statement_views.hh).
auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Count a counting statement", "Get started by using option --help");

        options.add_options("Program options") //
            ("help", "Display help information") //
            ("opb", "Write the tagged formula to this file", cxxopts::value<string>()) //
            ("render", "Print the statement as a compressed form of its solutions") //
            ("list", "Print every solution") //
            ("solution", "Print the solution with this index, from 0", cxxopts::value<string>()) //
            ("sample", "Print this many solutions, each uniformly at random", cxxopts::value<unsigned long>()) //
            ("seed", "Random seed for --sample", cxxopts::value<unsigned long>()->default_value("0")) //
            ("format", "With --pattern and --target, the format of the graph files, to print vertex names", cxxopts::value<string>()->default_value("auto")) //
            ("pattern", "The pattern graph, to print its vertex names", cxxopts::value<string>()) //
            ("target", "The target graph, to print its vertex names", cxxopts::value<string>());

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

        CountingStatementNames names;
        if (options_vars.count("pattern")) {
            auto pattern = read_file_format(options_vars["format"].as<string>(), options_vars["pattern"].as<string>());
            if (unsigned(pattern.size()) != statement.domain_sizes.size())
                throw CountingStatementError{"the pattern does not have one vertex per variable"};
            vector<string> pattern_names;
            for (int u = 0; u < pattern.size(); ++u)
                pattern_names.push_back(pattern.vertex_name(u));
            names.variable = [pattern_names](int u) { return pattern_names.at(u); };
        }
        if (options_vars.count("target")) {
            auto target = read_file_format(options_vars["format"].as<string>(), options_vars["target"].as<string>());
            vector<string> target_names;
            for (int t = 0; t < target.size(); ++t)
                target_names.push_back(target.vertex_name(t));
            names.value = [target_names](int t) { return target_names.at(t); };
        }

        if (options_vars.count("render"))
            render_counting_statement(cout, statement, names);

        if (options_vars.count("solution")) {
            loooong index;
            if (0 != mpz_set_str(index.raw(), options_vars["solution"].as<string>().c_str(), 10))
                throw CountingStatementError{"--solution needs a number"};
            write_counting_object(cout, unrank_counting_statement(statement, index), names);
        }

        if (options_vars.count("list"))
            for (loooong i{0}; i < statement.nodes.at(statement.root).count; ++i)
                write_counting_object(cout, unrank_counting_statement(statement, i), names);

        if (options_vars.count("sample")) {
            CountingStatementSampler sample{statement, options_vars["seed"].as<unsigned long>()};
            for (unsigned long i = 0; i < options_vars["sample"].as<unsigned long>(); ++i)
                write_counting_object(cout, sample(), names);
        }

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
