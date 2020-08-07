/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/bigraph.hh"
#include "homomorphism.hh"
#include "lackey.hh"
#include "symmetries.hh"
#include "proof.hh"
#include "restarts.hh"

#include <boost/program_options.hpp>

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include <unistd.h>

namespace po = boost::program_options;

using std::boolalpha;
using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::function;
using std::ifstream;
using std::localtime;
using std::make_pair;
using std::make_shared;
using std::make_unique;
using std::put_time;
using std::string;
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::operator""s;
using std::chrono::seconds;
using std::chrono::steady_clock;
using std::chrono::system_clock;

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                         "Display help information")
            ("timeout",            po::value<int>(),         "Abort after this many seconds");

        po::options_description problem_options{ "Problem options" };
        problem_options.add_options()
            ("count-solutions",                              "Count the number of solutions")
            ("print-all-solutions",                          "Print out every solution, rather than one");
        display_options.add(problem_options);

        po::options_description mangling_options{ "Advanced input processing options" };
        mangling_options.add_options()
            ("no-supplementals",                               "Do not use supplemental graphs")
            ("no-nds",                                         "Do not use neighbourhood degree sequences");
        display_options.add(mangling_options);

        po::options_description all_options{ "All options" };
        all_options.add_options()
            ("pattern-file", "specify the pattern file")
            ("target-file",  "specify the target file")
            ;

        all_options.add(display_options);

        po::positional_options_description positional_options;
        positional_options
            .add("pattern-file", 1)
            .add("target-file", 1)
            ;

        po::variables_map options_vars;
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(positional_options)
                .run(), options_vars);
        po::notify(options_vars);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << "Usage: " << argv[0] << " [options] pattern target" << endl;
            cout << endl;
            cout << display_options << endl;
            return EXIT_SUCCESS;
        }

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("pattern-file") || ! options_vars.count("target-file")) {
            cout << "Usage: " << argv[0] << " [options] pattern target" << endl;
            return EXIT_FAILURE;
        }

        /* Figure out what our options should be. */
        HomomorphismParams params;

        params.injectivity = Injectivity::Injective;
        params.induced = true;
        params.bigraph = true;
        params.count_solutions = options_vars.count("count-solutions") || options_vars.count("print-all-solutions");

        if (params.count_solutions)
            params.restarts_schedule = make_unique<NoRestartsSchedule>();
        else
            params.restarts_schedule = make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier);

        char hostname_buf[255];
        if (0 == gethostname(hostname_buf, 255))
            cout << "hostname = " << string(hostname_buf) << endl;
        cout << "commandline =";
        for (int i = 0 ; i < argc ; ++i)
            cout << " " << argv[i];
        cout << endl;

        auto started_at = system_clock::to_time_t(system_clock::now());
        cout << "started_at = " << put_time(localtime(&started_at), "%F %T") << endl;

        /* Read in the graphs */

        string pattern_filename = options_vars["pattern-file"].as<string>();
        string target_filename = options_vars["target-file"].as<string>();

        ifstream pattern_infile{ pattern_filename };
        if (! pattern_infile)
            throw GraphFileError{ pattern_filename, "unable to open pattern file", false };

        ifstream target_infile{ target_filename };
        if (! target_infile)
            throw GraphFileError{ target_filename, "unable to open target file", false };

        auto graphs = make_pair(
            read_pattern_bigraph(move(pattern_infile), pattern_filename),
            read_target_bigraph(move(target_infile), target_filename));

        cout << "pattern_file = " << pattern_filename << endl;
        cout << "target_file = " << target_filename << endl;

        if (options_vars.count("print-all-solutions")) {
            params.enumerate_callback = [&] (const VertexToVertexMapping & mapping) {
                cout << "mapping = ";
                for (auto v : mapping)
                    cout << "(" << graphs.first.vertex_name(v.first) << " -> " << graphs.second.vertex_name(v.second) << ") ";
                cout << endl;
            };
        }

        /* Prepare and start timeout */
        params.timeout = make_shared<Timeout>(options_vars.count("timeout") ? seconds{ options_vars["timeout"].as<int>() } : 0s);

        /* Start the clock */
        params.start_time = steady_clock::now();

        auto result = solve_homomorphism_problem(graphs.first, graphs.second, params);

        /* Stop the clock. */
        auto overall_time = duration_cast<milliseconds>(steady_clock::now() - params.start_time);

        cout << "status = ";
        if (params.timeout->aborted())
            cout << "aborted";
        else if ((! result.mapping.empty()) || (params.count_solutions && result.solution_count > 0))
            cout << "true";
        else
            cout << "false";
        cout << endl;

        if (params.count_solutions)
            cout << "solution_count = " << result.solution_count << endl;

        cout << "nodes = " << result.nodes << endl;
        cout << "propagations = " << result.propagations << endl;

        if (! result.mapping.empty() && ! options_vars.count("print-all-solutions")) {
            cout << "mapping = ";
            for (auto v : result.mapping)
                cout << "(" << graphs.first.vertex_name(v.first) << " -> " << graphs.second.vertex_name(v.second) << ") ";
            cout << endl;
        }

        cout << "runtime = " << overall_time.count() << endl;

        for (const auto & s : result.extra_stats)
            cout << s << endl;

        return EXIT_SUCCESS;
    }
    catch (const GraphFileError & e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
    catch (const po::error & e) {
        cerr << "Error: " << e.what() << endl;
        cerr << "Try " << argv[0] << " --help" << endl;
        return EXIT_FAILURE;
    }
    catch (const exception & e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}

