/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/read_file_format.hh"
#include "clique.hh"
#include "configuration.hh"
#include "proof.hh"

#include <boost/program_options.hpp>

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>

#include <unistd.h>

namespace po = boost::program_options;

using std::boolalpha;
using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::function;
using std::localtime;
using std::make_optional;
using std::make_pair;
using std::make_shared;
using std::make_unique;
using std::put_time;
using std::string;
using std::string_view;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::operator""s;
using std::chrono::seconds;
using std::chrono::steady_clock;
using std::chrono::system_clock;

auto colour_class_order_from_string(string_view s) -> ColourClassOrder
{
    if (s == "colour")
        return ColourClassOrder::ColourOrder;
    else if (s == "singletons-first")
        return ColourClassOrder::SingletonsFirst;
    else if (s == "sorted")
        return ColourClassOrder::Sorted;
    else
        throw UnsupportedConfiguration{ "Unknown colour class order '" + string(s) + "'" };
}

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                         "Display help information")
            ("timeout",            po::value<int>(),         "Abort after this many seconds")
            ("format",             po::value<string>(),      "Specify input file format (auto, lad, labelledlad, dimacs)")
            ("decide",             po::value<int>(),         "Solve this decision problem");

        po::options_description configuration_options{ "Advanced configuration options" };
        configuration_options.add_options()
            ("colour-ordering",    po::value<string>(),      "Specify colour-ordering (colour / singletons-first / sorted)")
            ("input-order",                                  "Use the input order for colouring (usually a bad idea)")
            ("restarts-constant",  po::value<int>(),         "How often to perform restarts (disabled by default)")
            ("geometric-restarts", po::value<double>(),      "Use geometric restarts with the specified multiplier (default is Luby)");
        display_options.add(configuration_options);

        po::options_description proof_logging_options{ "Proof logging options" };
        proof_logging_options.add_options()
            ("prove",               po::value<string>(),       "Write unsat proofs to this filename (suffixed with .opb and .veripb)")
            ("proof-names",                                    "Use 'friendly' variable names in the proof, rather than x1, x2, ...")
            ("compress-proof",                                 "Compress the proof using bz2");
        display_options.add(proof_logging_options);

        po::options_description all_options{ "All options" };
        all_options.add_options()
            ("graph-file", "Specify the graph file")
            ;

        all_options.add(display_options);

        po::positional_options_description positional_options;
        positional_options
            .add("graph-file", 1)
            ;

        po::variables_map options_vars;
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(positional_options)
                .run(), options_vars);
        po::notify(options_vars);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << "Usage: " << argv[0] << " [options] graph-file" << endl;
            cout << endl;
            cout << display_options << endl;
            return EXIT_SUCCESS;
        }

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("graph-file")) {
            cout << "Usage: " << argv[0] << " [options] graph-file" << endl;
            return EXIT_FAILURE;
        }

        /* Figure out what our options should be. */
        CliqueParams params;

        if (options_vars.count("decide"))
            params.decide = make_optional(options_vars["decide"].as<int>());

        if (options_vars.count("restarts-constant")) {
            if (options_vars.count("geometric-restarts")) {
                double initial_value = GeometricRestartsSchedule::default_initial_value;
                double multiplier = options_vars["geometric-restarts"].as<double>();
                initial_value = options_vars["restarts-constant"].as<int>();
                params.restarts_schedule = make_unique<GeometricRestartsSchedule>(initial_value, multiplier);
            }
            else {
                long long multiplier = LubyRestartsSchedule::default_multiplier;
                multiplier = options_vars["restarts-constant"].as<int>();
                params.restarts_schedule = make_unique<LubyRestartsSchedule>(multiplier);
            }
        }
        else
            params.restarts_schedule = make_unique<NoRestartsSchedule>();

        if (options_vars.count("colour-ordering"))
            params.colour_class_order = colour_class_order_from_string(options_vars["colour-ordering"].as<string>());
        params.input_order = options_vars.count("input-order");

#if !defined(_WIN32)
        char hostname_buf[255];
        if (0 == gethostname(hostname_buf, 255))
            cout << "hostname = " << string(hostname_buf) << endl;
#endif
        cout << "commandline =";
        for (int i = 0 ; i < argc ; ++i)
            cout << " " << argv[i];
        cout << endl;

        auto started_at = system_clock::to_time_t(system_clock::now());
        cout << "started_at = " << put_time(localtime(&started_at), "%F %T") << endl;

        /* Read in the graphs */
        string pattern_format_name = options_vars.count("format") ? options_vars["format"].as<string>() : "auto";
        auto graph = read_file_format(pattern_format_name, options_vars["graph-file"].as<string>());

        cout << "file = " << options_vars["graph-file"].as<string>() << endl;

        if (options_vars.count("prove")) {
            bool friendly_names = options_vars.count("proof-names");
            bool compress_proof = options_vars.count("compress-proof");
            string fn = options_vars["prove"].as<string>();
            string suffix = compress_proof ? ".bz2" : "";
            params.proof = make_unique<Proof>(fn + ".opb", fn + ".veripb", friendly_names, compress_proof);
            cout << "proof_model = " << fn << ".opb" << suffix << endl;
            cout << "proof_log = " << fn << ".veripb" << suffix << endl;
        }

        /* Prepare and start timeout */
        params.timeout = make_shared<Timeout>(options_vars.count("timeout") ? seconds{ options_vars["timeout"].as<int>() } : 0s);

        /* Start the clock */
        params.start_time = steady_clock::now();

        auto result = solve_clique_problem(graph, params);

        /* Stop the clock. */
        auto overall_time = duration_cast<milliseconds>(steady_clock::now() - params.start_time);

        params.timeout->stop();

        cout << "status = ";
        if (params.timeout->aborted())
            cout << "aborted";
        else if (! result.clique.empty())
            cout << "true";
        else
            cout << "false";
        cout << endl;

        cout << "nodes = " << result.nodes << endl;

        if (! result.clique.empty()) {
            cout << "omega = " << result.clique.size() << endl;
            cout << "clique =";
            for (auto v : result.clique)
                cout << " " << graph.vertex_name(v);
            cout << endl;
        }

        cout << "runtime = " << overall_time.count() << endl;

        for (const auto & s : result.extra_stats)
            cout << s << endl;

        return EXIT_SUCCESS;
    }
    catch (const GraphFileError & e) {
        cerr << "Error: " << e.what() << endl;
        if (e.file_at_least_existed())
            cerr << "Maybe try specifying --format?" << endl;
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

