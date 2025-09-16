#include <gss/clique.hh>
#include <gss/configuration.hh>
#include <gss/formats/read_file_format.hh>
#include <gss/utils/cout_formatting.hh>

#include <cxxopts.hpp>

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>

#include <unistd.h>

using namespace gss;

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
        throw UnsupportedConfiguration{"Unknown colour class order '" + string(s) + "'"};
}

auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Glasgow Clique Solver", "Get started by using option --help");

        options.add_options("Program options")
            ("help", "Display help information")
            ("timeout", "Abort after this many seconds", cxxopts::value<int>())
            ("format", "Specify input file format (auto, lad, labelledlad, dimacs)", cxxopts::value<string>())
            ("decide", "Solve this decision problem", cxxopts::value<int>())
            ("json", "Changes the stdout format to json");

        options.add_options("Advanced configuration options")
            ("colour-ordering", "Specify colour-ordering (colour / singletons-first / sorted)", cxxopts::value<string>())
            ("input-order", "Use the input order for colouring (usually a bad idea)")
            ("restarts-constant", "How often to perform restarts (disabled by default)", cxxopts::value<int>())
            ("geometric-restarts", "Use geometric restarts with the specified multiplier (default is Luby)", cxxopts::value<double>());

        options.add_options("Proof logging options")
            ("prove", "Write unsat proofs to this filename (suffixed with .opb and .pbp)", cxxopts::value<string>())
            ("verbose-proofs", "Write lots of comments to the proof, for tracing")
            ("recover-proof-encoding", "Recover the proof encoding, to work with verified encoders");

        options.add_options()
            ("graph-file", "Specify the graph file", cxxopts::value<string>());

        options.parse_positional({"graph-file"});

        auto options_vars = options.parse(argc, argv);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << options.help() << endl;
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
        if (options_vars.count("json"))
            params.json_output = true;

        if (params.json_output)
            cout << "{" << endl;
#if ! defined(_WIN32)
        char hostname_buf[255];
        if (0 == gethostname(hostname_buf, 255))
            format_cout_with_string_value("hostname", string(hostname_buf), params.json_output);
#endif
        string command;
        for (int i = 0; i < argc; ++i)
            command += argv[i];
        format_cout_with_string_value("command", command, params.json_output);

        auto started_at = system_clock::to_time_t(system_clock::now());
        std::ostringstream oss;
        oss << put_time(localtime(&started_at), "%F %T");
        format_cout_with_string_value("started_at", oss.str(), params.json_output);

        /* Read in the graphs */
        string pattern_format_name = options_vars.count("format") ? options_vars["format"].as<string>() : "auto";
        auto graph = read_file_format(pattern_format_name, options_vars["graph-file"].as<string>());

        format_cout_with_string_value("file", options_vars["graph-file"].as<string>(), params.json_output);

        if (options_vars.count("prove")) {
            string fn = options_vars["prove"].as<string>();
            ProofOptions proof_options{
                .opb_file = fn + ".opb",
                .log_file = fn + ".pbp",
                .recover_encoding = options_vars.contains("recover-proof-encoding"),
                .super_extra_verbose = options_vars.contains("verbose-proofs")};
            params.proof_options = proof_options;
            format_cout_with_string_value("proof_model", fn + ".opb", params.json_output);
            format_cout_with_string_value("proof_log", fn + ".pbp", params.json_output);
        }

        /* Prepare and start timeout */
        params.timeout = make_shared<Timeout>(options_vars.count("timeout") ? seconds{options_vars["timeout"].as<int>()} : 0s);

        /* Start the clock */
        params.start_time = steady_clock::now();

        auto result = solve_clique_problem(graph, params);

        /* Stop the clock. */
        auto overall_time = duration_cast<milliseconds>(steady_clock::now() - params.start_time);

        params.timeout->stop();

        const std::string status =
           params.timeout->aborted() ? "aborted" :
           ! result.clique.empty() ? "true" :
           "false";
        format_cout_with_string_value("status", status, params.json_output);

        format_cout_with_int_value("nodes", result.nodes, params.json_output);

        if (! result.clique.empty()) {
            format_cout_with_int_value("omega", result.clique.size(), params.json_output);
            string value;
            bool first = true;
            for (auto v : result.clique)
                if (first) {
                    value = graph.vertex_name(v);
                    first = false;
                }
                else
                    value += " " + graph.vertex_name(v);
        }

        format_cout_with_int_value("runtime", overall_time.count(), params.json_output);
        for (const auto & s : result.extra_stats)
            cout << s << endl;

        if (params.json_output)
            cout << "}" << endl;

        return EXIT_SUCCESS;
    }
    catch (const GraphFileError & e) {
        cerr << "Error: " << e.what() << endl;
        if (e.file_at_least_existed())
            cerr << "Maybe try specifying --format?" << endl;
        return EXIT_FAILURE;
    }
    catch (const cxxopts::exceptions::exception& e) {
        cerr << "Error: " << e.what() << endl;
        cerr << "Try " << argv[0] << " --help" << endl;
        return EXIT_FAILURE;
    }
    catch (const exception & e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}
