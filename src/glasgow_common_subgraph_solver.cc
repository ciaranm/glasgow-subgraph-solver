#include <gss/common_subgraph.hh>
#include <gss/formats/read_file_format.hh>

#include <cxxopts.hpp>

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <fstream>

#include <unistd.h>
#include <nlohmann/json.hpp>

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
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::operator""s;
using std::chrono::seconds;
using std::chrono::steady_clock;
using std::chrono::system_clock;

using json = nlohmann::json;

auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Glasgow Subgraph Solver - but different???", "Get started by using option --help");

        options.add_options("Program options")
            ("help", "Display help information")
            ("timeout", "Abort after this many seconds", cxxopts::value<int>())
            ("decide", "Solve this decision problem", cxxopts::value<int>())
            ("count-solutions", "Count the number of solutions (--decide only)")
            ("print-all-solutions", "Print out every solution, rather than one (--decide only)")
            ("connected", "Only find connected graphs")
            ("clique", "Use the clique solver")
            ("json-output", "Dumps results to json file - takes filename as arg", cxxopts::value<string>());

        options.add_options("Input file options")
            ("format", "Specify input file format (auto, lad, vertexlabelledlad, labelledlad, dimacs)",  cxxopts::value<string>())
            ("first-format", "Specify input file format just for the first graph", cxxopts::value<string>())
            ("second-format", "Specify input file format just for the second graph", cxxopts::value<string>());

        options.add_options("Proof Logging Options")
            ("prove", "Write unsat proofs to this filename (suffixed with .opb and .pbp)", cxxopts::value<string>())
            ("verbose-proofs", "Write lots of comments to the proof, for tracing")
            ("recover-proof-encoding", "Recover the proof encoding, to work with verified encoders");

        options.add_options()
            ("first-file", "Specify the first graph file", cxxopts::value<string>())
            ("second-file", "Specify the second graph file", cxxopts::value<string>());

        options.parse_positional({"first-file", "second-file"});

        auto options_vars = options.parse(argc, argv);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << options.help() << endl;
            return EXIT_SUCCESS;
        }

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("first-file") || ! options_vars.count("second-file")) {
            cout << "Usage: " << argv[0] << " [options] graph-file graph-file" << endl;
            return EXIT_FAILURE;
        }

        /* Figure out what our options should be. */
        CommonSubgraphParams params;

        if (options_vars.count("decide"))
            params.decide = make_optional(options_vars["decide"].as<int>());

        params.connected = options_vars.count("connected");
        params.count_solutions = options_vars.count("count-solutions") || options_vars.count("print-all-solutions");
        params.clique = options_vars.count("clique");
        if (options_vars.count("json-output")) {
            params.json_output = options_vars["json-output"].as<std::string>();
        }

#if ! defined(__WIN32)
        char hostname_buf[255];
        if (0 == gethostname(hostname_buf, 255))
            cout << "hostname = " << string(hostname_buf) << endl;
#endif
        cout << "commandline =";
        for (int i = 0; i < argc; ++i)
            cout << " " << argv[i];
        cout << endl;

        auto started_at = system_clock::to_time_t(system_clock::now());
        cout << "started_at = " << put_time(localtime(&started_at), "%F %T") << endl;

        /* Read in the graphs */
        string default_format_name = options_vars.count("format") ? options_vars["format"].as<string>() : "auto";
        string first_format_name = options_vars.count("first-format") ? options_vars["first-format"].as<string>() : default_format_name;
        string second_format_name = options_vars.count("second-format") ? options_vars["second-format"].as<string>() : default_format_name;
        auto first = read_file_format(first_format_name, options_vars["first-file"].as<string>());
        auto second = read_file_format(second_format_name, options_vars["second-file"].as<string>());

        cout << "first_file = " << options_vars["first-file"].as<string>() << endl;
        cout << "second_file = " << options_vars["second-file"].as<string>() << endl;

        auto describe = [&] (const InputGraph & g) {
            if (g.directed())
                cout << " directed";
            if (g.loopy())
                cout << " loopy";
            if (g.has_vertex_labels())
                cout << " vertex_labels";
            if (g.has_edge_labels())
                cout << " edge_labels";
            cout << endl;
        };

        cout << "first_properties =";
        describe(first);
        cout << "first_vertices = " << first.size() << endl;
        cout << "first_directed_edges = " << first.number_of_directed_edges() << endl;
        cout << "second_properties =";
        describe(second);
        cout << "second_vertices = " << second.size() << endl;
        cout << "second_directed_edges = " << second.number_of_directed_edges() << endl;

        vector<vector<vector<string>>> all_mappings;
        if (options_vars.count("print-all-solutions")) {

            params.enumerate_callback = [&](const VertexToVertexMapping & mapping) {
                if (options_vars.count("json-output")) {
                    vector<vector<string>> mapping_vec;
                    mapping_vec.reserve(mapping.size());
                    for (auto & v : mapping) {
                        mapping_vec.push_back({
                            first.vertex_name(v.first),
                            second.vertex_name(v.second)
                        });
                    }
                    all_mappings.push_back(std::move(mapping_vec));
                } else {
                    cout << "mapping = ";
                    for (auto & v : mapping)
                        cout << "(" << first.vertex_name(v.first) << " -> " << second.vertex_name(v.second) << ") ";
                    cout << endl;
                }
            };
        }

        if (options_vars.count("prove")) {
            string fn = options_vars["prove"].as<string>();
            ProofOptions proof_options{
                .opb_file = fn + ".opb",
                .log_file = fn + ".pbp",
                .recover_encoding = options_vars.contains("recover-proof-encoding"),
                .super_extra_verbose = options_vars.contains("verbose-proofs")};
            params.proof_options = proof_options;
            cout << "proof_model = " << fn << ".opb" << endl;
            cout << "proof_log = " << fn << ".pbp" << endl;
        }

        /* Prepare and start timeout */
        params.timeout = make_shared<Timeout>(options_vars.count("timeout") ? seconds{options_vars["timeout"].as<int>()} : 0s);

        /* Start the clock */
        params.start_time = steady_clock::now();

        auto result = solve_common_subgraph_problem(first, second, params);

        /* Stop the clock. */
        auto overall_time = duration_cast<milliseconds>(steady_clock::now() - params.start_time);

        params.timeout->stop();


        const std::string status =
            params.timeout->aborted() ? "aborted" :
            ((! result.mapping.empty()) || (params.count_solutions && result.solution_count > 0)) ? "true" :
            "false";
        cout << "status = " <<status << endl;

        if (params.count_solutions)
            cout << "solution_count = " << result.solution_count << endl;

        cout << "nodes = " << result.nodes << endl;

        if (! result.mapping.empty() && ! options_vars.count("print-all-solutions")) {
            cout << "size = " << result.mapping.size() << endl;
            cout << "mapping = ";
            for (auto v : result.mapping)
                cout << "(" << first.vertex_name(v.first) << " -> " << second.vertex_name(v.second) << ") ";
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
