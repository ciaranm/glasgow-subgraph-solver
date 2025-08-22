#include <gss/formats/read_file_format.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/lackey.hh>
#include <gss/innards/symmetries.hh>
#include <gss/innards/verify.hh>
#include <gss/restarts.hh>
#include <gss/sip_decomposer.hh>

#include <cxxopts.hpp>

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include <unistd.h>

using namespace gss;

using std::boolalpha;
using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::function;
using std::list;
using std::localtime;
using std::make_pair;
using std::make_shared;
using std::make_unique;
using std::optional;
using std::pair;
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
        cxxopts::Options options("Glasgow Subgraph Solver", "Get started by using option --help");

        options.add_options("Program options")
            ("help", "Display help information")
            ("timeout", "Abort after this many seconds", cxxopts::value<int>())
            ("parallel", "Use auto-configured parallel search (highly nondeterministic runtimes)");


        options.add_options("Program options")
            ("noninjective", "Drop the injectivity requirement")
            ("locally-injective", "Require only local injectivity")
            ("induced", "Find an induced mapping")
            ("count-solutions", "Count the number of solutions")
            ("print-all-solutions", "Print out every solution, rather than one")
            ("solution-limit", "Stop after finding this many solutions (only when --print-all-solutions)", cxxopts::value<unsigned long long>());

        options.add_options("Input file options")
            ("format", "Specify input file format (auto, lad, vertexlabelledlad, labelledlad, dimacs)", cxxopts::value<string>())
            ("pattern-format", "Specify input file format just for the pattern graph", cxxopts::value<string>())
            ("target-format", "Specify input file format just for the target graph", cxxopts::value<string>());

        options.add_options("Advanced search configuration options")
            ("restarts", "Specify restart policy (luby / geometric / timed / none)", cxxopts::value<string>())
            ("geometric-multiplier", "Specify multiplier for geometric restarts", cxxopts::value<double>())
            ("geometric-constant", "Specify starting constant for geometric restarts", cxxopts::value<double>())
            ("restart-interval", "Specify the restart interval in milliseconds for timed restarts", cxxopts::value<int>())
            ("restart-minimum", "Specify a minimum number of backtracks before a timed restart can trigger", cxxopts::value<int>())
            ("luby-constant", "Specify the starting constant / multiplier for Luby restarts", cxxopts::value<int>())
            ("value-ordering", "Specify value-ordering heuristic (biased / degree / antidegree / random / none)", cxxopts::value<string>())
            ("pattern-symmetries", "Eliminate pattern symmetries (requires Gap)")
            ("target-symmetries", "Eliminate target symmetries (requires Gap)");

        options.add_options("Advanced input processing options")
            ("no-clique-detection", "Disable clique / independent set detection")
            ("no-supplementals", "Do not use supplemental graphs")
            ("no-nds", "Do not use neighbourhood degree sequences");

        options.add_options("Advanced parallelism options")
            ("threads", "Use threaded search, with this many threads (0 to auto-detect)", cxxopts::value<int>())
            ("triggered-restarts", "Have one thread trigger restarts (more nondeterminism, better performance)")
            ("delay-thread-creation", "Do not create threads until after the first restart");

        vector<string> pattern_less_thans, target_occur_less_thans;
        options.add_options("Manual symmetry options")
            ("pattern-less-than", "Specify a pattern less than constraint, in the form v<w",
                cxxopts::value<vector<string>>(pattern_less_thans))
            ("pattern-automorphism-group-size", "Specify the size of the pattern graph automorphism group", cxxopts::value<string>())
            ("target-occurs-less-than", "Specify a target occurs less than constraint, in the form v<w",
                cxxopts::value<vector<string>>(target_occur_less_thans))
            ("target-automorphism-group-size", "Specify the size of the target graph automorphism group", cxxopts::value<string>());

        options.add_options("External constraint solver options")
            ("send-to-lackey", "Send candidate solutions to an external solver over this named pipe", cxxopts::value<string>())
            ("receive-from-lackey", "Receive responses from external solver over this named pipe", cxxopts::value<string>())
            ("send-partials-to-lackey", "Send partial solutions to the lackey")
            ("propagate-using-lackey", "Propagate using lackey (never / root / root-and-backjump / always)", cxxopts::value<string>());

        options.add_options("Proof logging options")
            ("prove", "Write unsat proofs to this filename (suffixed with .opb and .pbp)", cxxopts::value<string>())
            ("verbose-proofs", "Write lots of comments to the proof, for tracing")
            ("recover-proof-encoding", "Recover the proof encoding, to work with verified encoders");

        vector<string> shapes;
        vector<int> shape_counts, shape_injectives;
        options.add_options("Hidden")
            ("enumerate", "Alias for --count-solutions (backwards compatibility)")
            ("distance3", "Use distance 3 filtering (experimental)")
            ("k4", "Use 4-clique filtering (experimental)")
            ("n-exact-path-graphs", "Specify number of exact path graphs", cxxopts::value<int>())
            ("decomposition", "Use decomposition")
            ("cliques", "Use clique size constraints")
            ("cliques-on-supplementals", "Use clique size constraints on supplemental graphs too")
            ("shape", "Specify an extra shape graph (slow, experimental)", cxxopts::value<std::vector<std::string>>())
            ("shape-count", "Specify how many times the shape must occur", cxxopts::value<std::vector<int>>())
            ("shape-injective", "Specify whether the shape must occur injectively", cxxopts::value<std::vector<int>>());

        options.add_options()
            ("pattern-file", "specify the pattern file", cxxopts::value<std::string>())
            ("target-file", "specify the target file", cxxopts::value<std::string>());

        options.parse_positional({"pattern-file", "target-file"});

        auto options_vars = options.parse(argc, argv);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << options.help() << endl;
            return EXIT_SUCCESS;
        }

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("pattern-file") || ! options_vars.count("target-file")) {
            cout << "Usage: " << argv[0] << " [options] pattern target" << endl;
            return EXIT_FAILURE;
        }

        /* Figure out what our options should be. */
        HomomorphismParams params;

        if (options_vars.count("noninjective") && options_vars.count("locally-injective")) {
            cerr << "Cannot specify both --noninjective and --locally-injective" << endl;
            return EXIT_FAILURE;
        }
        else if (options_vars.count("noninjective"))
            params.injectivity = Injectivity::NonInjective;
        else if (options_vars.count("locally-injective"))
            params.injectivity = Injectivity::LocallyInjective;
        else
            params.injectivity = Injectivity::Injective;

        params.induced = options_vars.count("induced");
        params.count_solutions = options_vars.count("count-solutions") || options_vars.count("enumerate") || options_vars.count("print-all-solutions");

        params.triggered_restarts = options_vars.count("triggered-restarts") || options_vars.count("parallel");

        if (options_vars.count("threads"))
            params.n_threads = options_vars["threads"].as<unsigned>();
        else if (options_vars.count("parallel"))
            params.n_threads = 0;

        if (options_vars.count("delay-thread-creation") || options_vars.count("parallel"))
            params.delay_thread_creation = true;

        if (options_vars.count("restarts")) {
            string restarts_policy = options_vars["restarts"].as<string>();
            if (restarts_policy == "luby") {
                unsigned long long multiplier = LubyRestartsSchedule::default_multiplier;
                if (options_vars.count("luby-constant"))
                    multiplier = options_vars["luby-constant"].as<int>();
                params.restarts_schedule = make_unique<LubyRestartsSchedule>(multiplier);
            }
            else if (restarts_policy == "geometric") {
                double geometric_constant = GeometricRestartsSchedule::default_initial_value;
                double geometric_multiplier = GeometricRestartsSchedule::default_multiplier;
                if (options_vars.count("geometric-constant"))
                    geometric_constant = options_vars["geometric-constant"].as<double>();
                if (options_vars.count("geometric-multiplier"))
                    geometric_multiplier = options_vars["geometric-multiplier"].as<double>();
                params.restarts_schedule = make_unique<GeometricRestartsSchedule>(geometric_constant, geometric_multiplier);
            }
            else if (restarts_policy == "timed") {
                milliseconds duration = TimedRestartsSchedule::default_duration;
                unsigned long long minimum_backtracks = TimedRestartsSchedule::default_minimum_backtracks;
                if (options_vars.count("restart-interval"))
                    duration = milliseconds{options_vars["restart-interval"].as<int>()};
                if (options_vars.count("restart-minimum"))
                    minimum_backtracks = options_vars["restart-minimum"].as<int>();
                params.restarts_schedule = make_unique<TimedRestartsSchedule>(duration, minimum_backtracks);
            }
            else if (restarts_policy == "none") {
                params.restarts_schedule = make_unique<NoRestartsSchedule>();
            }
            else {
                cerr << "Unknown restarts policy '" << restarts_policy << "'" << endl;
                return EXIT_FAILURE;
            }
        }
        else {
            if (params.count_solutions && ! options_vars.count("parallel"))
                params.restarts_schedule = make_unique<NoRestartsSchedule>();
            else if (options_vars.count("parallel"))
                params.restarts_schedule = make_unique<TimedRestartsSchedule>(TimedRestartsSchedule::default_duration, TimedRestartsSchedule::default_minimum_backtracks);
            else
                params.restarts_schedule = make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier);
        }

        if (options_vars.count("value-ordering")) {
            string value_ordering_heuristic = options_vars["value-ordering"].as<string>();
            if (value_ordering_heuristic == "none")
                params.value_ordering_heuristic = ValueOrdering::None;
            else if (value_ordering_heuristic == "biased")
                params.value_ordering_heuristic = ValueOrdering::Biased;
            else if (value_ordering_heuristic == "degree")
                params.value_ordering_heuristic = ValueOrdering::Degree;
            else if (value_ordering_heuristic == "antidegree")
                params.value_ordering_heuristic = ValueOrdering::AntiDegree;
            else if (value_ordering_heuristic == "random")
                params.value_ordering_heuristic = ValueOrdering::Random;
            else {
                cerr << "Unknown value-ordering heuristic '" << value_ordering_heuristic << "'" << endl;
                return EXIT_FAILURE;
            }
        }

        params.clique_detection = ! options_vars.count("no-clique-detection");
        params.distance3 = options_vars.count("distance3");
        params.k4 = options_vars.count("k4");
        if (options_vars.count("n-exact-path-graphs"))
            params.number_of_exact_path_graphs = options_vars["n-exact-path-graphs"].as<int>();
        params.no_supplementals = options_vars.count("no-supplementals");
        params.no_nds = options_vars.count("no-nds");
        params.clique_size_constraints = options_vars.count("cliques");
        params.clique_size_constraints_on_supplementals = options_vars.count("cliques-on-supplementals");

        if (options_vars.count("shape")) {
            for (decltype(shapes.size()) s = 0; s != shapes.size(); ++s) {
                auto graph = make_unique<InputGraph>(read_file_format("csv", shapes[s]));
                params.extra_shapes.emplace_back(
                    std::move(graph), // weird - issue when using std::move at headers
                    s >= shape_injectives.size() ? true : shape_injectives[s],
                    s >= shape_counts.size() ? 1 : shape_counts[s]);
            }
        }

        string pattern_automorphism_group_size = "1", target_automorphism_group_size = "1";
        bool was_given_pattern_automorphism_group = false, was_given_target_automorphism_group = false;
        if (options_vars.count("pattern-automorphism-group-size")) {
            pattern_automorphism_group_size = options_vars["pattern-automorphism-group-size"].as<string>();
            was_given_pattern_automorphism_group = true;
        }

        if (options_vars.count("target-automorphism-group-size")) {
            target_automorphism_group_size = options_vars["target-automorphism-group-size"].as<string>();
            was_given_target_automorphism_group = true;
        }

        for (auto & s : pattern_less_thans) {
            auto p = s.find('<');
            if (p == string::npos) {
                cerr << "Invalid pattern less-than constraint '" << s << "'" << endl;
                return EXIT_FAILURE;
            }
            auto a = s.substr(0, p), b = s.substr(p + 1);
            params.pattern_less_constraints.emplace_back(a, b);
        }

        for (auto & s : target_occur_less_thans) {
            auto p = s.find('<');
            if (p == string::npos) {
                cerr << "Invalid target occurs-less-than constraint '" << s << "'" << endl;
                return EXIT_FAILURE;
            }
            auto a = s.substr(0, p), b = s.substr(p + 1);
            params.target_occur_less_constraints.emplace_back(a, b);
        }

        if (options_vars.count("send-to-lackey") ^ options_vars.count("receive-from-lackey")) {
            cerr << "Must specify both of --send-to-lackey and --receive-from-lackey" << endl;
            return EXIT_FAILURE;
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
        string pattern_format_name = options_vars.count("pattern-format") ? options_vars["pattern-format"].as<string>() : default_format_name;
        string target_format_name = options_vars.count("target-format") ? options_vars["target-format"].as<string>() : default_format_name;
        auto pattern = read_file_format(pattern_format_name, options_vars["pattern-file"].as<string>());
        auto target = read_file_format(target_format_name, options_vars["target-file"].as<string>());

        cout << "pattern_file = " << options_vars["pattern-file"].as<string>() << endl;
        cout << "target_file = " << options_vars["target-file"].as<string>() << endl;

        if (options_vars.count("send-to-lackey") && options_vars.count("receive-from-lackey")) {
            auto lackey_started_at = steady_clock::now();
            params.lackey = make_unique<innards::Lackey>(
                options_vars["send-to-lackey"].as<string>(),
                options_vars["receive-from-lackey"].as<string>(),
                pattern, target);
            auto lackey_time = duration_cast<milliseconds>(steady_clock::now() - lackey_started_at);
            cout << "lackey_init_time = " << lackey_time.count() << endl;
        }
        params.send_partials_to_lackey = options_vars.count("send-partials-to-lackey");
        if (options_vars.count("propagate-using-lackey")) {
            string propagate_using_lackey = options_vars["propagate-using-lackey"].as<string>();
            if (propagate_using_lackey == "always")
                params.propagate_using_lackey = PropagateUsingLackey::Always;
            else if (propagate_using_lackey == "root")
                params.propagate_using_lackey = PropagateUsingLackey::Root;
            else if (propagate_using_lackey == "root-and-backjump")
                params.propagate_using_lackey = PropagateUsingLackey::RootAndBackjump;
            else if (propagate_using_lackey == "never")
                params.propagate_using_lackey = PropagateUsingLackey::Never;
            else {
                cerr << "Unknown propagate-using-lackey option '" << propagate_using_lackey << "'" << endl;
                return EXIT_FAILURE;
            }
        }
        else
            params.propagate_using_lackey = PropagateUsingLackey::Never;

        optional<unsigned long long> solutions_remaining;
        if (options_vars.contains("solution-limit"))
            solutions_remaining = options_vars["solution-limit"].as<unsigned long long>();

        if (options_vars.count("print-all-solutions")) {
            params.enumerate_callback = [&](const VertexToVertexMapping & mapping) -> bool {
                cout << "mapping = ";
                for (auto v : mapping)
                    cout << "(" << pattern.vertex_name(v.first) << " -> " << target.vertex_name(v.second) << ") ";
                cout << endl;

                return (! solutions_remaining) || (0 != --*solutions_remaining);
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

        cout << "pattern_properties =";
        describe(pattern);
        cout << "pattern_vertices = " << pattern.size() << endl;
        cout << "pattern_directed_edges = " << pattern.number_of_directed_edges() << endl;
        cout << "target_properties =";
        describe(target);
        cout << "target_vertices = " << target.size() << endl;
        cout << "target_directed_edges = " << target.number_of_directed_edges() << endl;

        /* Prepare and start timeout */
        params.timeout = make_shared<Timeout>(options_vars.count("timeout") ? seconds{options_vars["timeout"].as<int>()} : 0s); // is this a default val assignment? think yes.

        /* Start the clock */
        params.start_time = steady_clock::now();

        if (options_vars.count("pattern-symmetries")) {
            auto gap_start_time = steady_clock::now();
            innards::find_symmetries(argv[0], pattern, params.pattern_less_constraints, pattern_automorphism_group_size);
            was_given_pattern_automorphism_group = true;
            cout << "pattern_symmetry_time = " << duration_cast<milliseconds>(steady_clock::now() - gap_start_time).count() << endl;
            cout << "pattern_less_constraints =";
            for (auto & [a, b] : params.pattern_less_constraints)
                cout << " " << a << "<" << b;
            cout << endl;
        }

        if (was_given_pattern_automorphism_group)
            cout << "pattern_automorphism_group_size = " << pattern_automorphism_group_size << endl;

        if (options_vars.count("target-symmetries")) {
            auto gap_start_time = steady_clock::now();
            innards::find_symmetries(argv[0], target, params.target_occur_less_constraints, target_automorphism_group_size);
            was_given_target_automorphism_group = true;
            cout << "target_symmetry_time = " << duration_cast<milliseconds>(steady_clock::now() - gap_start_time).count() << endl;
            cout << "target_occur_less_constraints =";
            for (auto & [a, b] : params.target_occur_less_constraints)
                cout << " " << a << "<" << b;
            cout << endl;
        }

        if (was_given_target_automorphism_group)
            cout << "target_automorphism_group_size = " << target_automorphism_group_size << endl;

        auto result = options_vars.count("decomposition") ? solve_sip_by_decomposition(pattern, target, params) : solve_homomorphism_problem(pattern, target, params);

        /* Stop the clock. */
        auto overall_time = duration_cast<milliseconds>(steady_clock::now() - params.start_time);

        cout << "status = ";
        if (params.timeout->aborted() || (solutions_remaining && 0 == *solutions_remaining))
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
                cout << "(" << pattern.vertex_name(v.first) << " -> " << target.vertex_name(v.second) << ") ";
            cout << endl;
        }

        cout << "runtime = " << overall_time.count() << endl;

        for (const auto & s : result.extra_stats)
            cout << s << endl;

        if (params.lackey) {
            cout << "lackey_calls = " << params.lackey->number_of_calls() << endl;
            cout << "lackey_checks = " << params.lackey->number_of_checks() << endl;
            cout << "lackey_deletions = " << params.lackey->number_of_deletions() << endl;
            cout << "lackey_propagations = " << params.lackey->number_of_propagations() << endl;
        }

        innards::verify_homomorphism(pattern, target, params.injectivity == Injectivity::Injective,
            params.injectivity == Injectivity::LocallyInjective, params.induced, result.mapping);

        return EXIT_SUCCESS;
    }
    catch (const GraphFileError & e) {
        cerr << "Error: " << e.what() << endl;
        if (e.file_at_least_existed())
            cerr << "Maybe try specifying one of --format, --pattern-format, or --target-format?" << endl;
        return EXIT_FAILURE;
    }
    catch (const cxxopts::exceptions::exception & e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Try " << argv[0] << " --help" << std::endl;
        return EXIT_FAILURE;
    }
    catch (const exception & e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}
