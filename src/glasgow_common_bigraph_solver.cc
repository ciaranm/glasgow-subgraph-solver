#include "formats/common_bigraph.hh"
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
            ("print-all-solutions",                          "Print out every solution, rather than one");
        display_options.add(problem_options);

        po::options_description all_options{ "All options" };
        all_options.add_options()
            ("bigraph-1", "specify the first bigraph file")
            ("bigraph-2",  "specify the second bigraph file")
            ;
        all_options.add(display_options);

        po::positional_options_description positional_options;
        positional_options
            .add("bigraph-1", 1)
            .add("bigraph-2", 1)
            ;

        po::variables_map options_vars;
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(positional_options)
                .run(), options_vars);
        po::notify(options_vars);

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("bigraph-1") || ! options_vars.count("bigraph-2")) {
            cout << "Usage: " << argv[0] << " [options] bigraph-1 bigraph-2" << endl;
            return EXIT_FAILURE;
        }

        /* Figure out what our options should be. */
        HomomorphismParams params;

        params.injectivity = Injectivity::Injective;
        params.induced = false;
        params.bigraph = true;
        params.count_solutions = false;
        params.use_bigraph_projection_nogoods = false;
        params.no_nds = false;
        params.no_supplementals = true;
        params.equality_check = false;
        params.directed = false;

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

        string bigraph_1_filename = options_vars["bigraph-1"].as<string>();
        string bigraph_2_filename = options_vars["bigraph-2"].as<string>();

        ifstream bigraph_1_infile{ bigraph_1_filename };
        if (! bigraph_1_infile)
            throw GraphFileError{ bigraph_1_filename, "unable to open bigraph-1 file", false };

        ifstream bigraph_2_infile{ bigraph_2_filename };
        if (! bigraph_2_infile)
            throw GraphFileError{ bigraph_2_filename, "unable to open target file", false };

        Bigraph big1 = free_all_entities(read_bigraph(move(bigraph_1_infile), bigraph_1_filename));
        Bigraph big2 = free_all_entities(read_bigraph(move(bigraph_2_infile), bigraph_2_filename)); 

        /* Prepare and start timeout */
        params.timeout = make_shared<Timeout>(options_vars.count("timeout") ? seconds{ options_vars["timeout"].as<int>() } : 0s);

        /* Start the clock */
        params.start_time = steady_clock::now();

        InputGraph big2_encoding = big2.encode(true);

        std::vector<Bigraph> components = full_decomp(big1);
        std::vector<Bigraph> candidates;

        for(auto c : components){ 
            auto result = solve_homomorphism_problem(c.encode(false), big2_encoding, params);
            if(! result.mapping.empty())
                candidates.push_back(c);
        }
        if (candidates.size() == 0) {
            cout << "Solutions found: 0";
            return EXIT_SUCCESS;
        }

        bool finished_flag = false;
        while(! finished_flag) {
            std::vector<Bigraph> new_candidates;
            for(int i=0;i<candidates.size();i++){
                for(int j=candidates[i].largest_component_index+1;j<components.size(); j++) {
                    auto new_comp = element_compose(candidates[i], components[j]);
                    if(new_comp.has_value()) {
                        auto result = solve_homomorphism_problem(new_comp.value().encode(false), big2_encoding, params);
                        if(! result.mapping.empty())
                            new_candidates.push_back(new_comp.value());
                    }
                }
            }
            if(new_candidates.size() == 0)
                finished_flag = true;
            else
                candidates = new_candidates;
        }

        cout << "Solutions found: " << candidates.size() << '\n';
        for (auto z : candidates)
            cout << z.toString();
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
    return 0;
}