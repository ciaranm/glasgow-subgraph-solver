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

struct MCB_solution {
    Bigraph bigraph;
    InputGraph encoding;
    HomomorphismResult result;
};

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                         "Display help information")
            ("timeout",            po::value<int>(),         "Abort after this many seconds");

        po::options_description problem_options{ "Problem options" };
        problem_options.add_options()
            ("print-all-solutions",                          "Print out every solution, rather than one")
            ("minimal-LTS",                                  "Find the miminal context for a grounded agent in a LTS");
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
        bool lts = options_vars.count("minimal-LTS");

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

        Bigraph big1;
        Bigraph big2;

        if(!lts) {
            big1 = free_hyperedges(free_sites(free_regions(read_bigraph(move(bigraph_1_infile), bigraph_1_filename))));
            big2 = free_hyperedges(free_sites(free_regions(read_bigraph(move(bigraph_2_infile), bigraph_2_filename)))); 
        }
        else {
            big1 = free_hyperedges(read_bigraph(move(bigraph_1_infile), bigraph_1_filename));
            big2 = free_hyperedges(read_bigraph(move(bigraph_2_infile), bigraph_2_filename));
        }

        /* Prepare and start timeout */
        params.timeout = make_shared<Timeout>(options_vars.count("timeout") ? seconds{ options_vars["timeout"].as<int>() } : 0s);

        /* Start the clock */
        params.start_time = steady_clock::now();

        InputGraph big2_encoding = big2.encode(true, false);

        std::vector<Bigraph> components = full_decomp(big1);

        std::vector<MCB_solution> candidates;
        std::vector<MCB_solution> prev_solutions;
        bool full_solution_exists = false;
        int matcher_calls = 0;

        for(auto c : components){ 
            if(! lts) {
                if(c.entities.size() > 0) {
                    auto encoding = c.encode(false, false);
                    auto result = solve_homomorphism_problem(encoding, big2_encoding, params);
                    matcher_calls++;
                    if(! result.mapping.empty()) {
                        full_solution_exists = true;
                        candidates.push_back(MCB_solution{c, encoding, result});
                    }
                }
            }
            else {
                if(c.entities.size() > 0 && c.entities[0].is_leaf) {
                    auto encoding = free_regions(c).encode(false, false);
                    auto result = solve_homomorphism_problem(encoding, big2_encoding, params);
                    matcher_calls++;
                    if(! result.mapping.empty()) {
                        auto new_encoding = c.encode(false, true);
                        auto new_result = solve_homomorphism_problem(new_encoding, big2_encoding, params);
                        matcher_calls++;
                        if(! new_result.mapping.empty()) {
                            full_solution_exists = true;
                        }
                        candidates.push_back(MCB_solution{c, new_encoding, new_result}); 
                    }
                }
            }
        }

        if (candidates.size() == 0) {
            cout << "Solutions found: 0\n";
            cout << "Matcher calls: " << matcher_calls << "\n";
            return EXIT_SUCCESS;
        }

        if (full_solution_exists)
            prev_solutions = candidates;
        
        bool finished_flag = false;
        vector<int> nogood_list;
        while(! finished_flag) {
            std::vector<MCB_solution> new_candidates;
            full_solution_exists = false;
            for(unsigned int i=0;i<candidates.size();i++){
                if(! lts) {
                    for(unsigned int j=candidates[i].bigraph.largest_component_index+1;j<components.size(); j++) {
                        auto new_comp = element_compose(candidates[i].bigraph, components[j], lts);
                        if(new_comp.has_value()) {
                            bool nogood_flag = false;
                            for (int n : nogood_list) {
                                if ((new_comp.value().nogood_id | n) == new_comp.value().nogood_id) {
                                    nogood_flag = true;
                                    break;
                                }
                            }
                            if (nogood_flag) continue;
                            InputGraph big1_encoding = new_comp.value().encode(false, false);
                            auto result = solve_homomorphism_problem(big1_encoding, big2_encoding, params);
                            matcher_calls++;
                            if(! result.mapping.empty()) {
                                new_candidates.push_back(MCB_solution{new_comp.value(), big1_encoding, result});
                                full_solution_exists = true;
                            }
                            else
                                nogood_list.push_back(new_comp.value().nogood_id);
                        }
                    }
                }
                else{
                    for(unsigned int j=0;j<components.size(); j++) {
                        auto new_comp = element_compose(candidates[i].bigraph, components[j], lts);
                        if(new_comp.has_value()) {
                            bool nogood_flag = false;
                            for (int n : nogood_list) {
                                if(new_comp.value().nogood_id == n) {
                                    nogood_flag = true;
                                    break;
                                }
                            }
                            if(nogood_flag) continue;
                            nogood_list.push_back(new_comp.value().nogood_id);
                            InputGraph big1_encoding = free_regions(new_comp.value()).encode(false, false);
                            auto result = solve_homomorphism_problem(big1_encoding, big2_encoding, params);
                            matcher_calls++;
                            if(! result.mapping.empty()) {
                                InputGraph big1_extra = new_comp.value().encode(false, true);
                                auto new_result = solve_homomorphism_problem(big1_extra, big2_encoding, params);
                                matcher_calls++;
                                if(! new_result.mapping.empty()) {
                                    full_solution_exists = true;
                                }
                                new_candidates.push_back(MCB_solution{new_comp.value(), big1_extra, new_result}); 
                            }
                        }
                    }
                }
            }
            if(new_candidates.size() == 0)
                finished_flag = true; 
            else {
                if(full_solution_exists)
                    prev_solutions = new_candidates;
                candidates = new_candidates;
                if(lts) nogood_list = {};
            }
        }

        string output = "";
        int count = 0;
        bool print_flag = false;
        for (auto z : prev_solutions)
            if(! z.result.mapping.empty()) {
                count++;
                if(!print_flag) {
                    output += "mapping = {";
                    bool lazy_flag = false;
                    for (auto v : z.result.mapping) {
                        if(z.encoding.vertex_name(v.first).find("C_LINK") != string::npos) break;
                        if(z.encoding.vertex_label(v.first) == "LINK") continue;
                        if(lazy_flag) output += ",";
                        lazy_flag = true;
                        if(z.encoding.vertex_name(v.first).find("ROOT") != string::npos){
                            output += "(" + z.encoding.vertex_name(v.first) + ", " + big2_encoding.vertex_name(v.second) + ")";
                        }
                        else{
                            output += "(" + std::to_string(z.bigraph.entities[std::stoi(z.encoding.vertex_name(v.first))].id);
                            output += ", " + big2_encoding.vertex_name(v.second) + ")";
                        }
                    }
                    output += "} -- {";
                    
                    lazy_flag = false;
                    for (auto v : z.result.mapping) { 
                        if(z.encoding.vertex_name(v.first).find("C_LINK") == string::npos) continue;
                        if(lazy_flag) cout << ",";
                        lazy_flag = true;
                        output += "(" + std::to_string(z.bigraph.closures[std::stoi(z.encoding.vertex_name(v.first).substr(7))].id);
                        output += ", " + big2_encoding.vertex_name(v.second).substr(7) + ")";
                    }              
                    output += "}\n---\n";
                    output += remove_redundant_sites(z.bigraph).toString() + "---\n";
                    if(!options_vars.count("print-all-solutions"))
                        print_flag = true;
                }
            }
        output = "Solutions found: " + std::to_string(count) + "\n" + "Matcher calls: " + std::to_string(matcher_calls) + "\n" + output;
        cout << output;
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