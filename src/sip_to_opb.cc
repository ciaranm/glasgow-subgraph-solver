/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/read_file_format.hh"

#include <boost/program_options.hpp>

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

namespace po = boost::program_options;

using std::cerr;
using std::copy;
using std::cout;
using std::endl;
using std::exception;
using std::greater;
using std::istreambuf_iterator;
using std::map;
using std::ostreambuf_iterator;
using std::pair;
using std::sort;
using std::string;
using std::string_view;
using std::stringstream;
using std::vector;

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                         "Display help information")
            ("format",             po::value<string>(),      "Specify input file format (auto, lad, vertexlabelledlad, labelledlad, dimacs)")
            ("degree",                                       "Apply degree preprocessing")
            ("nds",                                          "Also apply neighbourhood degree sequence preprocessing")
            ;

        po::options_description all_options{ "All options" };
        all_options.add_options()
            ("pattern-file", "Specify the pattern file")
            ("target-file",  "Specify the target file")
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

        bool degree = options_vars.count("degree");
        bool nds = options_vars.count("nds");

        /* Read in the graphs */
        string default_format_name = options_vars.count("format") ? options_vars["format"].as<string>() : "auto";
        auto pattern = read_file_format(default_format_name, options_vars["pattern-file"].as<string>());
        auto target = read_file_format(default_format_name, options_vars["target-file"].as<string>());

        vector<int> pattern_degrees(pattern.size());
        vector<int> target_degrees(target.size());

        vector<vector<int> > pattern_ndss(pattern.size());
        vector<vector<int> > target_ndss(target.size());

        vector<vector<int> > domains(pattern.size());
        vector<vector<int> > antidomains(target.size());
        map<pair<int, int>, int> numberings;

        if (degree) {
            pattern.for_each_edge([&] (int f, int t, string_view) {
                if (f != t) {
                    pattern_degrees[f]++;
                    pattern_ndss[f].push_back(t);
                }
            });

            for (int v = 0 ; v < pattern.size() ; ++v) {
                for (auto & v : pattern_ndss[v])
                    v = pattern_degrees[v];
                sort(pattern_ndss[v].begin(), pattern_ndss[v].end(), greater<int>());
            }

            target.for_each_edge([&] (int f, int t, string_view) {
                if (f != t) {
                    target_degrees[f]++;
                    target_ndss[f].push_back(t);
                }
            });

            for (int v = 0 ; v < target.size() ; ++v) {
                for (auto & v : target_ndss[v])
                    v = target_degrees[v];
                sort(target_ndss[v].begin(), target_ndss[v].end(), greater<int>());
            }
        }

        for (int i = 0 ; i < pattern.size() ; ++i) {
            for (int j = 0 ; j < target.size() ; ++j) {
                bool ok = true;

                if (pattern.has_vertex_labels() && pattern.vertex_label(i) != target.vertex_label(j))
                    ok = false;
                else if (pattern.adjacent(i, i) && ! target.adjacent(j, j))
                    ok = false;
                else if (degree && (pattern_degrees[i] > target_degrees[j]))
                    ok = false;
                else if (degree && nds) {
                    for (int x = 0 ; ok && x < pattern_degrees[i] ; ++x)
                        if (pattern_ndss[i][x] > target_ndss[j][x])
                            ok = false;
                }

                if (ok) {
                    domains.at(i).push_back(j);
                    antidomains.at(j).push_back(i);
                    numberings.emplace(pair{ i, j }, numberings.size() + 1);
                }
            }

            if (domains.at(i).empty()) {
                cerr << "Inconsistency detected during model creation: domain " << i << " empty" << endl;
                return EXIT_FAILURE;
            }
        }

        stringstream model;
        int nb_constraints = 0;

        // domains
        for (int i = 0 ; i < pattern.size() ; ++i) {
            for (auto & v : domains[i])
                model << "1 x" << numberings[{ i, v }] << " ";
            model << ">= 1 ;" << endl;
            for (auto & v : domains[i])
                model << "-1 x" << numberings[{ i, v }] << " ";
            model << ">= -1 ;" << endl;
            nb_constraints += 2;
        }

        // injectivity
        for (int i = 0 ; i < target.size() ; ++i) {
            if (antidomains.at(i).size() > 0) {
                for (auto & v : antidomains[i])
                    model << "-1 x" << numberings[{ v, i }] << " ";
                model << ">= -1 ;" << endl;
                ++nb_constraints;
            }
        }

        // edges
        for (int i = 0 ; i < pattern.size() ; ++i) {
            for (auto & v : domains[i]) {
                for (int j = 0 ; j < pattern.size() ; ++j) {
                    if (pattern.adjacent(i, j)) {
                        model << "1 ~x" << numberings[{ i, v }];
                        for (auto & k : domains[j])
                            if (k != v && target.adjacent(v, k))
                                model << " 1 x" << numberings[{ j, k }];
                        model << " >= 1 ;" << endl;
                        ++nb_constraints;
                    }
                }
            }
        }

        cout << "* #variable= " << numberings.size() << " #constraint= " << nb_constraints << endl;
        copy(istreambuf_iterator<char>{ model }, istreambuf_iterator<char>{}, ostreambuf_iterator<char>{ cout });

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


