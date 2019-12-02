/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <boost/program_options.hpp>

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace po = boost::program_options;

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::ifstream;
using std::make_unique;
using std::map;
using std::ofstream;
using std::string;
using std::stringstream;
using std::unique_ptr;
using std::vector;

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                         "Display help information");

        po::options_description all_options{ "All options" };
        all_options.add_options()
            ("instances-file",    po::value<string>(),          "Specify the instances file (first column specifies output name)")
            ("results-directory", po::value<vector<string> >(), "Directories which contain results")
            ;

        all_options.add(display_options);

        po::positional_options_description positional_options;
        positional_options
            .add("instances-file", 1)
            .add("results-directory", -1)
            ;

        po::variables_map options_vars;
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(positional_options)
                .run(), options_vars);
        po::notify(options_vars);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << "Usage: " << argv[0] << " [options] instances-file results-directory..." << endl;
            cout << endl;
            cout << display_options << endl;
            return EXIT_SUCCESS;
        }

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("instances-file") || ! options_vars.count("results-directory")) {
            cout << "Usage: " << argv[0] << " [options] instances-file results-directory..." << endl;
            return EXIT_FAILURE;
        }

        auto results_dirs = options_vars["results-directory"].as<vector<string> >();

        ifstream instances{ options_vars["instances-file"].as<string>() };
        if (! instances) {
            cerr << "Error reading instances file" << endl;
            return EXIT_FAILURE;
        }

        map<string, unique_ptr<ofstream> > output_files;
        output_files.emplace("runtime", make_unique<ofstream>("runtimes.data"));
        output_files.emplace("nodes", make_unique<ofstream>("nodes.data"));
        output_files.emplace("status", make_unique<ofstream>("statuses.data"));

        for (auto & [ _, f ] : output_files) {
            *f << "instance";
            for (auto & d : results_dirs)
                *f << " " << d;
            *f << endl;

            if (! *f) {
                cerr << "Error writing output file" << endl;
                return EXIT_FAILURE;
            }
        }

        string instances_line;
        while (getline(instances, instances_line)) {
            stringstream instances_line_s{ instances_line };
            string instance_name;
            if (! (instances_line_s >> instance_name))
                continue;

            for (auto & [ _, f ] : output_files)
                *f << instance_name;

            for (auto & d : results_dirs) {
                string instance_result_file_name{ d + "/" + instance_name + ".out" };
                ifstream result_file{ instance_result_file_name };
                if (! result_file) {
                    cerr << "Error reading " << instance_result_file_name << endl;
                    return EXIT_FAILURE;
                }

                map<string, string> keys;
                bool aborted = false;
                string key_line;
                while (getline(result_file, key_line)) {
                    auto pos = key_line.find('=');
                    if (string::npos == pos) {
                        cerr << "Couldn't parse '" << key_line << "' in " << instance_result_file_name << endl;
                        return EXIT_FAILURE;
                    }

                    auto k = key_line.substr(0, pos);
                    k.erase(k.find_last_not_of(" ") + 1, k.length());
                    auto v = key_line.substr(pos + 1);
                    v.erase(0, v.find_first_not_of(" "));
                    keys.emplace(k, v);

                    if (k == "status") {
                        if (v == "true") {
                        }
                        else if (v == "false") {
                        }
                        else if (v == "aborted") {
                            aborted = true;
                        }
                        else {
                            cerr << "Couldn't parse status value '" << v << " in " << instance_result_file_name << endl;
                            return EXIT_FAILURE;
                        }
                    }
                }

                for (auto & [ k, f ] : output_files) {
                    if (! keys.count(k)) {
                        cerr << "Missing key " << k << " in " << instance_result_file_name << endl;
                        return EXIT_FAILURE;
                    }

                    if (aborted && (k == "runtime" || k == "nodes"))
                        *f << " " << "NaN";
                    else
                        *f << " " << keys.find(k)->second;
                }
            }

            for (auto & [ _, f ] : output_files)
                *f << endl;
        }

        return EXIT_SUCCESS;
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


