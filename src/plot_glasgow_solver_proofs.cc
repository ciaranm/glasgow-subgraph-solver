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

auto read_instance_file(const string & instance_result_file_name, const string & prefix,
        map<string, string> & keys, bool & aborted, bool & excluded) -> void
{
    ifstream result_file{ instance_result_file_name };
    if (! result_file) {
        cerr << "Error reading " << instance_result_file_name << endl;
        throw 0;
    }
    string key_line;
    while (getline(result_file, key_line)) {
        auto pos = key_line.find('=');
        if (string::npos == pos) {
            cerr << "Couldn't parse '" << key_line << "' in '" << instance_result_file_name << endl;
            throw 0;
        }

        auto k = key_line.substr(0, pos);
        k.erase(k.find_last_not_of(" ") + 1, k.length());
        auto v = key_line.substr(pos + 1);
        v.erase(0, v.find_first_not_of(" "));
        keys.emplace(prefix + k, v);

        if (k == "status") {
            if (v == "true") {
            }
            else if (v == "false") {
            }
            else if (v == "aborted" || v == "tooslow") {
                aborted = true;
            }
            else if (v == "excluded") {
                excluded = true;
            }
            else {
                cerr << "Couldn't parse status value '" << v << "' in " << instance_result_file_name << endl;
                throw 0;
            }
        }
    }
}

auto read_veripb_file(const string & veripb_file_name, const string & prefix,
        map<string, string> & keys, bool & aborted, bool & failed) -> void
{
    ifstream result_file{ veripb_file_name };
    if (! result_file) {
        cerr << "Error reading " << veripb_file_name << endl;
        throw 0;
    }
    string key_line;

    while (getline(result_file, key_line)) {
        if (0 == key_line.compare(0, 22, "INFO:root:total time: ")) {
            key_line.erase(0, 22);
            keys[prefix + "runtime"] = key_line;
        }
        else if (0 == key_line.compare(0, 36, "INFO:root:time in OPBParser::parse: ")) {
            key_line.erase(0, 36);
            keys[prefix + "opbparsetime"] = key_line;
        }
        else if (0 == key_line.compare(0, 41, "INFO:root:time in RuleParserBase::parse: ")) {
            key_line.erase(0, 41);
            keys[prefix + "ruleparsetime"] = key_line;
        }
        else if (0 == key_line.compare(0, 38, "INFO:root:time in propEngine::attach: ")) {
            key_line.erase(0, 38);
            keys[prefix + "attachtime"] = key_line;
        }
        else if (0 == key_line.compare(0, 40, "INFO:root:time in LoadFormula::compute: ")) {
            key_line.erase(0, 40);
            keys[prefix + "loadformulatime"] = key_line;
        }
        else if (0 == key_line.compare(0, 50, "INFO:root:time in ReversePolishNotation::compute: ")) {
            key_line.erase(0, 50);
            keys[prefix + "rpntime"] = key_line;
        }
        else if (0 == key_line.compare(0, 56, "INFO:root:time in ConstraintImpliesGetImplied::compute: ")) {
            key_line.erase(0, 56);
            keys[prefix + "impliedtime"] = key_line;
        }
        else if (0 == key_line.compare(0, 51, "INFO:root:time in ReverseUnitPropagation::compute: ")) {
            key_line.erase(0, 51);
            keys[prefix + "ruptime"] = key_line;
        }
        else if (0 == key_line.compare(0, 37, "INFO:root:time in Solution::compute: ")) {
            key_line.erase(0, 37);
            keys[prefix + "solutiontime"] = key_line;
        }
        else if (0 == key_line.compare(0, 22, "used database memory: ")) {
            key_line.erase(0, 22);
            keys[prefix + "used_db_memory"] = key_line;
        }
        else if (0 == key_line.compare(0, 28, "cumulative database memory: ")) {
            key_line.erase(0, 28);
            keys[prefix + "cumulative_db_memory"] = key_line;
        }
        else if (0 == key_line.compare(0, 30, "maixmal used database memory: ")) {
            key_line.erase(0, 30);
            keys[prefix + "maximal_db_memory"] = key_line;
        }
        else if (0 == key_line.compare(0, 7, "WARNING")) {
            cerr << "Warning in " << veripb_file_name << endl;
            throw 0;
        }
        else if (0 == key_line.compare(0, 5, "ERROR")) {
            cerr << "Error in " << veripb_file_name << endl;
            throw 0;
        }
        else if (key_line == "Verification succeeded.") {
            keys[prefix + "status"] = "success";
        }
        else if (key_line == "Verification failed.") {
            keys[prefix + "status"] = "failed";
            failed = true;
        }
        else if (key_line == "return code exited failure" && ! failed) {
            keys[prefix + "status"] = "aborted";
            aborted = true;
        }
    }
}

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                         "Display help information");

        po::options_description all_options{ "All options" };
        all_options.add_options()
            ("instances-file",    po::value<string>(),          "Specify the instances file (first column specifies output name)")
            ("results-directory", po::value<string>(),          "Directory for results without proofs")
            ("proofs-directory",  po::value<string>(),          "Directory for results with proofs")
            ;

        all_options.add(display_options);

        po::positional_options_description positional_options;
        positional_options
            .add("instances-file", 1)
            .add("results-directory", 1)
            .add("proofs-directory", 1)
            ;

        po::variables_map options_vars;
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(positional_options)
                .run(), options_vars);
        po::notify(options_vars);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << "Usage: " << argv[0] << " [options] instances-file results-directory proofs-directory" << endl;
            cout << endl;
            cout << display_options << endl;
            return EXIT_SUCCESS;
        }

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("instances-file") || ! options_vars.count("results-directory")
                || ! options_vars.count("proofs-directory")) {
            cout << "Usage: " << argv[0] << " [options] instances-file results-directory proofs-directory" << endl;
            return EXIT_FAILURE;
        }

        auto results_dir = options_vars["results-directory"].as<string>();
        auto proofs_dir = options_vars["proofs-directory"].as<string>();

        ifstream instances{ options_vars["instances-file"].as<string>() };
        if (! instances) {
            cerr << "Error reading instances file" << endl;
            return EXIT_FAILURE;
        }

        ofstream out("proofs.data");

        out << "instance status pattern_vertices pattern_directed_edges target_vertices target_directed_edges "
            << "runtime search_nodes solution_count proof_status proof_runtime proof_opbsize proof_logsize "
            << "verify_status verify_runtime "
            << "verify_opbparsetime verify_attachtime verify_ruleparsetime verify_loadformulatime verify_rpntime "
            << "verify_impliedtime verify_ruptime verify_solutiontime verify_used_db_memory verify_cumulative_db_memory "
            << "verify_maximal_db_memory"
            << endl;
        if (! out) {
            cerr << "Error writing output file" << endl;
            return EXIT_FAILURE;
        }

        vector<string> show_keys = { "instance", "status", "pattern_vertices", "pattern_directed_edges",
            "target_vertices", "target_directed_edges", "runtime", "nodes", "solution_count",
            "proof_status", "proof_runtime", "proof_opbsize", "proof_logsize", "verify_status", "verify_runtime",
            "verify_opbparsetime", "verify_attachtime", "verify_ruleparsetime", "verify_loadformulatime", "verify_rpntime",
            "verify_impliedtime", "verify_ruptime", "verify_solutiontime", "verify_used_db_memory",
            "verify_cumulative_db_memory", "verify_maximal_db_memory" };

        string instances_line;
        while (getline(instances, instances_line)) {
            stringstream instances_line_s{ instances_line };
            string instance_name;
            if (! (instances_line_s >> instance_name))
                continue;

            string instance_result_file_name{ results_dir + "/" + instance_name + ".out" };

            map<string, string> keys;
            bool aborted = false, excluded = false, proof_aborted = false, proof_excluded = false;
            read_instance_file(instance_result_file_name, "", keys, aborted, excluded);

            keys["instance"] = instance_name;
            if (aborted) {
                keys["status"] = "aborted";
                keys["nodes"] = "NaN";
                keys["runtime"] = "NaN";
                keys["solution_count"] = "NaN";
                keys["proof_status"] = "excluded";
                keys["proof_runtime"] = "NaN";
                keys["proof_opbsize"] = "0";
                keys["proof_logsize"] = "0";
                keys["verify_status"] = "excluded";
                for (auto & k : show_keys)
                    if (0 == k.compare(0, 7, "verify_"))
                        keys.emplace(k, "NaN");
            }
            else {
                keys["status"] = "success";

                string proofs_result_file_name{ proofs_dir + "/" + instance_name + ".out" };
                read_instance_file(proofs_result_file_name, "proof_", keys, proof_aborted, proof_excluded);
                if (proof_aborted || proof_excluded) {
                    keys["proof_status"] = (proof_aborted ? "aborted" : "excluded");
                    keys["proof_runtime"] = "NaN";
                    keys["proof_opbsize"] = "0";
                    keys["proof_logsize"] = "0";
                    keys["verify_status"] = "excluded";
                    for (auto & k : show_keys)
                        if (0 == k.compare(0, 7, "verify_"))
                            keys.emplace(k, "NaN");
                }
                else {
                    keys["proof_status"] = "success";

                    string ls_file_name{ proofs_dir + "/" + instance_name + ".ls" };
                    ifstream ls{ ls_file_name };
                    string line;
                    while (getline(ls, line)) {
                        stringstream line_s{ line };
                        string ignore, size_word, filename_word;
                        line_s >> ignore >> ignore >> ignore >> ignore >> size_word >> ignore >> ignore >> ignore >> filename_word;
                        if (0 == filename_word.compare(filename_word.length() - 4, 4, ".veripb", 0, 4))
                            keys["proof_logsize"] = size_word;
                        else if (0 == filename_word.compare(filename_word.length() - 4, 4, ".opb", 0, 4))
                            keys["proof_opbsize"] = size_word;
                    }

                    bool veripb_aborted = false, veripb_failed = false;
                    string veripb_file_name{ proofs_dir + "/" + instance_name + ".veripb" };
                    read_veripb_file(veripb_file_name, "verify_", keys, veripb_aborted, veripb_failed);
                    if (veripb_aborted || veripb_failed) {
                        keys["verify_runtime"] = "NaN";
                        for (auto & k : show_keys)
                            if (0 == k.compare(0, 7, "verify_"))
                                keys.emplace(k, "NaN");
                    }
                    else if (! keys.count("verify_solutiontime"))
                        keys.emplace("verify_solutiontime", "0");
                }
            }

            bool first = true;
            for (auto & k : show_keys) {
                if (! keys.count(k)) {
                    cerr << "Missing key " << k << " in " << instance_result_file_name << endl;
                    return EXIT_FAILURE;
                }
                if (! first)
                    out << " ";
                first = false;

                for (auto p = keys.find(k)->second.find(' ') ; p != string::npos ; p = keys.find(k)->second.find(' '))
                    keys.find(k)->second[p] = '_';
                out << keys.find(k)->second;
            }

            out << endl;
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

