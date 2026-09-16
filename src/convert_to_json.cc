#include <gss/formats/json_graph.hh>
#include <gss/formats/read_file_format.hh>

#include <iostream>

#include <cxxopts.hpp>

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::string;

auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Convert graph to the gss-graph JSON format", "Get started by using option --help");

        options.add_options("Program options") //
            ("help", "Display help information") //
            ("format", "Specify input file format (auto, lad, vertexlabelledlad, labelledlad, directedlad, dimacs, csv, json)", cxxopts::value<string>());

        options.add_options() //
            ("graph-file", "Specify the graph file", cxxopts::value<string>());

        options.parse_positional({"graph-file"});

        auto options_vars = options.parse(argc, argv);

        if (options_vars.count("help")) {
            cout << options.help() << endl;
            return EXIT_SUCCESS;
        }

        if (! options_vars.count("graph-file")) {
            cout << "Usage: " << argv[0] << " [options] graph-file" << endl;
            return EXIT_FAILURE;
        }

        string format_name = options_vars.count("format") ? options_vars["format"].as<string>() : "auto";
        auto graph = read_file_format(format_name, options_vars["graph-file"].as<string>());

        // Unlike convert_to_lad, nothing has to be refused here: the format carries
        // directedness, names and both kinds of label, which is the point of it.
        write_json_graph(cout, graph);

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
