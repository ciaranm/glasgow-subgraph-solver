/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/read_file_format.hh"

#include <boost/program_options.hpp>

#include <iostream>
#include <unistd.h>

namespace po = boost::program_options;

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::string;

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                         "Display help information")
            ("format",             po::value<string>(),      "Specify input file format (auto, lad, vertexlabelledlad, labelledlad, dimacs)")
            ;

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
            cout << "Usage: " << argv[0] << " [options] file" << endl;
            cout << endl;
            cout << display_options << endl;
            return EXIT_SUCCESS;
        }

        /* No algorithm or no input file specified? Show a message and exit. */
        if (! options_vars.count("graph-file")) {
            cout << "Usage: " << argv[0] << " [options] file" << endl;
            return EXIT_FAILURE;
        }

        /* Read in the graphs */
        string default_format_name = options_vars.count("format") ? options_vars["format"].as<string>() : "auto";
        auto graph = read_file_format(default_format_name, options_vars["graph-file"].as<string>());

        cout << graph.size() << endl;
        for (int v = 0 ; v < graph.size() ; ++v) {
            cout << graph.degree(v);
            for (int w = 0 ; w < graph.size() ; ++w)
                if (graph.adjacent(v, w))
                    cout << " " << w;
            cout << endl;
        }

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


