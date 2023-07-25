#include <iostream>
#include <random>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::mt19937;
using std::uniform_real_distribution;

auto main(int argc, char * argv[]) -> int
{
    try {
        po::options_description display_options{"Program options"};
        display_options.add_options()("help", "Display help information");

        po::options_description graph_options{"Graph options"};
        graph_options.add_options()                             //
            ("seed", po::value<int>(), "Specify a random seed") //
            ("directed", "Generate a directed graph")           //
            ("loops", po::value<double>(), "Generate loops with this probability");
        display_options.add(graph_options);

        po::options_description all_options{"All options"};
        all_options.add_options()                                            //
            ("vertices", po::value<int>(), "Specify the number of vertices") //
            ("edge-probability", po::value<double>(), "Specify the edge probability");

        all_options.add(display_options);

        po::positional_options_description positional_options;
        positional_options
            .add("vertices", 1)
            .add("edge-probability", 1);

        po::variables_map options_vars;
        po::store(po::command_line_parser(argc, argv)
                      .options(all_options)
                      .positional(positional_options)
                      .run(),
            options_vars);
        po::notify(options_vars);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << "Usage: " << argv[0] << " [options] number-of-vertices edge-probability" << endl;
            cout << endl;
            cout << display_options << endl;
            return EXIT_SUCCESS;
        }

        if (! options_vars.count("vertices") || ! options_vars.count("edge-probability")) {
            cout << "Usage: " << argv[0] << " [options] number-of-vertices edge-probability" << endl;
            return EXIT_FAILURE;
        }

        int seed = 0;
        if (options_vars.count("seed"))
            seed = options_vars["seed"].as<int>();

        int vertices = options_vars["vertices"].as<int>();
        double density = options_vars["edge-probability"].as<double>();
        double loops = options_vars.count("loops") ? options_vars["loops"].as<double>() : 0;

        bool directed = options_vars.count("directed");

        mt19937 rand;
        rand.seed(seed);
        uniform_real_distribution<double> dist(0.0, 1.0);

        for (int v = 0; v < vertices; ++v) {
            cout << "v" << v << "," << endl;
            if (loops > dist(rand))
                cout << "v" << v << ","
                     << "v" << v << endl;
            for (int w = (directed ? 0 : v + 1); w < vertices; ++w)
                if (v != w && density > dist(rand))
                    cout << "v" << v << (directed ? ">" : ",") << "v" << w << endl;
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
