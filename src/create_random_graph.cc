#include <iostream>
#include <random>

#include <cxxopts.hpp>

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::mt19937;
using std::uniform_real_distribution;

auto main(int argc, char * argv[]) -> int
{
    try {
        cxxopts::Options options("Create a random graph", "Get started by using option --help");

        options.add_options("Program options")
            ("help", "Display help information");

        options.add_options("Graph options")
            ("seed", "Specify a random seed", cxxopts::value<int>())
            ("directed", "Generate a directed graph")
            ("loops", "Generate loops with this probability", cxxopts::value<double>());

        options.add_options()
            ("vertices", "Specify the number of vertices", cxxopts::value<int>())
            ("edge-probability", "Specify the edge probability", cxxopts::value<double>());

        options.parse_positional({"vertices", "edge-probability"});

        auto options_vars = options.parse(argc, argv);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            cout << options.help() << endl;
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
    catch (const cxxopts::exceptions::exception& e) {
        cerr << "Error: " << e.what() << endl;
        cerr << "Try " << argv[0] << " --help" << endl;
        return EXIT_FAILURE;
    }
    catch (const exception & e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
}
