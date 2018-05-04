/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "dimacs.hh"
#include "graph.hh"
#include "graph_file_error.hh"

#include <boost/regex.hpp>
#include <fstream>

using std::getline;
using std::ifstream;
using std::stoi;
using std::string;

using boost::regex;
using boost::smatch;

auto read_dimacs(const string & filename) -> Graph
{
    Graph result{ 0 };

    ifstream infile{ filename };
    if (! infile)
        throw GraphFileError{ filename, "unable to open file" };

    string line;
    while (getline(infile, line)) {
        if (line.empty())
            continue;

        /* Lines are comments, a problem description (contains the number of
         * vertices), or an edge. */
        static const regex
            comment{ R"(c(\s.*)?)" },
            problem{ R"(p\s+(edge|col)\s+(\d+)\s+(\d+)?\s*)" },
            edge{ R"(e\s+(\d+)\s+(\d+)\s*)" };

        smatch match;
        if (regex_match(line, match, comment)) {
            /* Comment, ignore */
        }
        else if (regex_match(line, match, problem)) {
            /* Problem. Specifies the size of the graph. Must happen exactly
             * once. */
            if (0 != result.size())
                throw GraphFileError{ filename, "multiple 'p' lines encountered" };
            result.resize(stoi(match.str(2)));
        }
        else if (regex_match(line, match, edge)) {
            /* An edge. DIMACS files are 1-indexed. We assume we've already had
             * a problem line (if not our size will be 0, so we'll throw). */
            int a{ stoi(match.str(1)) }, b{ stoi(match.str(2)) };
            if (0 == a || 0 == b || a > result.size() || b > result.size())
                throw GraphFileError{ filename, "line '" + line + "' edge index out of bounds" };
            result.add_edge(a - 1, b - 1);
        }
        else
            throw GraphFileError{ filename, "cannot parse line '" + line + "'" };
    }

    if (! infile.eof())
        throw GraphFileError{ filename, "error reading file" };

    return result;
}

