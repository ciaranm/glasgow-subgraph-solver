/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/read_file_format.hh"
#include "formats/dimacs.hh"
#include "formats/lad.hh"

#include <boost/regex.hpp>
#include <fstream>
#include <sstream>
#include <vector>

using std::ifstream;
using std::ios;
using std::move;
using std::stoi;
using std::string;
using std::stringstream;
using std::to_string;
using std::vector;

using boost::regex;
using boost::smatch;

auto detect_format(ifstream & infile, const string & filename) -> string
{
    string line;
    if (! getline(infile, line) || line.empty())
        throw GraphFileError{ filename, "unable to read file to detect file format" };

    static const regex
        dimacs_comment{ R"(c(\s.*)?)" },
        dimacs_problem{ R"(p\s+(edge|col)\s+(\d+)\s+(\d+)?\s*)" },
        lad_header{ R"(\d+)" },
        lad_zero_labelled_line{ R"(0 \d+)" },
        lad_zero_unlabelled_line{ R"(0)" },
        lad_line{ R"([1-9]\d*\s+(\d+\s+)*\d+\s*)" };

    smatch match;
    if (regex_match(line, match, dimacs_comment)) {
        while (regex_match(line, match, dimacs_comment)) {
            // looks like a DIMACS comment, ignore
            if (! getline(infile, line) || line.empty())
                throw GraphFileError{ filename, "unable to auto-detect file format (entirely c lines?)" };
        }
        if (! regex_match(line, match, dimacs_problem))
            throw GraphFileError{ filename, "unable to auto-detect file format (c line not followed by a p line?)" };
        return "dimacs";
    }
    else if (regex_match(line, match, dimacs_problem))
        return "dimacs";
    else if (regex_match(line, match, lad_header)) {
        if ("0" == line)
            return "lad";

        // got to figure out whether we're labelled or not
        if (! getline(infile, line) || line.empty())
            throw GraphFileError{ filename, "unable to auto-detect file format (number followed by nothing)" };
        if (regex_match(line, match, lad_zero_labelled_line))
            return "labelledlad";
        else if (regex_match(line, match, lad_zero_unlabelled_line))
            return "lad";
        else if (regex_match(line, match, lad_line)) {
            stringstream line_stream{ line };
            vector<string> words;
            string word;
            while (line_stream >> word)
                words.emplace_back(word);

            if (words.size() < 2)
                throw GraphFileError{ filename, "unable to auto-detect file format (not enough words in a lad line)" };

            unsigned items = stoi(words.at(1));
            if (words.size() == items + 1)
                return "lad";
            else if (words.size() == (2 * items) + 2)
                return "labelledlad";
            else
                throw GraphFileError{ filename, "unable to auto-detect file format (looks like lad, but got " +
                    to_string(words.size()) + " items on a line with degree " + to_string(items) + ")" };
        }
        else
            throw GraphFileError{ filename, "unable to auto-detect file format (looks like lad, but no edge line found)" };
    }

    throw GraphFileError{ filename, "unable to auto-detect file format (no recognisable header found)" };
}

auto read_file_format(const string & format, const string & filename) -> InputGraph
{
    ifstream infile{ filename };
    if (! infile)
        throw GraphFileError{ filename, "unable to open file" };

    auto actual_format = format;
    if (actual_format == "auto") {
        actual_format = detect_format(infile, filename);
        infile.clear();
        if (! infile.seekg(0, ios::beg))
            throw GraphFileError{ filename, "unable to seek on input file (try specifying file format explicitly)" };
    }

    if (actual_format == "dimacs")
        return read_dimacs(move(infile), filename);
    else if (actual_format == "lad")
        return read_lad(move(infile), filename);
    else if (actual_format == "labelledlad")
        return read_labelled_lad(move(infile), filename);
    else
        throw GraphFileError{ filename, "Unknown file format '" + format + "'" };
}

