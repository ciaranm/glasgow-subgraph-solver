#include <gss/formats/csv.hh>
#include <gss/formats/dimacs.hh>
#include <gss/formats/lad.hh>
#include <gss/formats/read_file_format.hh>
#include <gss/formats/vfmcs.hh>

#include <fstream>
#include <regex>
#include <sstream>
#include <vector>

using std::ifstream;
using std::ios;
using std::move;
using std::regex;
using std::smatch;
using std::stoi;
using std::string;
using std::stringstream;
using std::to_string;
using std::vector;

auto detect_format(ifstream & infile, const string & filename) -> string
{
    string line;
    if (! getline(infile, line) || line.empty())
        throw GraphFileError{filename, "unable to read file to detect file format", true};

    static const regex
        dimacs_comment{R"(c(\s.*)?)"},
        dimacs_problem{R"(p\s+(edge|col)\s+(\d+)\s+(\d+)?\s*)"},
        lad_header{R"(\d+)"},
        lad_zero_labelled_line{R"(0 \d+)"},
        lad_zero_unlabelled_line{R"(0)"},
        lad_line{R"(\d+\s+(\d+\s+)*\d+\s*)"},
        csv_problem{R"(\S+[,>]\S+)"};

    smatch match;
    if (regex_match(line, match, dimacs_comment)) {
        while (regex_match(line, match, dimacs_comment)) {
            // looks like a DIMACS comment, ignore
            if (! getline(infile, line) || line.empty())
                throw GraphFileError{filename, "unable to auto-detect file format (entirely c lines?)", true};
        }
        if (! regex_match(line, match, dimacs_problem))
            throw GraphFileError{filename, "unable to auto-detect file format (c line not followed by a p line?)", true};
        return "dimacs";
    }
    else if (regex_match(line, match, dimacs_problem))
        return "dimacs";
    else if (regex_match(line, match, lad_header)) {
        if ("0" == line)
            return "lad";

        // got to figure out whether we're labelled or not
        if (! getline(infile, line) || line.empty())
            throw GraphFileError{filename, "unable to auto-detect file format (number followed by nothing)", true};
        if (regex_match(line, match, lad_zero_labelled_line))
            return "labelledlad";
        else if (regex_match(line, match, lad_zero_unlabelled_line))
            return "lad";
        else if (regex_match(line, match, lad_line)) {
            stringstream line_stream{line};
            vector<string> words;
            string word;
            while (line_stream >> word)
                words.emplace_back(word);

            if (words.size() < 2)
                throw GraphFileError{filename, "unable to auto-detect file format (not enough words in a lad line)", true};

            unsigned items_if_lad = stoi(words.at(0)), items_if_labelledlad = stoi(words.at(1));
            if (words.size() == items_if_lad + 1 && ! (words.size() == (2 * items_if_labelledlad) + 2))
                return "lad";
            else if (! (words.size() == items_if_lad + 1) && (words.size() == (2 * items_if_labelledlad) + 2))
                return "labelledlad";
            else if ((words.size() == items_if_lad + 1) && (words.size() == (2 * items_if_labelledlad) + 2))
                throw GraphFileError{filename, "unable to auto-detect between lad or labelledlad", true};
            else
                throw GraphFileError{filename, "unable to auto-detect file format (looks like lad, but got " + to_string(words.size()) + " items on a line with first two elements " + to_string(items_if_lad) + " and " + to_string(items_if_labelledlad) + ")", true};
        }
        else
            throw GraphFileError{filename, "unable to auto-detect file format (looks like lad, but no edge line found)", true};
    }
    else if (regex_match(line, match, csv_problem))
        return "csv";

    throw GraphFileError{filename, "unable to auto-detect file format (no recognisable header found)", true};
}

auto read_file_format(const string & format, const string & filename) -> InputGraph
{
    ifstream infile{filename};
    if (! infile)
        throw GraphFileError{filename, "unable to open file", false};

    auto actual_format = format;
    if (actual_format == "auto") {
        actual_format = detect_format(infile, filename);
        infile.clear();
        if (! infile.seekg(0, ios::beg))
            throw GraphFileError{filename, "unable to seek on input file (try specifying file format explicitly)", true};
    }

    if (actual_format == "dimacs")
        return read_dimacs(move(infile), filename);
    else if (actual_format == "lad")
        return read_lad(move(infile), filename);
    else if (actual_format == "directedlad")
        return read_directed_lad(move(infile), filename);
    else if (actual_format == "labelledlad")
        return read_labelled_lad(move(infile), filename);
    else if (actual_format == "vertexlabelledlad")
        return read_vertex_labelled_lad(move(infile), filename);
    else if (actual_format == "csv")
        return read_csv(move(infile), filename);
    else if (actual_format == "vfmcs")
        return read_unlabelled_undirected_vfmcs(move(infile), filename);
    else if (actual_format == "vfmcsv")
        return read_vertex_labelled_undirected_vfmcs(move(infile), filename);
    else if (actual_format == "vfmcsvd")
        return read_vertex_labelled_directed_vfmcs(move(infile), filename);
    else if (0 == actual_format.compare(0, 8, "csvname:"))
        return read_csv_name(move(infile), filename, actual_format.substr(8));
    else
        throw GraphFileError{filename, "Unknown file format '" + format + "'", true};
}
