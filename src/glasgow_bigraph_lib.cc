#include "glasgow_bigraph_lib.hh"

#include "formats/bigraph.hh"
#include "homomorphism.hh"
#include "lackey.hh"
#include "symmetries.hh"
#include "proof.hh"
#include "restarts.hh"

#include <chrono>
#include <cstdlib>
#include <string>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <regex>

#include <unistd.h>

using std::boolalpha;
using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::function;
using std::ifstream;
using std::localtime;
using std::make_pair;
using std::make_shared;
using std::make_unique;
using std::put_time;
using std::string;
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::operator""s;
using std::chrono::seconds;
using std::chrono::steady_clock;
using std::chrono::system_clock;

// Global variables to maintain results/graphs, cleared as required
static Results res;
static bool isPat = true;
static InputGraph patG;
static InputGraph tarG;

// Compile regex ahead of time (and globally accessible) for performance
const std::regex linkL { R"(L(\d+)_.*)" };
const std::regex linkAny { R"((L|C)(\d+)_.*)" };

VertexToVertexMapping gbs_nextsol() {
    if (res.mapping.empty() || res.mapping.size() <= res.next) {
        return {};
    }


    auto r = res.mapping[res.next];
    res.next++;
    return r;
}

std::map<int, int> gbs_getEdges(const VertexToVertexMapping & mapping) {
    std::map<int, int> res;
    for (auto v : mapping) {
        if(patG.vertex_name(v.first).find("C_LINK") != string::npos) {
            int k = std::stoi(patG.vertex_name(v.first).substr(7));
            int val = std::stoi(tarG.vertex_name(v.second).substr(7));
            res[k] = val;
        }
    }
    return res;
}


std::map<int, int> gbs_getNodes(const VertexToVertexMapping & mapping) {
    std::map<int, int> res;
    for (auto v : mapping) {
        if(patG.vertex_name(v.first).find("C_LINK") == string::npos
           && patG.vertex_label(v.first) != "LINK") {
            int k = std::stoi(patG.vertex_name(v.first));
            int val = std::stoi(tarG.vertex_name(v.second));
            res[k] = val;
        }
    }
    return res;
}

std::vector<std::pair<int,int>> gbs_getHyp(const VertexToVertexMapping & mapping) {
    std::vector<std::pair<int, int>> res;

    std::smatch match;
    for (auto v : mapping) {
        if(patG.vertex_label(v.first) == "LINK") {
            std::string str = patG.vertex_name(v.first);
            if (regex_match(str, match, linkL)) {
                int l1 = stoi(match.str(1));

                str = tarG.vertex_name(v.second);
                if (regex_match(str, match, linkAny)) {
                    int l2 = stoi(match.str(2));
                    res.emplace_back(std::make_pair(l1,l2));
                }
            }
        }
    }
    return res;
}

// Assumes target/pattern setup before calling
void gbs_match_one() {
    res.clear();

    if (tarG.size() == 0 || patG.size() == 0) {
        std::cerr << "Target or Pattern not available in match\n";
    }

    HomomorphismParams params;
    params.injectivity = Injectivity::Injective;
    params.induced = false;
    params.bigraph = true;
    params.count_solutions = false;
    params.no_supplementals = true;
    params.use_bigraph_projection_nogoods = false;

    /* Prepare and start timeout */
    params.restarts_schedule = make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier);
    params.timeout = make_shared<Timeout>(0s);

    auto result = solve_homomorphism_problem(patG, tarG, params);

    if(!result.mapping.empty() && params.bigraph) {
        res.mapping.push_back(result.mapping);
    }
}

void gbs_match_all() {
    res.clear();

    if (tarG.size() == 0 || patG.size() == 0) {
        std::cerr << "Target or Pattern not available in match\n";
    }

    HomomorphismParams params;
    params.injectivity = Injectivity::Injective;
    params.induced = false;
    params.bigraph = true;
    params.count_solutions = true;
    params.no_supplementals = true;
    params.use_bigraph_projection_nogoods = false;

    /* Prepare and start timeout */
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.timeout = make_shared<Timeout>(0s);

    params.enumerate_callback = [&](auto m) { res.mapping.push_back(m); };

    auto result = solve_homomorphism_problem(patG, tarG, params);
}

int gbs_count_sols() {
    res.clear();

    HomomorphismParams params;
    params.injectivity = Injectivity::Injective;
    params.induced = false;
    params.bigraph = true;
    params.count_solutions = true;
    params.no_supplementals = true;
    params.use_bigraph_projection_nogoods = false;

    /* Prepare and start timeout */
    params.restarts_schedule = make_unique<NoRestartsSchedule>();
    params.timeout = make_shared<Timeout>(0s);

    auto result = solve_homomorphism_problem(patG, tarG, params);

    return (int) result.solution_count;
}

auto gbs_equal() -> bool {
    res.clear();

    HomomorphismParams params;
    params.injectivity = Injectivity::Injective;
    params.induced = false;
    params.bigraph = true;
    params.equality_check = true;
    // params.bigraph_equal = true;
    // Get all in case the first solution doesn't do the hyperedges properly
    // TODO constriain this in search
    params.count_solutions = true;
    params.no_supplementals = true;

    params.restarts_schedule = make_unique<NoRestartsSchedule>();

    if(params.bigraph) {
        params.enumerate_callback = [&](auto m) {res.mapping.push_back(m);};
    }

    /* Prepare and start timeout */
    params.timeout = make_shared<Timeout>(0s);

    auto result = solve_homomorphism_problem(patG, tarG, params);

    if (res.mapping.empty()) {
        return false;
    }

    // Just need one to be equal
    std::smatch match;
    for (auto m : res.mapping) {
        bool failure = false;
        for (auto v : m) {
            if (patG.vertex_name(v.first).find("ROOT") != string::npos) {
                int l = std::stoi(patG.vertex_name(v.first).substr(4));
                int r = std::stoi(tarG.vertex_name(v.second).substr(4));
                if (l != r) { failure = true; break; } // Roots are not identity
            }

            if(patG.vertex_label(v.first) == "LINK") {
                int l1, l2;
                std::string str = patG.vertex_name(v.first);
                if (regex_match(str, match, linkL)) {
                    l1 = stoi(match.str(1));

                    str = tarG.vertex_name(v.second);
                    if (regex_match(str, match, linkAny)) {
                        l2 = stoi(match.str(2));
                        if (l1 != l2) {
                            failure = true;
                            break; } // Hyperedge not identity
                    }
                }
            }
        }

        if (!failure) {
            return true;
        }

    }

    return false;
}


void gbs_add_node(const int i, const char* lbl, const char* name,
                  const std::vector<int> indeg, const std::vector<int> outdeg) {
    auto &ig = isPat? patG : tarG;
    ig.set_vertex_label(i, lbl);
    ig.set_vertex_name(i, name);

    if (indeg.size() > 0) {
        ig.set_child_of_root(i);
        for (const auto & j : indeg) {
            ig.add_pattern_root_edge(j, i);
        }
    }
    if (outdeg.size() > 0) {
        ig.set_parent_of_site(i);
        for (const auto & j : outdeg) {
            ig.add_pattern_site_edge(i, j);
        }
    }
}

void gbs_add_edge(const int i, const int j) {
    auto &ig = isPat? patG : tarG;
    ig.add_directed_edge(i,j,"dir");
}

void gbs_start_pattern(int size) {
    isPat = true;
    patG = InputGraph(size, true, true, true);
}

void gbs_start_target(int size) {
    isPat = false;
    tarG = InputGraph(size, true, true, true);
}
