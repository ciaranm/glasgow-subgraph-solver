#include <cstdint>
#include "homomorphism.hh"
#include "formats/input_graph.hh"
#include <map>
#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>
#include <bitset>
#include <set>

using std::istream;
using std::string;
using std::set;

class Closure {
    public:
        int id;
        int port_count;
        string name;
        std::vector<int> adjacencies;
};

class Entity {
    public:
        int id;
        std::vector<int> child_indices;
        std::set<int> regions;
        std::set<int> sites;

        string control;
        int arity;
        int parent_index = -1;

        bool is_leaf;

        Entity(int i, string ctrl, int ar);
        Entity();

        auto copy() -> Entity;
};

class Bigraph {
    public:
        int nogood_id = 0;
        std::set<int> regions;
        std::set<int> sites;
        std::vector<Entity> entities;
        int original_size = 0;

        std::vector<std::pair<string, std::vector<int>>> hyperedges;
        std::vector<Closure> closures;

        std::vector<std::vector<bool>> reachability;
        int largest_component_index;

        std::vector<std::pair<bool, std::vector<int>>> mappings;
        
        Bigraph();

        auto decomp(std::vector<int> values) -> Bigraph;

        auto toString() const -> string;

        auto encode(bool target, bool special_lts_case) const -> InputGraph;

        auto copy() -> Bigraph;
};


auto read_bigraph(istream && infile, const string &) -> Bigraph;

auto full_decomp(Bigraph big) -> std::vector<Bigraph>;

auto split_regions(Bigraph a) -> Bigraph;

auto free_hyperedges(Bigraph a) -> Bigraph;

auto make_place_RPO(Bigraph big1, Bigraph big2, Bigraph solution, std::vector<int> mapping) -> Bigraph;

auto element_compose(Bigraph a, Bigraph b, Bigraph c, bool lts) -> std::pair<int, Bigraph>;

