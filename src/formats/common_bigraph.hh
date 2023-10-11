#include <cstdint>
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

class Entity {
    public:
        int id;
        std::vector<Entity> children;
        std::set<int> regions;
        std::set<int> sites;

        string control;
        int arity;
        Entity *parent = NULL;

        Entity(int i, string ctrl, int ar);
        Entity();

        auto copy() -> Entity;
};

class Bigraph {
    public:
        std::set<int> regions;
        std::set<int> sites;
        std::vector<Entity> entities;
        std::vector<std::vector<bool>> reachability;
        int largest_component_index;

        Bigraph();

        auto decomp(std::vector<int> values) -> Bigraph;

        auto toString() const -> string;

        auto encode(bool target) const -> InputGraph;

        auto copy() -> Bigraph;
};

auto read_bigraph(istream && infile, const string &) -> Bigraph;

auto full_decomp(Bigraph big) -> std::vector<Bigraph>;

auto element_compose(Bigraph a, Bigraph b) -> std::optional<Bigraph>;

