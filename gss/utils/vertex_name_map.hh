#ifndef GLASGOW_SUBGRAPH_SOLVER_VERTEX_NAME_MAP_HH
#define GLASGOW_SUBGRAPH_SOLVER_VERTEX_NAME_MAP_HH

#include <unordered_map>
#include <string>

struct BiMap {
    std::unordered_map<int, std::string> id_to_name;
    std::unordered_map<std::string, int> name_to_id;

    void insert(const int id, const std::string& name) {
        erase(id);
        erase(name);

        id_to_name[id] = name;
        name_to_id[name] = id;
    }

    void erase(const int id) {
        const auto it = id_to_name.find(id);
        if (it != id_to_name.end()) {
            name_to_id.erase(it->second);
            id_to_name.erase(it);
        }
    }

    void erase(const std::string& name) {
        const auto it = name_to_id.find(name);
        if (it != name_to_id.end()) {
            id_to_name.erase(it->second);
            name_to_id.erase(it);
        }
    }

    auto find_left(const int id) const {
        return id_to_name.find(id);
    }

    auto find_right(const std::string& name) const {
        return name_to_id.find(name);
    }
};

#endif // GLASGOW_SUBGRAPH_SOLVER_VERTEX_NAME_MAP_HH
