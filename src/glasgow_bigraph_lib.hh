#ifndef GLASGOW_BIGRAPH_LIB_H_
#define GLASGOW_BIGRAPH_LIB_H_

#include <vector>
#include <map>
#include "homomorphism.hh"

class Results {
    public:
        std::vector<VertexToVertexMapping> mapping;
        unsigned next = 0;
        int count = 0;

        Results () { mapping.clear(); };

        void clear () { mapping.clear(); next = 0; }
        bool match_found () { return !mapping.empty(); }
};

#ifdef __cplusplus
extern "C" {
#endif


// Simple interface for now, using the strings. Later we will read/write directly
void gbs_match_all(const char* pat, const char* tar);
void gbs_match_one(const char* pat, const char* tar);
int gbs_count_sols(const char* pat, const char* tar);
bool gbs_equal(const char* pat, const char* tar);

void gbs_start_pattern(int size);
void gbs_start_target(int size);

void gbs_add_node(const int i, const char* lbl, const char* name,
                  const std::vector<int> indeg, const std::vector<int> outdeg);
void gbs_add_edge(const int i, const int j);

void gbs_match_one_flat();
void gbs_match_all_flat();
int gbs_count_sols_flat();
bool gbs_equal_flat();

VertexToVertexMapping gbs_nextsol();
std::map<int, int> gbs_getEdges(const VertexToVertexMapping & mapping);
std::map<int, int> gbs_getNodes(const VertexToVertexMapping & mapping);
std::vector<std::pair<int, int>> gbs_getHyp(const VertexToVertexMapping & mapping);

#ifdef __cplusplus
}
#endif

#endif // GLASGOW_BIGRAPH_LIB_H_
