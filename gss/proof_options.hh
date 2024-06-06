#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_PROOF_OPTIONS_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_PROOF_OPTIONS_HH

#include <string>

namespace gss
{
    struct ProofOptions
    {
        std::string opb_file;
        std::string log_file;
        bool recover_encoding = false;
        bool super_extra_verbose = false;
    };
}

#endif
