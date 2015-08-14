#ifndef TCC_MAX_H
#define TCC_MAX_H

#include "hds.h"

namespace TCC_MAX
{
    bool compute_TCC_max_stencils(HDS &hds, uint64_t fI, std::v_double &Fs, std::v_double &Es, std::v_double &Vs);
    bool compute_Tknot_insertion(HDS &hds, std::v_double &VV, std::v_double &EV, std::v_size_t &Ktip);
    void subdivide(HDS &hds, size_t nSubdivs);
    void linear_subdivide(HDS &hds, size_t nSubdivs);
}
#endif