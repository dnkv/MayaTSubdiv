#ifndef TCC_H
#define TCC_H

#include "hds.h"


bool compute_TCC_stencils(HDS &hds, uint64_t fI, std::v_double &Fs, std::v_double &Es, std::v_double &Vs);
bool compute_Tknot_insertion(HDS &hds, std::v_double &VV, std::v_double &EV, std::v_size_t &Ktip);
void tccsubdiv(HDS &hds, size_t nSubdivs);

#endif