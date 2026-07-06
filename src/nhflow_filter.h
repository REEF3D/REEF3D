/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_FILTER_H_
#define NHFLOW_FILTER_H_

#include"increment.h"
#include<vector>

class lexer;
class ghostcell;
class slice;

class nhflow_filter : public increment
{
public:
    // poly_order : polynomial degree of the local fit (2 is a good default)
    // half_width : window half-width in cells (2 -> 5x5 window)
    nhflow_filter(lexer*);
    virtual ~nhflow_filter();

    // Savitzky-Golay filter smooth f in place; iterations>1 applies the filter repeatedly
    void filter_sg(lexer*, ghostcell*, slice& f);

    // Savitzky-Golay filter smooth "in" into a separate output slice "out" (single pass)
    void filter_sg(lexer*, ghostcell*, slice& in, slice& out);
    
    // linear filter with predictor corrector step
    void filter_pc(lexer*, ghostcell*, slice&);

    // dump the convolution kernel (rank 0) for inspection
    void print_kernel(lexer*);

private:
    void build_coeffs(lexer*);
    void invert(std::vector<std::vector<double> > M,
                std::vector<std::vector<double> >& Minv, int n);

    int order;                    // polynomial degree
    int hw;                       // window half-width (cells)
    int nw;                       // window edge length = 2*hw+1
    std::vector<double> kernel;   // size nw*nw, ordered dv outer, du inner
    std::vector<double> buffer;   // scratch for in-place filtering
};

#endif
