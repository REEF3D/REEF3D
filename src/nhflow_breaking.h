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

#ifndef NHFLOW_BREAKING_H_
#define NHFLOW_BREAKING_H_

#include"nhflow_fsf.h"
#include"sliceint4.h"
#include"slice4.h"

using namespace std;

class ioflow;

class nhflow_breaking : public increment 
{
public:
	nhflow_breaking(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_breaking();
    
    void breaking(lexer*,fdm_nhf*,ghostcell*,slice&,slice&,slice&,double);
    
    void breaking_baquet(lexer*,fdm_nhf*,ghostcell*,slice&,slice&,slice&,double);

    // the parts of breaking_baquet: detection (bx, by, bd), flags (brkflag), viscosity (vb)
    void breaking_detect(lexer*,fdm_nhf*,ghostcell*,slice&,slice&,slice&,double);
    void breaking_flags(lexer*,fdm_nhf*,ghostcell*);
    void breaking_vb(lexer*,fdm_nhf*,ghostcell*,slice&);

    // mesh refinement (nhflow_amr): with brk_defer breaking() only keeps its arguments, the stage
    // runner calls the parts on all grids (part 0, 1, 2) and exchanges the detection and the flags
    // of the cells around the patches between sibling patches in between
    bool brk_defer = false;
    void breaking_part(lexer*,fdm_nhf*,ghostcell*,int);
    sliceint4& brk_bx() { return bx; }
    sliceint4& brk_by() { return by; }
    sliceint4& brk_bd() { return bd; }
    sliceint4& brk_flag() { return brkflag; }

    void filter(lexer*, fdm_nhf*, ghostcell*, slice&);
    
    // A 550 3: eddy viscosity of Kennedy et al. (2000)
    void breaking_kennedy(lexer*,fdm_nhf*,ghostcell*,slice&);
    
    // no detection in the relaxation zones of the wave generation and the numerical beach (ioflow::relax_boxes);
    // called by the momentum schemes before breaking(), builds the mask once
    void breaking_zones(lexer*,ioflow*);

    
private:

    double visc;
    int count_n;
    
    sliceint4 bx,by,bd,brkflag;
    sliceint4 brkzone;          // 1 in the relaxation zones
    bool brkzone_built = false;
    slice4 brkage,brkage_n;     // A 550 3: age of the breaking event, at the end of the stage and of the last time step
    slice4 brkvb;               // A 550 3: work array of the band
    slice *brk_eta = nullptr, *brk_eta_n = nullptr, *brk_WL = nullptr;
    double brk_alpha = 1.0;
};

#endif
