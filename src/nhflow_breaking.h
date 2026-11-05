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
#include<functional>

using namespace std;

class ioflow;
class breaking_crest3d;

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
    slice* brk_wl() { return brk_WL; }      // the water level of the deferred breaking (A 552: filtered)

    void filter(lexer*, fdm_nhf*, ghostcell*, slice&);
    
    // A 550 3: eddy viscosity of Kennedy et al. (2000)
    void breaking_kennedy(lexer*,fdm_nhf*,ghostcell*,slice&);
    // its parts (multi-block grid: nhflow_amr::breaking_level0 runs them on all blocks): 0 the ages of the last step,
    // 1 the columns and the band width kn_band, 2 the ghost cells of a pass of the band, 3 the pass qn, 4 the flags
    void kennedy_part(lexer*,fdm_nhf*,ghostcell*,slice&,int,int);
    void kennedy_run(lexer *p, fdm_nhf *d, ghostcell *pgc, int part, int qn) { kennedy_part(p,d,pgc,*brk_WL,part,qn); }
    int kn_band = 0;            // the band width of the grid in cells (part 1)
    int kn_count = 0;           // its breaking columns (part 1)
    slice4& brk_age_n() { return brkage_n; }
    slice4& brk_k_n() { return brkk_n; }
    slice4& brk_r_n() { return brkr_n; }
    
    // A 550 4: onset from crest tracking B = u_s/C (Barthelemy et al. 2018), the seed of breaking_kennedy
    void breaking_barthelemy(lexer*,fdm_nhf*,ghostcell*,slice&);
    
    // no detection in the relaxation zones of the wave generation and the numerical beach (ioflow::relax_boxes);
    // called by the momentum schemes before breaking(), builds the mask once
    void breaking_zones(lexer*,ioflow*);
    // a block of a multi-block grid (nhflow_amr): the zone columns of its edges, in(i, j) for every column with i, j set
    void breaking_zone_mask(lexer*,const std::function<void(const std::function<void(int,int)>&)>&);

    
private:

    double visc;
    int count_n;
    
    sliceint4 bx,by,bd,brkflag;
    sliceint4 brkzone;          // 1 in the relaxation zones
    bool brkzone_built = false;
    slice4 brkage,brkage_n;     // A 550 3: age of the breaking event, at the end of the stage and of the last time step
    slice4 brkvb;               // A 550 3: work array of the band
    breaking_crest3d *bt = nullptr;   // A 550 4
    slice4 bt_us, bt_vs;        // A 550 4: horizontal velocity at the surface
    slice4 bt_vbs;              // A 550 4: viscosity of the tracking (not used)
    sliceint4 brkseed;          // A 550 4: seeded columns
    slice4 brkkseed, brkk, brkk_n;   // A 550 4: wavenumber of the seed, of the event, at the end of the last step
    slice4 brkr, brkr_n;        // A 550 4: initial threshold of the event / c, at the end of the stage and of the last step
    int bt_count = -1;
    slice *brk_eta = nullptr, *brk_eta_n = nullptr, *brk_WL = nullptr;
    double brk_alpha = 1.0;
};

#endif
