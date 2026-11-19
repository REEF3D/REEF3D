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

#ifndef SFLOW_RECONSTRUCT_WENO_H_
#define SFLOW_RECONSTRUCT_WENO_H_

#include"sflow_reconstruct.h"
#include"weno_nug_func.h"
#include"slice4.h"
#include<vector>

class lexer;
class ghostcell;
class fdm2D;
class slice;
class patchBC_interface;

using namespace std;

class sflow_reconstruct_weno final : public sflow_reconstruct, public weno_nug_func
{
public:
	sflow_reconstruct_weno(lexer*,patchBC_interface*);
	sflow_reconstruct_weno(lexer*,patchBC_interface*,int);   // own WENO tables (SFLOW AMR patch)
	virtual ~sflow_reconstruct_weno();

    void reconstruct_x(lexer*,ghostcell*,fdm2D*,slice&,slice&,slice&) override final;
    void reconstruct_y(lexer*,ghostcell*,fdm2D*,slice&,slice&,slice&) override final;
    void reconstruct_WL(lexer*,ghostcell*,fdm2D*) override final;
    void prepare(lexer*) override final;
    

private:
    inline void iqmin(lexer*, slice&);
    inline void iqmax(lexer*, slice&);
    inline void jqmin(lexer*, slice&);
    inline void jqmax(lexer*, slice&);

    void ytab_fill(lexer*);

    int qq;

    // one row of left and right face values; the WENO coefficients of the y-faces, flat (ytab_fill)
    std::vector<double> rowL,rowR,ytab;
    
    // dry cells in the stencil of a wet cell (prepare), per direction: 0 WENO, 1 the cell value, 2 the central stencil
    std::vector<unsigned char> drx,dry;
    bool drs;
    
    patchBC_interface *pBC;
};

#endif
