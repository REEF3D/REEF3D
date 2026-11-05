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

#ifndef PATCHBC_2D_H_
#define PATCHBC_2D_H_

#include"patchBC_core.h"

using namespace std;

// patch boundary conditions ::SFLOW (cell centred velocities)
class patchBC_2D final : public patchBC_core
{
public:
	patchBC_2D(lexer*,ghostcell*);
	virtual ~patchBC_2D();
    
    void patchBC_ini(lexer*, ghostcell*) override final;
    
    // BC update ::CFD
    void patchBC_ioflow(lexer*, fdm*, ghostcell*, field&,field&,field&) override final;
    void patchBC_rkioflow(lexer*, fdm*, ghostcell*, field&,field&,field&) override final;
    void patchBC_discharge(lexer*, fdm*, ghostcell*) override final;
    void patchBC_pressure(lexer*, fdm*, ghostcell*, field&) override final;
    void patchBC_waterlevel(lexer*, fdm*, ghostcell*, field&) override final;
    
    // BC update ::SFLOW
    void patchBC_ioflow2D(lexer*, ghostcell*, slice&, slice&, slice&, slice&) override final;
    void patchBC_discharge2D(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&) override final;
    void patchBC_waterlevel2D(lexer*, fdm2D*, ghostcell*, slice&) override final;
    
private:
    // ghost cell q (1..3) of the face (i,j,cs)
    void ghost(int cs, int q, int &ii, int &jj);
    
    // code of the patch face matching a staggered slice entry (same side, same tangential index,
    // normal index within one cell), 0 if none
    int staggered_flag(int ii, int jj, int cs);
};

#endif
