/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"patchBC_interface.h"

using namespace std;

class patchBC_2D : public patchBC_interface, public increment
{
public:
	patchBC_2D(lexer*,ghostcell*);
	virtual ~patchBC_2D();
    
    virtual void patchBC_ini(lexer *p, ghostcell *pgc);
    
    // BC update ::CFD
    virtual void patchBC_ioflow(lexer*, fdm*, ghostcell*, field&,field&,field&);
    virtual void patchBC_rkioflow(lexer*, fdm*, ghostcell*, field&,field&,field&);
    virtual void patchBC_discharge(lexer*, fdm*, ghostcell*);
    virtual void patchBC_pressure(lexer*, fdm*, ghostcell*, field&);
    virtual void patchBC_waterlevel(lexer*, fdm*, ghostcell*, field&);
    
    // BC update ::SFLOW
    virtual void patchBC_ioflow2D(lexer*, ghostcell*, slice&, slice&, slice&, slice&);
    virtual void patchBC_rkioflow2D(lexer*, ghostcell*, slice&, slice&, slice&, slice&);
    virtual void patchBC_discharge2D(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&);
    virtual void patchBC_pressure2D(lexer*, ghostcell*, slice&);
    virtual void patchBC_pressure2D_ugrad(lexer*, fdm2D*, slice&,slice&);
    virtual void patchBC_pressure2D_vgrad(lexer*, fdm2D*, slice&, slice&);
    virtual void patchBC_waterlevel2D(lexer*, fdm2D*,  ghostcell*, slice&);

    virtual void patchBC_loop2D(lexer*, fdm2D*, int&, int&, int&, int&);
        
private:
     // ini
    void patchBC_gcb_count(lexer *p, ghostcell *pgc);
    void patchBC_IDcount(lexer *p, ghostcell *pgc);
    void patchBC_fillobj(lexer *p, ghostcell *pgc);
    
    void patchBC_hydrograph_Q_read(lexer *p, ghostcell *pgc,int,int);
    double patchBC_hydrograph_Q_ipol(lexer *p, ghostcell *pgc,int,int);
    
    void patchBC_hydrograph_FSF_read(lexer *p, ghostcell *pgc,int,int);
    double patchBC_hydrograph_FSF_ipol(lexer *p, ghostcell *pgc,int,int);
    
    
    int q,n,qn,qq,count,ID_count;
    int istart,iend,jstart,jend,kstart,kend;
    
    int *inflow_ID;
    int *outflow_ID;
    
    int geo_count;
    int *ID_array;
    
};

#endif
