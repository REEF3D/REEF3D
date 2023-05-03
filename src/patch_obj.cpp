/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"patch_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

patch_obj::patch_obj(lexer *p, int ID_ini) 
{
    ID = ID_ini;
    
    gcb_count=0;
    
    
    // BC options ini
    
    Q_flag=0;
    Q=0.0;
    Uq=0.0;
    
    velocity_flag=0;
    velocity=0.0;
    
    pressure_flag=0;
    pio_flag=0;
    pressure=0.0;
    
    waterlevel_flag=0;
    waterlevel=0.0;
    
    Uio_flag=0;
    Uio=0.0;
    
    velcomp_flag=0;
    U=V=W=0.0;
    
    flowangle_flag=0;
    alpha=0.0;
    cosalpha=0.0;
    sinalpha=1.0;
    
    flownormal_flag=0;
    Nx=Ny=Nz=0.0;
    
    hydroQ_flag=0;
    hydroFSF_flag=0;
    
    counter=0;
    
    
    // measurements
    
    Q0 = 0.0;
    U0 = 0.0;
    A0 = 0.0;
    h0 = 0.0;
    
    // gcbflags
    gcb_uflag=1;
    gcb_pressflag=1;
    gcb_phiflag=1;
    
    gcb_flag=111;
    
}

patch_obj::~patch_obj()
{
}

void patch_obj::patch_obj_ini(lexer *p, ghostcell *pgc)
{
}

void patch_obj::patch_obj_gcb_generate(lexer *p, ghostcell *pgc)
{
    p->Iarray(gcb,gcb_count,4);
}