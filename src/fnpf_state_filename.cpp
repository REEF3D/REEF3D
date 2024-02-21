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

#include"fnpf_state.h"
#include"lexer.h"

void fnpf_state::filename_single(lexer *p, fdm_fnpf *c, ghostcell *pgc, int num)
{
    if(p->P14==0)
    {
    sprintf(name,"REEF3D_FNPF-State-%08i-%05i.r3d",num,p->mpirank+1);
    }
    
    if(p->P14==1)
    {
    sprintf(name,"./REEF3D_FNPF_STATE/REEF3D_FNPF-State-%08i-%05i.r3d",num,p->mpirank+1);
    }
}

void fnpf_state::filename_continuous(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->P14==0)
    {
    sprintf(name,"REEF3D_FNPF-State-%05i.r3d",p->mpirank+1);
    }
    
    if(p->P14==1)
    {
    sprintf(name,"./REEF3D_FNPF_STATE/REEF3D_FNPF-State-%05i.r3d",p->mpirank+1);
    }
}

void fnpf_state::filename_header(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->P14==0)
    {
	sprintf(name,"REEF3D-FNPF-State-Header-%05i.r3d",p->mpirank+1);
    }
    
    if(p->P14==1)
    {
	sprintf(name,"./REEF3D_FNPF_STATE/REEF3D-FNPF-State-Header-%05i.r3d",p->mpirank+1);
    }
}



