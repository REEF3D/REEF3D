/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"fnpf_breaking_log.h"
#include"lexer.h"

void fnpf_breaking_log::filename(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    
    if(p->P14==0)
    {
        if(p->mpirank<9)
        sprintf(name,"REEF3D_FNPF-Breaking-Log-0000%i.r3d",p->mpirank+1);

        if(p->mpirank<99&&p->mpirank>8)
        sprintf(name,"REEF3D_FNPF-Breaking-Log-000%i.r3d",p->mpirank+1);

        if(p->mpirank<999&&p->mpirank>98)
        sprintf(name,"REEF3D_FNPF-Breaking-Log-00%i.r3d",p->mpirank+1);

        if(p->mpirank<9999&&p->mpirank>998)
        sprintf(name,"REEF3D_FNPF-Breaking-Log-0%i.r3d",p->mpirank+1);

        if(p->mpirank>9998)
        sprintf(name,"REEF3D_FNPF-Breaking-Log-%i.r3d",p->mpirank+1);
    }
    
    if(p->P14==1)
    {
        if(p->mpirank<9)
        sprintf(name,"./REEF3D_FNPF_Breaking_Log/REEF3D_FNPF-Breaking-Log-0000%i.r3d",p->mpirank+1);

        if(p->mpirank<99&&p->mpirank>8)
        sprintf(name,"./REEF3D_FNPF_Breaking_Log/REEF3D_FNPF-Breaking-Log-000%i.r3d",p->mpirank+1);

        if(p->mpirank<999&&p->mpirank>98)
        sprintf(name,"./REEF3D_FNPF_Breaking_Log/REEF3D_FNPF-Breaking-Log-00%i.r3d",p->mpirank+1);

        if(p->mpirank<9999&&p->mpirank>998)
        sprintf(name,"./REEF3D_FNPF_Breaking_Log/REEF3D_FNPF-Breaking-Log-0%i.r3d",p->mpirank+1);

        if(p->mpirank>9998)
        sprintf(name,"./REEF3D_FNPF_Breaking_Log/REEF3D_FNPF-Breaking-Log-%i.r3d",p->mpirank+1);

    }

}



