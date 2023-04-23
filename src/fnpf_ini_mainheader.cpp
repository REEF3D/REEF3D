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

#include"fnpf_ini.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"


void fnpf_ini::fnpf_restart_mainheader(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    ifstream mainhead;
    int ii1,ii2;
    
    sprintf(name,"./REEF3D_FNPF_STATE/REEF3D-FNPF_State_Mainheader.r3d");
    
	mainhead.open(name, ios::binary);
    
    // count numiter
    mainhead.read((char*)&iin, sizeof (int));
	//numprocs=iin;
    
    mainhead.read((char*)&iin, sizeof (int));
	//jdir=iin;
    
    mainhead.read((char*)&iin, sizeof (int));
	//NGx=iin;
    
    mainhead.read((char*)&iin, sizeof (int));
	//NGy=iin;
    
    mainhead.read((char*)&iin, sizeof (int));
	//NGz=iin;
    
    mainhead.read((char*)&iin, sizeof (int));
	//file_version=iin;
    
    mainhead.read((char*)&iin, sizeof (int));
	file_type=iin;
    
    mainhead.close();

    
}