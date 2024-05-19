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

#include"cfd_state.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void cfd_state::ini_mainheader(lexer *p, fdm *a, ghostcell *pgc)
{
    ofstream mainout;
    
    // open file
	mainout.open("./REEF3D_CFD_STATE/REEF3D-CFD_State_Mainheader.r3d", ios::binary);

    // ini write
    iin=p->M10;
    mainout.write((char*)&iin, sizeof (int));

    iin=p->j_dir;
    mainout.write((char*)&iin, sizeof (int));

    iin=ie_global-is_global;
    mainout.write((char*)&iin, sizeof (int));

    iin=je_global-js_global;
    mainout.write((char*)&iin, sizeof (int));

    iin=p->gknoz;
    mainout.write((char*)&iin, sizeof (int));

    iin=file_version;
    mainout.write((char*)&iin, sizeof (int));
    
    iin=file_type;
    mainout.write((char*)&iin, sizeof (int));

    // flag: is process within P43 bounds
    for(int qn=0;qn<p->M10;++qn)
    {
    iin = flag_all[qn];
    mainout.write((char*)&iin, sizeof (int));
    }

    mainout.close();
}

void cfd_state::write_mainheader(lexer *p, fdm *a, ghostcell *pgc)
{
    ofstream mainout;
    
    // open file
	mainout.open("./REEF3D_CFD_STATE/REEF3D-CFD_State_Mainheader.r3d", ios::binary | ios::app);

    iin=p->count;
    mainout.write((char*)&iin, sizeof (int));

	ddn=p->simtime;
    mainout.write((char*)&ddn, sizeof (double));

    mainout.close();
}
