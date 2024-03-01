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

#include"flowfile_out.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"


flowfile_out::flowfile_out(lexer *p, fdm* a, ghostcell *pgc) : probenum(p->P230), eps(1.0e-10*p->DXM)
{
	//------------------------------
    initialize(p,a,pgc);
	//------------------------------
	ini_location(p,a,pgc);
	//------------------------------
    header_file_ini(p,a,pgc);
    //------------------------------
}

flowfile_out::~flowfile_out()
{
	for(n=0;n<probenum;++n)
    fileout[n].close();
}

void flowfile_out::start(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{
	int num;
	
	if(p->P15==1)
    num = filecount;

    if(p->P15==2)
    num = p->count;
    
    for(n=0;n<p->P230;++n)
    if(p->mpirank==0)
	{
        // header file
        header_file(p,a,pgc);
    }

	for(n=0;n<p->P230;++n)
    if(p->P230_x[n]>=p->originx && p->P230_x[n]<p->endx)
	{

        // filename
        filename(p,a,pgc);
        fileout[n].open(name, ios::binary);
        
        // write
        write_data(p,a,pgc);
	}
	
	++filecount;
}


