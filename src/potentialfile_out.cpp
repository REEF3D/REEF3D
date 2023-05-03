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

#include"potentialfile_out.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"turbulence.h"

potentialfile_out::potentialfile_out(lexer *p, fdm_fnpf *c, ghostcell *pgc) : probenum(p->P230), eps(1.0e-10*p->DXM)
{
	//------------------------------
    initialize(p,c,pgc);
	//------------------------------
	ini_location(p,c,pgc);
	//------------------------------
    header_file_ini(p,c,pgc);
    //------------------------------
}

potentialfile_out::~potentialfile_out()
{
	for(n=0;n<probenum;++n)
    fileout[n].close();
}

void potentialfile_out::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, turbulence *pturb)
{
	int num;
	
	if(p->P15==1)
    num = filecount;

    if(p->P15==2)
    num = p->count;


	for(n=0;n<p->P230;++n)
    if(p->P230_x[n]>=p->originx && p->P230_x[n]<p->endx)
	{

        
        
        // write
        write_data(p,c,pgc);
	}
	
	++filecount;
}


