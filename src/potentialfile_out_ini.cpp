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
#include"potentialfile_out.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void potentialfile_out::initialize(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    filecount=0;
    
    if(p->mpirank==0)
	mkdir("./REEF3D_PotentialFile",0777);
	
	if(p->mpirank==0 && p->P240>0)
	cout<<"PotentialFile: "<<probenum<<endl;

	fileout = new ofstream[p->P240];
    
    p->Iarray(iloc,p->P240);
    
    for(n=0;n<p->P240;++n)
    iloc[n] = p->posf_i(p->P240_x[n]);
    
    
    // filename
        filename(p,c,pgc);
        fileout[n].open(name, ios::binary);
        
    // start file
    i=iloc[n];
    j=0;
    
    for(n=0;n<p->P240;++n)
    if(p->P240_x[n]>=p->originx && p->P240_x[n]<p->endx)
	{

        ffn = float(c->bed(i,j));
        fileout[n].write((char*)&ffn, sizeof (float));
        
        iin=p->knoz;
        fileout[n].write((char*)&iin, sizeof (int));
        
        for(qn=0;qn<p->knoz+1;++qn)
        {
        ffn = float(p->ZN[KP]);
        fileout[n].write((char*)&ffn, sizeof (float));
        }
    }
    
}

void potentialfile_out::ini_location(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{

}

