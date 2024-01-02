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


void potentialfile_out::write_data(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    
    /*
    bedlevel SWL
    number of vertical points n
    sig_1….sig_n
    t  eta  U1……U2  W1…..W2  	*/

    i=iloc[n];
    j=0;

            
    // PRINT DATA
    /*
    ffn = float(p->simtime);
    fileout[n].write((char*)&ffn, sizeof (float));
    
    ffn = float(c->eta(i,j));
    fileout[n].write((char*)&ffn, sizeof (float));
    
    for(qn=0;qn<p->knoz;++qn)
    {
    ffn = float(c->u(i,j,k));
    fileout[n].write((char*)&ffn, sizeof (float));
    }
    
    for(qn=0;qn<p->knoz;++qn)
    {
    ffn = float(c->w(i,j,k));
    fileout[n].write((char*)&ffn, sizeof (float));
    }*/
    
}








