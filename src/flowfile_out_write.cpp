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

#include"flowfile_out.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void flowfile_out::write_data(lexer *p, fdm *a, ghostcell *pgc)
{
    
    // PRINT DATA
            i=iloc[n];
            j=0;
        
        
            for(k=0;k<Nk;++k)
            {
            ffn = float(p->XP[IP]);
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
            for(k=0;k<Nk;++k)
            {
            ffn = float(p->YP[JP]);
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
            for(k=0;k<Nk;++k)
            {
            ffn = float(p->ZP[KP]);
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
            for(k=0;k<Nk;++k)
            {
            if(flag[n][k]==1)
            ffn = float(a->u(i,j,k));//0.5*(a->u(i,j,k) + a->u(i+1,j,k);
            
            if(flag[n][k]==0)
            ffn=0.0;
            
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
            for(k=0;k<Nk;++k)
            {
            if(flag[n][k]==1)
            ffn = float(a->v(i,j,k));//0.5*(a->v(i,j,k) + a->v(i,j+1,k);
            
            if(flag[n][k]==0)
            ffn=0.0;
            
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
            for(k=0;k<Nk;++k)
            {
            if(flag[n][k]==1)
            ffn = float(a->w(i,j,k));//0.5*(a->w(i,j,k) + a->w(i,j,k+1));
            
            if(flag[n][k]==0)
            ffn=0.0;
            
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
            for(k=0;k<Nk;++k)
            {
            if(flag[n][k]==1)
            ffn = float(a->press(i,j,k));
            
            if(flag[n][k]==0)
            ffn=0.0;
            
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
            for(k=0;k<Nk;++k)
            {
            if(flag[n][k]==1)
            ffn = float(a->phi(i,j,k));
            
            if(flag[n][k]==0)
            ffn=0.0;
            
            fileout[n].write((char*)&ffn, sizeof (float));
            }
            
		fileout[n].close();
		
    
}








