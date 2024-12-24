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

#include"fnpf_print_kinematics.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void fnpf_print_kinematics::print_ini(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    // open force_ale file
    sprintf(name,"./REEF3D_FNPF_Kinematics/REEF3D_FNPF-Kinematics-%i.dat",ID+1);
        
    fout.open(name);
    
    /*ddn=p->P88_x[ID];
	fout.write((char*)&ddn, sizeof (double));
    
    ddn=p->P88_y[ID];
	fout.write((char*)&ddn, sizeof (double));
    
    iin=p->knoz+1;
	fout.write((char*)&iin, sizeof (int));*/
    
    fout<<p->P88_x[ID]<<" "<<p->P88_y[ID]<<" "<<p->knoz+1<<endl;
}

void fnpf_print_kinematics::print_kinematics(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    // time
    ddn=p->simtime;
	fout.write((char*)&ddn, sizeof (double));
    
    // results: z ; U ; dUdt ; V ; dVdt
    for(k=0;k<p->knoz+1;++k)
    {
    ddn=p->ZSN[FIJK];
	fout.write((char*)&ddn, sizeof (double));
    
    ddn=c->U[FIJK];
	fout.write((char*)&ddn, sizeof (double));
    
    ddn=ax[k];
	fout.write((char*)&ddn, sizeof (double));
    
    ddn=c->V[FIJK];
	fout.write((char*)&ddn, sizeof (double));
    
    ddn=ay[k];
	fout.write((char*)&ddn, sizeof (double));
    
    //if(ID==0)
    //cout<<k<<" "<<p->ZSN[FIJK]<<" "<<c->U[FIJK]<<" "<<ax[k]<<endl;
    }
}