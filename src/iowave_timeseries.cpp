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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void iowave::timeseries(lexer *p, ghostcell* pgc)
{
    ofstream pout;
    char name[100];
    
    double time0=p->simtime;
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_Log-Wave",0777);
    
    if(p->P58==0)
    {
    p->P58=1;
        
    p->Darray(p->P58_x,p->P58);
	p->Darray(p->P58_y,p->P58);
    p->Darray(p->P58_T,p->P58);
    
    p->P58_x[0] = 0.0;
    p->P58_y[0] = 0.0;
    p->P58_T[0] = 12800.0;
    }
    
	

    for(int n=0; n<p->P58; ++n)
    {
		if(p->P14==0)
		sprintf(name,"REEF3D-Wave-Timeseries-%i.dat",n+1);
		
		if(p->P14==1)
		sprintf(name,"./REEF3D_Log-Wave/REEF3D-Wave-Timeseries-%i.dat",n+1);
		
		pout.open(name);

	    pout<<"Wave Timeseries ID:  "<<n<<endl<<endl;
		pout<<"x_coord     y_coord     T"<<endl;
		
		pout<<"\t "<<p->P58_x[n]<<"\t "<<p->P58_y[n]<<"\t "<<p->P58_T[n]<<endl<<endl;
        
        pout<<"t \t eta"<<endl<<endl;
        
        p->simtime=0.0;
        do
        {
        pout<<p->simtime<<" "<<wave_eta(p,pgc,p->P58_x[n],p->P58_y[n])<<endl;
            
        p->simtime+=0.1;
        
        }while(p->simtime<=p->P58_T[n]);


    pout.close();
    }
	

    p->simtime=time0;
    
}