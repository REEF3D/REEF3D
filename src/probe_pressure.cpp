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

#include"probe_pressure.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<sys/stat.h>
#include<sys/types.h>

probe_pressure::probe_pressure(lexer *p, fdm* a, ghostcell *pgc) : probenum(p->P64)
{
    p->Iarray(iloc,probenum);
	p->Iarray(jloc,probenum);
	p->Iarray(kloc,probenum);
	p->Iarray(flag,probenum);
	
    //cout<<p->mpirank<<" pressure probepoint_num: "<<probenum<<endl;
    
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_PressureProbe",0777);
	
	pout = new ofstream[probenum+1];
	
    if(p->mpirank==0 && probenum>0)
    {
		// open file
		for(n=0;n<probenum;++n)
		{     
		if(p->P14==0)
		sprintf(name,"REEF3D-CFD-Probe-Pressure-%i.dat",n+1);
		
		if(p->P14==1)
		sprintf(name,"./REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-%i.dat",n+1);
        
        
		pout[n].open(name,std::fstream::out);
        
        //cout<<pout[n].is_open()<<" "<<n+1<<endl;
		}
        
        
		for(n=0;n<probenum;++n)
		{
	    pout[n]<<"Point Probe ID:  "<<n+1<<endl<<endl;
		pout[n]<<"x_coord     y_coord     z_coord"<<endl;
		
		pout[n]<<n+1<<"\t "<<p->P64_x[n]<<"\t "<<p->P64_y[n]<<"\t "<<p->P64_z[n]<<endl;

		pout[n]<<endl<<endl;
		}
    }
	
    ini_location(p,a,pgc);
}

probe_pressure::~probe_pressure()
{
	for(n=0;n<probenum;++n)
    pout[n].close();
}

void probe_pressure::start(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{
	double xp,yp,zp;
	
	for(n=0;n<probenum;++n)
	{
	uval=vval=wval=pval=kval=eval=edval=-1.0e20;
	
		if(flag[n]>0)
		{
		xp=p->P64_x[n];
		yp=p->P64_y[n];
		zp=p->P64_z[n];
		
		pval = p->ccipol4_a(a->press, xp, yp, zp) - p->pressgage;
		}
	
	pval=pgc->globalmax(pval);
	
	if(p->mpirank==0)
	pout[n]<<setprecision(9)<<p->simtime<<" \t "<<pval<<endl;
	}		
}

void probe_pressure::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void probe_pressure::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check;

    for(n=0;n<probenum;++n)
    {
    check=0;
    
    iloc[n]=p->posc_i(p->P64_x[n]);
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n]=p->posc_j(p->P64_y[n]);
    
	kloc[n]=p->posc_k(p->P64_z[n]);

    check=boundcheck(p,a,iloc[n],jloc[n],kloc[n],0);

    if(check==1)
    flag[n]=1;
    }
}


