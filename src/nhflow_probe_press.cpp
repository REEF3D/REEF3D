/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"nhflow_probe_press.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_probe_press::nhflow_probe_press(lexer *p, fdm_nhf *d) : probenum(p->P64)
{

    p->Iarray(iloc,probenum);
	p->Iarray(jloc,probenum);
	p->Iarray(kloc,probenum);
	p->Iarray(flag,probenum);
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_NHFLOW_ProbePoint",0777);
	
	pout = new ofstream[probenum];
	
    if(p->mpirank==0 && probenum>0)
    {
		cout<<"probepoint_num: "<<probenum<<endl;
		// open file
		for(n=0;n<probenum;++n)
		{
		sprintf(name,"./REEF3D_NHFLOW_ProbePoint/REEF3D-NHFLOW-Probe_Press-%i.dat",n+1);
		
		pout[n].open(name);
        
        //cout<<pout[n].is_open()<<" "<<n+1<<endl;

	    pout[n]<<"Vel Probe ID:  "<<n+1<<endl<<endl;
		pout[n]<<"x_coord     y_coord     z_coord"<<endl;
		
		pout[n]<<n+1<<"\t "<<p->Xout(p->P64_x[n],p->P64_y[n])<<"\t "<<p->Yout(p->P64_x[n],p->P64_y[n])<<"\t "<<p->P64_z[n]<<endl;

		pout[n]<<endl<<endl;
		
		pout[n]<<"t \t U \t V \t W "<<endl;
		
		}
    }
	
}

nhflow_probe_press::~nhflow_probe_press()
{
    for(n=0;n<probenum;++n)
    pout[n].close();
    
    delete [] pout;
}

void nhflow_probe_press::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double xp,yp,zp;
    
    ini_location(p,d);
	
	for(n=0;n<probenum;++n)
	{
	pval=-1.0e20;
	
		if(flag[n]>0)
		{
		xp=p->P64_x[n];
		yp=p->P64_y[n];
		zp=p->P64_z[n];
    
		pdyn = p->ccipol4V(d->P, d->WL, d->bed, xp, yp, zp);
         etaval = p->ccslipol4(d->eta,xp,yp);  
         phs = (p->wd + etaval - zp)*p->W1*fabs(p->W22);
         
         pval = pdyn + phs;
         
         if(zp > p->wd + etaval)
         pval = 0.0;
		}
	
	pval=pgc->globalmax(pval);


	
	if(p->mpirank==0)
	pout[n]<<setprecision(9)<<p->simtime<<" \t "<<pval<<endl;
	}	
}

void nhflow_probe_press::ini_location(lexer *p, fdm_nhf *d)
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
    
	//kloc[n]=p->posf_sig(iloc[n],jloc[n],p->P64_z[n]);
    //check=boundcheck(p,iloc[n],jloc[n],kloc[n],0);
    
    if(iloc[n]>=0 && iloc[n]<p->knox)
    if((jloc[n]>=0 && jloc[n]<p->knoy) || p->j_dir==0)
    check=1;
    
    //cout<<p->mpirank<<" PROBE check: "<<check<<" i: "<<iloc[n]<<" j: "<<jloc[n]<<" k: "<<kloc[n]<<" ZSN: "<<p->ZSN[10+marge]<<endl;
    if(check==1)
    flag[n]=1;
    }
    
}
