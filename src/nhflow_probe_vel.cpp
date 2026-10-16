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

#include"probe_core.h"
#include"nhflow_probe_vel.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_probe_vel::nhflow_probe_vel(lexer *p, fdm_nhf *d) : probenum(p->P65)
{

    p->Iarray(iloc,probenum);
	p->Iarray(jloc,probenum);
	p->Iarray(kloc,probenum);
	p->Iarray(flag,probenum);
    p->Darray(val,probenum*3);
	
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
		sprintf(name,"./REEF3D_NHFLOW_ProbePoint/REEF3D-NHFLOW-Vel-Probe-%i.dat",n+1);
		
		pout[n].open(name);
        
        //cout<<pout[n].is_open()<<" "<<n+1<<endl;

	    pout[n]<<"Vel Probe ID:  "<<n+1<<endl<<endl;
		pout[n]<<"x_coord     y_coord     z_coord"<<endl;
		
		pout[n]<<n+1<<"\t "<<p->Xout(p->P65_x[n],p->P65_y[n])<<"\t "<<p->Yout(p->P65_x[n],p->P65_y[n])<<"\t "<<p->P65_z[n]<<endl;

		pout[n]<<endl<<endl;
		
		pout[n]<<"t \t U \t V \t W "<<endl;
		
		}
    }
	
}

nhflow_probe_vel::~nhflow_probe_vel()
{
    for(n=0;n<probenum;++n)
    pout[n].close();
    
    delete [] pout;
}

void nhflow_probe_vel::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double xp,yp,zp;
    
    ini_location(p,d);

    // the mesh refinement: the finest grid at the probe (one rank writes it)
    amr_points ap;
    if(amr!=nullptr)
    ap.find(amr,p,probenum,p->P65_x,p->P65_y,flag);
	
	for(n=0;n<probenum;++n)
	{
	uval=vval=wval=-1.0e20;
	
		if(amr==nullptr ? flag[n]>0 : ap.mine[n]!=0)
		{
		xp=p->P65_x[n];
		yp=p->P65_y[n];
		zp=p->P65_z[n];

        lexer *q = p;
        fdm_nhf *g = d;
        if(amr!=nullptr && ap.f[n]!=nullptr)
        {
            q = ap.q[n];
            g = static_cast<fdm_nhf*>(ap.f[n]);
        }
    
		uval = q->ccipol4V(g->U, g->WL, g->bed, xp, yp, zp);
		vval = q->ccipol4V(g->V, g->WL, g->bed, xp, yp, zp);
		wval = q->ccipol4V(g->W, g->WL, g->bed, xp, yp, zp);
        
        etaval = q->ccslipol4(g->eta,xp,yp);
        
        if(zp > p->wd + etaval)
         uval = vval = wval = 0.0;
		}
	
	val[n*3+0]=uval;
	val[n*3+1]=vval;
	val[n*3+2]=wval;
	}

    probe_rows(p,pgc,pout,val,probenum,3);	
}

void nhflow_probe_vel::ini_location(lexer *p, fdm_nhf *d)
{
    int check;

    for(n=0;n<probenum;++n)
    {
    check=0;
    
    iloc[n]=p->posc_i(p->P65_x[n]);
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n]=p->posc_j(p->P65_y[n]);
    
	//kloc[n]=p->posf_sig(iloc[n],jloc[n],p->P65_z[n]);
    //check=boundcheck(p,iloc[n],jloc[n],kloc[n],0);
    
    if(iloc[n]>=0 && iloc[n]<p->knox)
    if((jloc[n]>=0 && jloc[n]<p->knoy) || p->j_dir==0)
    check=1;
    
    //cout<<p->mpirank<<" PROBE check: "<<check<<" i: "<<iloc[n]<<" j: "<<jloc[n]<<" k: "<<kloc[n]<<" ZSN: "<<p->ZSN[10+marge]<<endl;
    if(check==1)
    flag[n]=1;
    }
    
}
