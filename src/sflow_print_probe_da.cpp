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

#include"sflow_print_probe_da.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_print_probe_da::sflow_print_probe_da(lexer *p, fdm2D *b, ghostcell *pgc) : probenum(p->P63)
{
    p->Iarray(iloc,probenum);
	p->Iarray(jloc,probenum);
	p->Iarray(flag,probenum);
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_SFLOW_ProbePoint",0777);
	
	pout = new ofstream[probenum];
	
    if(p->mpirank==0 && probenum>0)
    {
		cout<<"probepoint_num: "<<probenum<<endl;
		// open file
		for(n=0;n<probenum;++n)
		{
		sprintf(name,"./REEF3D_SFLOW_ProbePoint/REEF3D-SFLOW-Probe-Point-%i.dat",n+1);
		
		pout[n].open(name);

	    pout[n]<<"Depth Averaged Point Probe ID:  "<<n<<endl<<endl;
		pout[n]<<"x_coord     y_coord"<<endl;
		
		pout[n]<<n+1<<"\t "<<p->P63_x[n]<<"\t "<<p->P63_y[n]<<endl;

		pout[n]<<endl<<endl;
		
		pout[n]<<"t \t Um \t Vm \t Wm \t Pm \t eta"<<endl;
		}
    }
	
    ini_location(p,b,pgc);
}

sflow_print_probe_da::~sflow_print_probe_da()
{
	for(n=0;n<probenum;++n)
    pout[n].close();
}

void sflow_print_probe_da::start(lexer *p, fdm2D *b, ghostcell *pgc)
{
	double xp,yp;
	
	for(n=0;n<probenum;++n)
	{
	uval=vval=wval=pval=eval=-1.0e20;
	
		if(flag[n]>0)
		{
		xp=p->P63_x[n];
		yp=p->P63_y[n];
        
		uval = p->ccslipol1(b->P, xp, yp);
		vval = p->ccslipol2(b->Q, xp, yp);
		wval = p->ccslipol4(b->ws,xp,yp); 
		pval = p->ccslipol4(b->press,xp,yp);
        eval = p->ccslipol4(b->eta,xp,yp); 
		}
	
	uval=pgc->globalmax(uval);
	vval=pgc->globalmax(vval);
	wval=pgc->globalmax(wval);
	pval=pgc->globalmax(pval);
    eval=pgc->globalmax(eval);
	
	if(p->mpirank==0)
	pout[n]<<setprecision(9)<<p->simtime<<" \t "<<uval<<" \t "<<vval<<" \t "<<wval<<" \t "<<pval<<" \t "<<eval<<endl;
	}
			
}

void sflow_print_probe_da::write(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

void sflow_print_probe_da::ini_location(lexer *p, fdm2D *b, ghostcell *pgc)
{
    int check;

    for(n=0;n<probenum;++n)
    {
    check=0;
    
    iloc[n]=p->posc_i(p->P63_x[n]);
    
    if(p->j_dir==0)
    jloc[n]=0;
  
    if(p->j_dir==1)
    jloc[n]=p->posc_j(p->P63_y[n]);
    
    if(iloc[n]>=0 && iloc[n]<p->knox)
    if(jloc[n]>=0 && jloc[n]<p->knoy)
    check=1;
    
    if(check==1)
    {
    i = iloc[n];
    j = jloc[n];
    
    if(p->flagslice4[IJ]<0)
    check=0;
    }

    if(check==1)
    flag[n]=1;
    }
}


