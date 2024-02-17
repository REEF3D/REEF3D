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

#include"probe_point.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<sys/stat.h>
#include<sys/types.h>

probe_point::probe_point(lexer *p, fdm* a, ghostcell *pgc) : probenum(p->P61)
{
    p->Iarray(iloc,probenum);
	p->Iarray(jloc,probenum);
	p->Iarray(kloc,probenum);
	p->Iarray(flag,probenum);
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_ProbePoint",0777);
	
	pout = new ofstream[probenum];
	
    if(p->mpirank==0 && probenum>0)
    {
		cout<<"probepoint_num: "<<probenum<<endl;
		// open file
		for(n=0;n<probenum;++n)
		{
		sprintf(name,"./REEF3D_CFD_ProbePoint/REEF3D-CFD-Probe-Point-%i.dat",n+1);
		
		pout[n].open(name);
        
        //cout<<pout[n].is_open()<<" "<<n+1<<endl;

	    pout[n]<<"Point Probe ID:  "<<n<<endl<<endl;
		pout[n]<<"x_coord     y_coord     z_coord"<<endl;
		
		pout[n]<<n+1<<"\t "<<p->P61_x[n]<<"\t "<<p->P61_y[n]<<"\t "<<p->P61_z[n]<<endl;

		pout[n]<<endl<<endl;
		
		pout[n]<<"t \t U \t V \t W \t P \t Kin \t Eps/Omega \t Eddyv"<<endl;
		
		}
    }
	
    ini_location(p,a,pgc);
}

probe_point::~probe_point()
{
	for(n=0;n<probenum;++n)
    pout[n].close();
}

void probe_point::start(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{
	double xp,yp,zp;
	
	for(n=0;n<probenum;++n)
	{
	uval=vval=wval=pval=kval=eval=edval=-1.0e20;
	
		if(flag[n]>0)
		{
		xp=p->P61_x[n];
		yp=p->P61_y[n];
		zp=p->P61_z[n];
		
		uval = p->ccipol1(a->u, xp, yp, zp);
		vval = p->ccipol2(a->v, xp, yp, zp);
		wval = p->ccipol3(a->w, xp, yp, zp);
		pval = p->ccipol4_a(a->press, xp, yp, zp) - p->pressgage;
		kval = pturb->ccipol_kinval(p, pgc, xp, yp, zp);
		eval = pturb->ccipol_epsval(p, pgc, xp, yp, zp);
		edval = p->ccipol4_a(a->eddyv, xp, yp, zp);
		}
	
	uval=pgc->globalmax(uval);
	vval=pgc->globalmax(vval);
	wval=pgc->globalmax(wval);
	pval=pgc->globalmax(pval);
	kval=pgc->globalmax(kval);
	eval=pgc->globalmax(eval);
	edval=pgc->globalmax(edval);

	
	if(p->mpirank==0)
	pout[n]<<setprecision(9)<<p->simtime<<" \t "<<uval<<" \t "<<vval<<" \t "<<wval<<" \t "<<pval<<" \t "<<kval<<" \t "<<eval<<" \t "<<edval<<endl;
	}			
}

void probe_point::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void probe_point::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check;

    for(n=0;n<probenum;++n)
    {
    check=0;
    
    iloc[n]=p->posc_i(p->P61_x[n]);
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n]=p->posc_j(p->P61_y[n]);
    
	kloc[n]=p->posc_k(p->P61_z[n]);

    check=boundcheck(p,a,iloc[n],jloc[n],kloc[n],0);

    if(check==1)
    flag[n]=1;
    }
}


