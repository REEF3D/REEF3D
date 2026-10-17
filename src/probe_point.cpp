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

#include"probe_point.h"
#include<vector>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"amr_grid_lookup.h"
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
    std::vector<double> val(7*probenum,-1.0e20);   // u v w p k eps eddyv per probe, one reduction for all

    // the probes of this rank: on its level-0 grid; with the placement of the CFD patches (G 40)
    // the rank of the finest grid there
    std::vector<char> mine;
    amr_cfd_mine(amr,probenum,p->P61_x,p->P61_y,p->P61_z,flag,mine);
	
	for(n=0;n<probenum;++n)
	if(mine[n])
	{
		xp=p->P61_x[n];
		yp=p->P61_y[n];
		zp=p->P61_z[n];
		
        // the grid of the values: level 0, or the finest patch of a mesh refinement there
        lexer *q = p;
        fdm *fa = a;
        turbulence *t = pturb;
        if(amr!=nullptr)
        amr->finest_grid(xp,yp,zp,q,fa,t);

		val[7*n+0] = q->ccipol1(fa->u, xp, yp, zp);
		val[7*n+1] = q->ccipol2(fa->v, xp, yp, zp);
		val[7*n+2] = q->ccipol3(fa->w, xp, yp, zp);
		val[7*n+3] = q->ccipol4a(fa->press, xp, yp, zp) - p->pressgage;
		val[7*n+4] = t->ccipol_kinval(q, pgc, xp, yp, zp);
		val[7*n+5] = t->ccipol_epsval(q, pgc, xp, yp, zp);
		val[7*n+6] = q->ccipol4a(fa->eddyv, xp, yp, zp);
	}
	
	pgc->globalmax(val.data(),7*probenum);

	if(p->mpirank==0)
	for(n=0;n<probenum;++n)
	pout[n]<<setprecision(9)<<p->simtime<<" \t "<<val[7*n+0]<<" \t "<<val[7*n+1]<<" \t "<<val[7*n+2]<<" \t "<<val[7*n+3]<<" \t "<<val[7*n+4]<<" \t "<<val[7*n+5]<<" \t "<<val[7*n+6]<<endl;
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

    check=boundcheck(p,iloc[n],jloc[n],kloc[n],0);

    if(check==1)
    flag[n]=1;
    }
}


