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

#include"fnpf_vel_probe_theory.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ioflow.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_vel_probe_theory::fnpf_vel_probe_theory(lexer *p, fdm_fnpf *c) : probenum(p->P66)
{

    p->Iarray(iloc,probenum);
	p->Iarray(jloc,probenum);
	p->Iarray(kloc,probenum);
	p->Iarray(flag,probenum);
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_FNPF_ProbePoint",0777);
	
	pout = new ofstream[probenum];
	
    if(p->mpirank==0 && probenum>0)
    {
		cout<<"probepoint_num: "<<probenum<<endl;
		// open file
		for(n=0;n<probenum;++n)
		{
		sprintf(name,"./REEF3D_FNPF_ProbePoint/REEF3D-FNPF-Vel-Probe-Theory-%i.dat",n+1);
		
		pout[n].open(name);
        
        //cout<<pout[n].is_open()<<" "<<n+1<<endl;

	    pout[n]<<"Vel Probe Theory ID:  "<<n<<endl<<endl;
		pout[n]<<"x_coord     y_coord     z_coord"<<endl;
		
		pout[n]<<n+1<<"\t "<<p->P66_x[n]<<"\t "<<p->P66_y[n]<<"\t "<<p->P66_z[n]<<endl;

		pout[n]<<endl<<endl;
		
		pout[n]<<"t \t U \t V \t W "<<endl;
		
		}
    }
	
}

fnpf_vel_probe_theory::~fnpf_vel_probe_theory()
{
    for(n=0;n<probenum;++n)
    pout[n].close();
    
    delete [] pout;
}

void fnpf_vel_probe_theory::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow)
{
    double xp,yp,zp;
    
    ini_location(p,c);
	
	for(n=0;n<probenum;++n)
	{
	uval=vval=wval=-1.0e20;
	
		if(flag[n]>0)
		{
		xp=p->P66_x[n];
		yp=p->P66_y[n];
		zp=p->P66_z[n];
    
		uval = pflow->wave_xvel(p,pgc, xp, yp, zp);
		vval = pflow->wave_yvel(p,pgc, xp, yp, zp);
		wval = pflow->wave_zvel(p,pgc, xp, yp, zp);
		}
	
	uval=pgc->globalmax(uval);
	vval=pgc->globalmax(vval);
	wval=pgc->globalmax(wval);

	
	if(p->mpirank==0)
	pout[n]<<setprecision(9)<<p->simtime<<" \t "<<uval<<" \t "<<vval<<" \t "<<wval<<endl;
	}	
}

void fnpf_vel_probe_theory::ini_location(lexer *p, fdm_fnpf *c)
{
    int check;

    for(n=0;n<probenum;++n)
    {
    check=0;
    
    iloc[n]=p->posc_i(p->P66_x[n]);
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n]=p->posc_j(p->P66_y[n]);
    
	kloc[n]=p->posf_sig(iloc[n],jloc[n],p->P66_z[n]);

    check=boundcheck(p,iloc[n],jloc[n],kloc[n],0);
    //cout<<p->mpirank<<" PROBE check: "<<check<<" i: "<<iloc[n]<<" j: "<<jloc[n]<<" k: "<<kloc[n]<<" ZSN: "<<p->ZSN[10+marge]<<endl;
    if(check==1)
    flag[n]=1;
    }
}
