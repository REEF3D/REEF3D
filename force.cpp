/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include "force.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "gradient.h"
#include<sys/stat.h>
#include<sys/types.h>

force::force(lexer *p, fdm *a, ghostcell *pgc, int ID) : wave_interface(p,pgc), epsi(0.0*p->dx)
{	
    
	if(p->P81>0)
	{
	is = p->posc_i(p->P81_xs[ID]);
    ie = p->posc_i(p->P81_xe[ID]);

    js = p->posc_j(p->P81_ys[ID]);
    je = p->posc_j(p->P81_ye[ID]);
	
	ks = p->posc_k(p->P81_zs[ID]);
    ke = p->posc_k(p->P81_ze[ID]);
	
	xs = p->P81_xs[ID];
	xe = p->P81_xe[ID];
	
	ys = p->P81_ys[ID];
	ye = p->P81_ye[ID];
	
	zs = p->P81_zs[ID];
	ze = p->P81_ze[ID];
	
	xm = xs + (xe-xs)*0.5;
	ym = ys + (ye-ys)*0.5;
	zm = zs + (ze-zs)*0.5;
	
	fsfmargin = p->P83;
	}
	
	
	if(p->P85>0)
	{
	is = p->posc_i(p->P85_xs[ID]);
    ie = p->posc_i(p->P85_xe[ID]);

    js = p->posc_j(p->P85_ys[ID]);
    je = p->posc_j(p->P85_ye[ID]);
	
	ks = p->posc_k(p->P85_zs[ID]);
    ke = p->posc_k(p->P85_ze[ID]);
	
	xs = p->P85_xs[ID];
	xe = p->P85_xe[ID];
	
	ys = p->P85_ys[ID];
	ye = p->P85_ye[ID];
	
	zs = p->P85_zs[ID];
	ze = p->P85_ze[ID];
	
	xm = xe;// + (xe-xs)*0.5;
	ym = ye;// + (ye-ys)*0.5;
	zm = ze;// + (ze-zs)*0.5;
	
	fsfmargin = p->P87;
	}
	
	
	x_center = xs + 0.5*(xe-xs);
	y_center = ys + 0.5*(ye-ys);
	z_center = zs + 0.5*(ze-zs);
	
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_Force",0777);
	
    if(p->mpirank==0)
    {
    // open force surf file
	if(p->P14==0)
	sprintf(name,"REEF3D_Force-%d.dat",ID+1);
    
	if(p->P14==1)
	sprintf(name,"./REEF3D_Force/REEF3D_Force-%d.dat",ID+1);
	
	fout.open(name);

    fout<<"x_start xend     y_start y_end     z_start z_end"<<endl;

    fout<<xs<<" "<<xe<<" . "<<ys<<" "<<ye<<" . "<<zs<<" "<<ze<<endl;
    fout<<endl<<endl;
    
    if(p->P93==0)
    fout<<"it \t time \t Fx \t Fz ";
    
    if(p->P93==1)
    fout<<"it \t time \t Fx \t Fz \t Fmor \t FmorRCT \t tnorm \t Fz_Norm \t Fx_norm \t Fmor_norm \t Fx_Press \t Fx_Shear ";
    
    fout<<endl;
	}
	
	for(int qn=0;qn<p->P81;++qn)
	ini(p,a,pgc);
	
	for(int qn=0;qn<p->P85;++qn)
	ini(p,a,pgc);
}

force::~force()
{
}

void force::ini(lexer *p, fdm *a, ghostcell *pgc)
{	
    // ini forcesurf
	if(p->P81>0)
	{
		count=0;
		for(n=0;n<p->facetnum;++n)
		{
		i=p->facet[n][0];
		j=p->facet[n][1];
		k=p->facet[n][2];

		if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke)
		++count;
		}

		fnum=count;
	}
	
	if(p->P85>0)
	{
		count=0;
		GC4LOOP
		{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];

		if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke && p->gcb4[n][4]==21)
		++count;
		}

		fnum=count;
	}

    p->Darray(ccpt,fnum,6,3);
    p->Darray(farea,fnum);
    p->Darray(fn,fnum,3);
	p->Darray(fx,fnum);
	p->Darray(fy,fnum);
	p->Darray(fz,fnum);
	p->Iarray(fcheck,fnum);
    p->Iarray(fid,fnum,3);
    p->Iarray(surfnum,fnum);
	
	if(p->P81>0)
	{
		count=0;
		for(n=0;n<p->facetnum;++n)
		{
		i=p->facet[n][0];
		j=p->facet[n][1];
		k=p->facet[n][2];

			if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke)
			{
			if(a->topo(i,j,k)>=-1.0e-7)
			fcheck[count]=1;

			++count;
			}
		}
	}
	
	if(p->P85>0)
	{
		count=0;
		GC4LOOP
		{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];

			if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke && p->gcb4[n][4]==21)
			{
			if(a->topo(i,j,k)>=-1.0e-7)
			fcheck[count]=1;
			
			++count;
			}
		}
	}

//-------------------------
	if(p->P81>0)
	{
    cellnodes(p,a,pgc);
    surfarea(p,a,pgc);
	}
//-------------------------
	if(p->P85>0)
	{
    cellnodes_gcb(p,a,pgc);
    surfarea_gcb(p,a,pgc);
	}
//-------------------------

    un1=un2=vn1=vn2=dt1=dt2=0.0;
    udt1=udt2=0.0;
	Fvert0=0.0;
	
}

void force::start(lexer *p, fdm *a, ghostcell *pgc)
{
	
    force_surface(p,a,pgc);
    crossarea(p,a,pgc);
    
    if(p->P93==1)
    {
    velocity(p,a,pgc);
    coefficients(p,a,pgc);
    morison(p,a,pgc);
    }
    
    if(p->mpirank==0)
    {
    cout<<"Fx: "<<FDs<<" Fz: "<<Fvert<<endl;
    print(p,a,pgc);
    }
}
