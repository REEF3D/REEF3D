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
#include"print_wsf.h"
#include"lexer.h"
#include"amr_grid_lookup.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>
#include"runlog.h"

print_wsf::print_wsf(lexer *p, fdm* a, ghostcell *pgc, int num) : fileFlushMaxCount(100)
{
	gauge_num = p->P51;
	x = p->P51_x;
	y = p->P51_y;
	
    if(p->P51>0 && num==0)
	{
	gauge_num = p->P51;
	x = p->P51_x;
	y = p->P51_y;
	}
	
	if(p->P351>0 && num==1)
	{
	gauge_num = p->P351;
	x = p->P351_x;
	y = p->P351_y;
	}
	
	if(p->P352>0 && num==2)
	{
	gauge_num = p->P352;
	x = p->P352_x;
	y = p->P352_y;
	}
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_WSF",0777);
	
    if(p->mpirank==0 && p->P51>0 && num==0)
    {
    // open file
	wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG.dat");
	if(p->plog)
	p->plog->table_file(p,"wsf","gauges","./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }
	
	//-------------------
	
	if(p->mpirank==0 && p->P351>0 && num==1)
    {
    // open file
	wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-1.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }
	
	//-------------------
	
	if(p->mpirank==0 && p->P352>0 && num==2)
    {
    // open file
	wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-2.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }
	
	p->Iarray(iloc,gauge_num);
	p->Iarray(jloc,gauge_num);
	p->Iarray(flag,gauge_num);
	p->Darray(wsf,gauge_num);

    ini_location(p,a,pgc);
}

print_wsf::~print_wsf()
{
    wsfout.close();
}

void print_wsf::height_gauge(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    double zval=0.0;

    for(n=0;n<gauge_num;++n)
    wsf[n]=-1.0e20;

	if(p->A10==6 && p->F80!=4)
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
        KLOOP
        PCHECK
        {
            if(f(i,j,k)>=0.0 && f(i,j,k+1)<0.0)
            wsf[n]=MAX(wsf[n],-(f(i,j,k)*p->DZP[KP])/(f(i,j,k+1)-f(i,j,k)) + p->pos_z());
        }

        // CFD mesh refinement: when the finest grid at the free surface of level 0 is a patch, the
        // free surface of the patch column that holds the gauge (its topmost zero, as above)
        if(amr!=nullptr && wsf[n]>-1.0e19)
        amr_gauge(p,n);
    }
    
    if(p->A10==6 && p->F80==4)
    {
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
        KLOOP
        {
            if(f(i,j,k)>p->F94 && f(i,j,k+1)<p->F93)
                wsf[n]=MAX(wsf[n],p->pos_z()+0.5*p->DZN[KP]);
            else if(f(i,j,k)<=p->F94 && f(i,j,k)>=p->F93)
            {
                /*if(a->Alpha(i,j,k)<1E05)
                {
                    if(a->nZ(i,j,k)>1E-06 || a->nZ(i,j,k)<-1E-06)
                    {
                        if(fabs(a->Alpha(i,j,k)/a->nZ(i,j,k))<0.5*p->DZN[KP])
                            wsf[n]=MAX(wsf[n],p->pos_z()+a->Alpha(i,j,k)/a->nZ(i,j,k));
                        else
                        {
                            if(a->nZ(i,j,k)>0.0)
                                wsf[n]=MAX(wsf[n],p->pos_z()-0.5*p->DZN[KP]);
                            else
                                wsf[n]=MAX(wsf[n],p->pos_z()+0.5*p->DZN[KP]);
                        }
                    }
                    else
                        wsf[n]=MAX(wsf[n],p->pos_z()+0.5*p->DZN[KP]);
                }
                else
                {
                    cout<<"surface cell in wave gauge does not have a plane"<<endl;*/
                    wsf[n]=MAX(wsf[n],(p->pos_z()-0.5*p->DZN[KP])+f(i,j,k)*p->DZN[KP]);
               // }
            }
        }
    }
    }
    
    if(p->A10==5)
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
			wsf[n] = a->eta(i,j);

    }
	
    pgc->globalmax(wsf,gauge_num);

    // write to file
    gauge_row(p,wsfout,wsf,gauge_num,fileFlushMaxCount);
}

// the free surface of gauge n on the finest grid at (x, y, the free surface of level 0): the
// topmost zero of the level set in the column of that grid that holds the gauge; none in that
// column (the free surface left the patch vertically): level 0
void print_wsf::amr_gauge(lexer *p, int n)
{
    lexer *q = nullptr;
    fdm *fa = nullptr;
    if(!amr->finest_grid(x[n],y[n],wsf[n],q,fa))
    return;

    auto cell = [&](const double *N, int nc, double v) -> int
    {
        for(int c=0; c<nc; ++c)
        if(v>=N[c+marge] && v<N[c+1+marge])
        return c;
        return -1;
    };
    const int ii = cell(q->XN,q->knox,x[n]);
    const int jj = (q->j_dir==1) ? cell(q->YN,q->knoy,y[n]) : 0;
    if(ii<0 || jj<0)
    return;

    double zs = -1.0e20;
    for(int kk=0; kk<q->knoz; ++kk)
    {
        const double f0 = fa->phi(ii,jj,kk), f1 = fa->phi(ii,jj,kk+1);
        if(f0>=0.0 && f1<0.0)
        zs = MAX(zs,-(f0*q->DZP[kk+marge])/(f1-f0) + q->ZP[kk+marge]);
    }
    if(zs>-1.0e19)
    wsf[n] = zs;
}

void print_wsf::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void print_wsf::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check;

    for(n=0;n<gauge_num;++n)
    {
    iloc[n] = p->posc_i(x[n]); 
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n] = p->posc_j(y[n]); 

    check=ij_boundcheck(p,iloc[n],jloc[n],0);

    if(check==1)
    flag[n]=1;
    
    //cout<<p->mpirank<<" n: "<<n<<" flag: "<<flag[n]<<" x: "<<x[n]<<" y: "<<y[n]<<" iloc: "<<iloc[n]<<" jloc: "<<jloc[n]<<endl;
    }
}


