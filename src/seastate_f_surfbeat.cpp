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
Architect: Hans Bihs
--------------------------------------------------------------------*/


#include"seastate_f.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"seastate_surfbeat.h"
#include"seastate_roller.h"
#include"seastate_source.h"
#include"seastate_implicit.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

/*--------------------------------------------------------------------
Surfbeat (A 770 1): the wave-group action balance on one representative
frequency. See seastate_f.h and seastate_surfbeat.h.
--------------------------------------------------------------------*/

void seastate_f::surfbeat_input(lexer *p, ghostcell *pgc)
{
    // boundary spectrum on the multi-frequency grid A 701-703
    gin = new seastate_grid(p->A701,p->A702_fmin,p->A702_fmax,p->A703);

    if(!gin->valid())
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  "<<gin->message()<<endl<<endl;

        pgc->final(true);
    }

    boundary_spectrum(p,pgc,*gin,Nin);

    seastate_param sp;
    sp.compute(*gin,Nin.data());

    if(!(sp.m0>0.0))
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  A 770 1: the boundary spectrum has no energy"<<endl<<endl;

        pgc->final(true);
    }

    trep = p->A771>0.0 ? p->A771 : sp.Tm10;

    if(p->mpirank==0)
    cout<<"SEASTATE surfbeat: representative period "<<setprecision(5)<<trep<<" s ("<<(p->A771>0.0 ? "A 771" : "Tm-1,0 of the boundary spectrum")
        <<"), boundary Hm0 "<<sp.Hs<<" m"<<endl;
}

void seastate_f::surfbeat_ini(lexer *p, ghostcell *pgc)
{
    const seastate_grid &g = *e->grid;
    const int ndir = g.ndir;

    sb = new seastate_surfbeat(*gin,Nin.data(),p->A772,p->A773);

    // boundary rows of this rank on the x- side: y and water depth of the boundary cells
    std::vector<double> y, d;

    if(p->origin_i==0)
    for(j=0; j<p->knoy; ++j)
    {
    y.push_back(p->YP[JP]);
    d.push_back(std::max(e->depth(0,j),std::max(p->A705,1.0e-3)));
    }

    sb->series(y,d,p->A774==1);

    Nbx.assign(size_t(std::max(p->knoy,1))*ndir,0.0f);
    Nbx0 = Nbx;

    // mean spectrum for the other sides with A 712 1
    Nb.assign(ndir,0.0f);
    for(int m=0; m<ndir; ++m)
    Nb[m] = float(sb->m0*sb->Dbar[m]/(g.sig[0]*g.dsig[0]));

    if(p->mpirank==0)
    {
    cout<<"SEASTATE surfbeat boundary (x-): "<<sb->K<<" components, df "<<setprecision(4)<<sb->df<<" Hz, record "<<sb->trec<<" s, time step of the series "<<sb->dtbc
        <<" s, seed "<<p->A773<<", long waves "<<(p->A774==1 ? "bound (Herbers 1994)" : "none")<<endl;
    cout<<"SEASTATE surfbeat: Roelvink breaking "<<(p->A740==2 ? "on" : "off")<<" (alpha "<<p->A741_alpha<<", gamma "<<p->A741_gamma<<", n "<<p->A746<<"), H/h <= "<<p->A747
        <<", roller "<<(p->A748==1 && p->A740==2 ? "on" : "off")<<" (beta "<<p->A749<<")"<<endl<<endl;
    cout<<setprecision(6);
    }

    surfbeat_boundary(p,p->simtime);

    // components and the boundary series of the first row (global j = 0)
    if(p->origin_i==0 && p->origin_j==0)
    {
    mkdir("./REEF3D_SEASTATE_Log",0777);

    ofstream cf("./REEF3D_SEASTATE_Log/REEF3D_SEASTATE_surfbeat_components.dat");
    cf<<"REEF3D::SEASTATE surfbeat boundary components: eta = sum a cos(2 pi f t - k (cos(theta) x + sin(theta) y) + phase), x = 0 at the x- side"<<endl;
    cf<<"T_rep "<<setprecision(10)<<trep<<" s, T_rec "<<sb->trec<<" s, df "<<sb->df<<" Hz, m0 "<<sb->m0<<" m^2, seed "<<p->A773<<endl;
    cf<<"#f [Hz] \t #a [m] \t #theta [rad] \t #phase [rad]"<<endl;
    for(int m=0; m<sb->K; ++m)
    cf<<setprecision(12)<<sb->fk[m]<<" \t "<<sb->ak[m]<<" \t "<<sb->thk[m]<<" \t "<<sb->phk[m]<<endl;
    cf.close();

    ofstream bf("./REEF3D_SEASTATE_Log/REEF3D_SEASTATE_surfbeat_boundary.dat");
    bf<<"REEF3D::SEASTATE surfbeat boundary series of the first row (y = "<<setprecision(10)<<p->YP[0+marge]<<" m, depth "<<e->depth(0,0)<<" m), without ramp"<<endl;
    bf<<"#t [s] \t #E [m^2] \t #zeta_b [m] \t #q_bx [m^2/s] \t #q_by [m^2/s]"<<endl;
    const int nt = int(std::lround(sb->trec/sb->dtbc));
        for(int n=0; n<nt; ++n)
        {
        double E, z, qx, qy;
        sb->at(0,n*sb->dtbc,E,z,qx,qy);
        bf<<setprecision(10)<<n*sb->dtbc<<" \t "<<E<<" \t "<<z<<" \t "<<qx<<" \t "<<qy<<endl;
        }
    bf.close();
    }
}

void seastate_f::surfbeat_boundary(lexer *p, double t)
{
    if(sb==nullptr || p->origin_i!=0)
    return;

    const seastate_grid &g = *e->grid;
    const int ndir = g.ndir;
    const double r = (p->A752>0.0) ? std::min(1.0,std::max(t,0.0)/p->A752) : 1.0;

    // rows of the previous time level (Crank-Nicolson)
    Nbx0 = Nbx;

    for(int jj=0; jj<p->knoy; ++jj)
    {
    double E, z, qx, qy;
    sb->at(jj,t,E,z,qx,qy);

        for(int m=0; m<ndir; ++m)
        Nbx[size_t(jj)*ndir+m] = float(r*E*sb->Dbar[m]/(g.sig[0]*g.dsig[0]));
    }
}

bool seastate_f::longwave(lexer *p, int jj, double t, double &zeta, double &qx, double &qy) const
{
    zeta = qx = qy = 0.0;

    if(sb==nullptr || p->origin_i!=0 || jj<0 || jj>=p->knoy || p->A774!=1)
    return false;

    double E;
    sb->at(jj,t,E,zeta,qx,qy);

    const double r = (p->A752>0.0) ? std::min(1.0,std::max(t,0.0)/p->A752) : 1.0;
    zeta *= r;
    qx *= r;
    qy *= r;

    return true;
}

void seastate_f::surfbeat_cap(lexer *p)
{
    // H = sqrt(8 E) <= A 747 h (the waves cannot be higher than the water is deep)
    const seastate_grid &g = *e->grid;
    const int ndir = g.ndir;
    const double gmax = p->A747;

    if(!(gmax>0.0))
    return;

    IMALOOP
    JMALOOP
    if(e->wet(i,j)==1)
    {
    float *N = e->N->spec(i,j);
    double E = 0.0;

        for(int m=0; m<ndir; ++m)
        E += g.sig[0]*double(N[m])*g.dsig[0]*g.dtheta;

    const double Hm = gmax*e->depth(i,j);

        if(8.0*E>Hm*Hm && E>0.0)
        {
        const double f = Hm*Hm/(8.0*E);

        for(int m=0; m<ndir; ++m)
        N[m] = float(f*double(N[m]));
        }
    }
}

void seastate_f::surfbeat_step(lexer *p, ghostcell *pgc)
{
    surfbeat_cap(p);

    if(proll!=nullptr)
    proll->step(p,pgc,e,psrc,dtw,iter_max);
}
