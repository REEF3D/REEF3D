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


#include"seastate_coupling.h"
#include"seastate_f.h"
#include"seastate_source.h"
#include"seastate_roller.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_dispersion.h"
#include"seastate_vtp.h"
#include"fdm_seastate.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<vector>

seastate_coupling::seastate_coupling(lexer *p) : pwave(nullptr), pnet(nullptr),
                Fw(p),Gw(p),divM(p),Sxx(p),Sxy(p),Syy(p),Mx(p),My(p),Fdx(p),Fdy(p),Q(p),B(p),
                tlast(0.0), tnext(0.0), tprint(0.0), nstep(0)
{
}

seastate_coupling::~seastate_coupling()
{
    delete pnet;
    delete pwave;
}

void seastate_coupling::ini_wave(lexer *p, ghostcell *pgc, const char *host)
{
    if(p->mpirank==0)
    {
    cout<<endl<<"REEF3D::SEASTATE coupled with "<<host<<" (A 750 1): ";
    if(p->A770==1)
    cout<<"surfbeat (A 770 1), wave step every "<<host<<" step";
    else
    cout<<"coupling interval (wave time step) "<<p->A706<<" s";
    cout<<", wave forcing "<<(p->A751==1 ? "radiation stress" : "vortex force")<<", ramp-up "<<p->A752<<" s, feedback "
        <<(p->A753==0 ? "none" : (p->A753==1 ? "water level" : "water level and currents"))<<endl;
    }

    pwave = new seastate_f(p,pgc);
    pwave->ini_coupled(p,pgc);

    // wave force and mass flux in the SEASTATE VTP files (from the first wave step on)
    pwave->printer()->add_field("Fx",&Fw);
    pwave->printer()->add_field("Fy",&Gw);
    pwave->printer()->add_field("MStx",&Mx);
    pwave->printer()->add_field("MSty",&My);

    // net source terms without wind input for the vortex-force dissipation term
    // (with the roller, the roller dissipation replaces the breaking dissipation)
    if(p->A751==2 && pwave->source()!=nullptr)
    {
    seastate_source_param prm = pwave->source()->param();
    prm.wind = false;

    if(pwave->roller()!=nullptr)
    prm.breaking = false;

    if(prm.any())
    pnet = new seastate_source(*pwave->e->grid,prm);
    }

    tlast = p->simtime;
    tnext = p->simtime + p->A706;
    tprint = p->simtime + (p->P182>0.0 ? p->P182 : 0.0);
    nstep = 0;
}

double seastate_coupling::ramp(lexer *p) const
{
    return (p->A752>0.0) ? std::min(1.0,std::max(p->simtime,0.0)/p->A752) : 1.0;
}

bool seastate_coupling::due(lexer *p, double &dt) const
{
    dt = p->simtime - tlast;

    if(pwave->surfbeat())
    return dt>0.0;

    return !(p->simtime < tnext - 1.0e-9*p->A706);
}

void seastate_coupling::wave_step(lexer *p, ghostcell *pgc, double dt)
{
    pwave->step_coupled(p,pgc,dt);

    tlast = p->simtime;
    tnext = tlast + p->A706;
    ++nstep;

    if(p->mpirank==0 && (!pwave->surfbeat() || p->count%p->P12==0))
    cout<<"SEASTATE wave step "<<nstep<<" at simtime "<<p->simtime<<" (dt "<<dt<<" s), forcing ramp "<<ramp(p)<<endl;
}

void seastate_coupling::print(lexer *p, ghostcell *pgc)
{
    if(p->P10!=1)
    return;

    // coupling interval: every wave step
    if(!pwave->surfbeat())
    {
    pwave->printer()->print2D(p,pwave->e,pgc);
    return;
    }

    // surfbeat: own print interval P 182 [s] or P 181 [iterations]
    if(p->P182>0.0)
    {
        if(p->simtime>=tprint-1.0e-9)
        {
        pwave->printer()->print2D(p,pwave->e,pgc);
        tprint += p->P182;
        }
    }
    else if(p->P181>0 && p->count%p->P181==0)
    pwave->printer()->print2D(p,pwave->e,pgc);
}

void seastate_coupling::environment(lexer *p, slice &eta, slice &WL, slice &U, slice &V)
{
    fdm_seastate *e = pwave->e;

    IMALOOP
    JMALOOP
    {
    const bool inside = i+p->origin_i>=0 && i+p->origin_i<p->gknox && j+p->origin_j>=0 && j+p->origin_j<p->gknoy;

        if(!inside || e->wet0(i,j)==0)
        continue;

        // water level: total water depth of the host
        if(p->A753>=1)
        {
        e->eta(i,j)   = eta(i,j);
        e->depth(i,j) = std::max(WL(i,j),0.0);
        }

        // current: depth-averaged Eulerian velocity
        if(p->A753==2)
        {
        double u = U(i,j), v = V(i,j);

            // radiation stress: the host velocity is the mass transport velocity
            if(p->A751==1 && WL(i,j)>p->A705)
            {
            u -= Mx(i,j)/WL(i,j);
            v -= My(i,j)/WL(i,j);
            }

        e->U(i,j) = u;
        e->V(i,j) = v;
        }
    }
}

// derivatives between wave-active cells: central, one-sided next to dry cells and at the domain edge
double seastate_coupling::ddx(lexer *p, slice &f, int ii, int jj)
{
    fdm_seastate *e = pwave->e;
    i = ii;
    j = jj;

    const bool w = e->wet(i-1,j)==1, o = e->wet(i+1,j)==1;

    if(w && o)
    return (f(i+1,j)-f(i-1,j))/(p->XP[IP1]-p->XP[IM1]);

    if(o)
    return (f(i+1,j)-f(i,j))/(p->XP[IP1]-p->XP[IP]);

    if(w)
    return (f(i,j)-f(i-1,j))/(p->XP[IP]-p->XP[IM1]);

    return 0.0;
}

double seastate_coupling::ddy(lexer *p, slice &f, int ii, int jj)
{
    fdm_seastate *e = pwave->e;
    i = ii;
    j = jj;

    if(p->j_dir==0)
    return 0.0;

    const bool s = e->wet(i,j-1)==1, n = e->wet(i,j+1)==1;

    if(s && n)
    return (f(i,j+1)-f(i,j-1))/(p->YP[JP1]-p->YP[JM1]);

    if(n)
    return (f(i,j+1)-f(i,j))/(p->YP[JP1]-p->YP[JP]);

    if(s)
    return (f(i,j)-f(i,j-1))/(p->YP[JP]-p->YP[JM1]);

    return 0.0;
}

void seastate_coupling::forces(lexer *p, ghostcell *pgc, slice &U, slice &V)
{
    fdm_seastate *e = pwave->e;
    const seastate_grid &g = *e->grid;
    const seastate_roller *proll = pwave->roller();
    const double grav = seastate_gravity;
    const bool vf = (p->A751==2);

    std::vector<double> P(g.nbin), D(g.nbin);

    // spectral integrals per wave-active cell (incl. the ghost cells: the spectra are exchanged)
    IMALOOP
    JMALOOP
    {
    Sxx(i,j)=Sxy(i,j)=Syy(i,j)=Mx(i,j)=My(i,j)=Fdx(i,j)=Fdy(i,j)=Q(i,j)=B(i,j)=0.0;

        if(e->wet(i,j)==0 || e->N->spec(i,j)==nullptr)
        continue;

    const float *N  = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cc = e->cg->spec(i,j);
    const double d  = e->depth(i,j);

        if(vf && pnet!=nullptr)
        pnet->compute(N,d,kc,cc,P.data(),D.data());

        for(int l=0; l<g.nsig; ++l)
        {
        const double sig = g.sig[l];
        const double k = double(kc[l]);

            if(!(k>0.0))
            continue;

        const double n = double(cc[l])*k/sig;
        const double A = seastate_refraction(sig,k,d);          // sig/sinh(2kd)
        const double w = g.dsig[l]*g.dtheta;

            for(int m=0; m<g.ndir; ++m)
            {
            const int bb = g.bin(l,m);
            const double cs = g.costh[m], sn = g.sinth[m];
            const double E = sig*double(N[bb])*w;               // [m^2]

            Mx(i,j) += grav*E*k/sig*cs;
            My(i,j) += grav*E*k/sig*sn;

                if(!vf)
                {
                Sxx(i,j) += grav*E*(n*(cs*cs+1.0)-0.5);
                Sxy(i,j) += grav*E*n*cs*sn;
                Syy(i,j) += grav*E*(n*(sn*sn+1.0)-0.5);
                }
                else
                {
                Q(i,j) += grav*E*(n-0.5);
                B(i,j) += grav*E/sig*k*A;

                    if(pnet!=nullptr)
                    {
                    const double S = (P[bb] - D[bb]*double(N[bb]))*w;   // net source dN/dt dsig dtheta
                    Fdx(i,j) -= grav*S*k*cs;
                    Fdy(i,j) -= grav*S*k*sn;
                    }
                }
            }
        }

        // roller (surfbeat, one frequency): stress, mass flux 2 g R/c, dissipation force
        if(proll!=nullptr && proll->spec(i,j)!=nullptr)
        {
        const float *r = proll->spec(i,j);
        const double sig = g.sig[0];
        const double k = double(kc[0]);
        const double c = k>0.0 ? sig/k : 0.0;

            if(c>0.0)
            for(int m=0; m<g.ndir; ++m)
            {
            const double cs = g.costh[m], sn = g.sinth[m];
            const double Rm = double(r[m])*g.dtheta;             // [m^2]

            Mx(i,j) += 2.0*grav*Rm/c*cs;
            My(i,j) += 2.0*grav*Rm/c*sn;

                if(!vf)
                {
                Sxx(i,j) += grav*Rm*cs*cs;
                Sxy(i,j) += grav*Rm*cs*sn;
                Syy(i,j) += grav*Rm*sn*sn;
                }
                else
                {
                const double Dr = 2.0*grav*proll->beta/c*Rm;      // roller dissipation [m^2/s]
                Fdx(i,j) += grav*Dr/sig*k*cs;
                Fdy(i,j) += grav*Dr/sig*k*sn;
                }
            }
        }
    }

    // forces in the host cells
    SLICELOOP4
    {
    Fw(i,j)=Gw(i,j)=divM(i,j)=0.0;

        if(e->wet(i,j)==0)
        continue;

        if(!vf)
        {
        Fw(i,j) = -(ddx(p,Sxx,i,j) + ddy(p,Sxy,i,j));
        Gw(i,j) = -(ddx(p,Sxy,i,j) + ddy(p,Syy,i,j));
        }
        else
        {
        const int ic=i, jc=j;

        // div(M) in flux form: face values between wave-active cells, no Stokes transport
        // through the domain edges, walls or the shoreline, so the mass source integrates to zero
        const double fe = (e->wet(i+1,j)==1) ? 0.5*(Mx(i,j)+Mx(i+1,j)) : 0.0;
        const double fw = (e->wet(i-1,j)==1) ? 0.5*(Mx(i,j)+Mx(i-1,j)) : 0.0;
        const double fn = (p->j_dir==1 && e->wet(i,j+1)==1) ? 0.5*(My(i,j)+My(i,j+1)) : 0.0;
        const double fs = (p->j_dir==1 && e->wet(i,j-1)==1) ? 0.5*(My(i,j)+My(i,j-1)) : 0.0;
        const double dM  = (fe-fw)/p->DXN[IP] + (fn-fs)/p->DYN[JP];

        const double chi = ddx(p,V,ic,jc) - ddy(p,U,ic,jc);
        const double dQx = ddx(p,Q,ic,jc), dQy = ddy(p,Q,ic,jc);
        i = ic;
        j = jc;

        divM(i,j) = dM;
        Fw(i,j) = Fdx(i,j) + B(i,j)*e->ddx(i,j) - dQx + chi*My(i,j) - U(i,j)*dM;
        Gw(i,j) = Fdy(i,j) + B(i,j)*e->ddy(i,j) - dQy - chi*Mx(i,j) - V(i,j)*dM;
        }
    }
}

bool seastate_coupling::longwave_bc(lexer *p, int jj, double h, double eta1, double v1, double &etag, double &ug, double &vg)
{
    if(pwave==nullptr || !pwave->surfbeat() || p->origin_i!=0 || !(h>p->A705))
    return false;

    double z=0.0, qx=0.0, qy=0.0;
    pwave->longwave(p,jj,p->simtime,z,qx,qy);

    etag = eta1;
    ug = qx/h - std::sqrt(std::fabs(p->W22)/h)*(eta1 - z);
    vg = v1;

    return true;
}
