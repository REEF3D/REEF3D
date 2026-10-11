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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#include"cavitation.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"hires.h"
#include<cmath>
#include<cstring>
#include<iomanip>
#include<iostream>
#include<sys/stat.h>
#include<sys/types.h>

cavitation::cavitation(lexer *p, fdm *a, ghostcell *pgc) : alpha(p), Cs(p), Sv(p), Dv(p)
{
    if(p->F30>0 || p->F80>0)
    {
        if(p->mpirank==0)
        std::cout<<"W 60: cavitation is single phase, use F 30 0 and F 80 0"<<std::endl;
        pgc->final(true);
    }

    rhol = p->W1;
    nul = p->W2;
    rhov = p->W62;
    nuv = p->W63;
    pv = p->W61;
    n0 = p->W64;
    Cv = p->W66_v;
    Cc = p->W66_c;
    const double vnuc = n0*M_PI*p->W65*p->W65*p->W65/6.0;
    alnuc = vnuc/(1.0 + vnuc);
    eps = 1.0;              // [Pa], regularises 1/sqrt|p - p_v| at p = p_v
    gcval_alpha = 1;
    Smax = 0.0;

    pconvec = new hires(p,11);  // van Leer

    LOOP
    {
    alpha(i,j,k) = 0.0;
    Cs(i,j,k) = 0.0;
    Sv(i,j,k) = 0.0;
    Dv(i,j,k) = 0.0;

    // initial vapour: box W 67, sphere W 68 (cell centres)
    const double xc = p->pos_x(), yc = p->pos_y(), zc = p->pos_z();
    if(p->W67==1 && xc>=p->W67_xs && xc<p->W67_xe && yc>=p->W67_ys && yc<p->W67_ye && zc>=p->W67_zs && zc<p->W67_ze)
    alpha(i,j,k) = 1.0;
    if(p->W68==1 && (xc-p->W68_x)*(xc-p->W68_x) + (yc-p->W68_y)*(yc-p->W68_y) + (zc-p->W68_z)*(zc-p->W68_z) < p->W68_r*p->W68_r)
    alpha(i,j,k) = 1.0;
    }
    pgc->start4(p,alpha,gcval_alpha);
    properties(p,a,pgc);

    active = this;

    if(p->mpirank==0)
    {
        std::cout<<"cavitation: Schnerr-Sauer, p_v "<<pv<<" rho_v "<<rhov<<" nu_v "<<nuv<<" n0 "<<n0
                 <<" alpha_nuc "<<alnuc<<" C_v "<<Cv<<" C_c "<<Cc<<std::endl;
        if(p->B30>0)
        std::cout<<"cavitation: B 30 shifts the pressure level after each projection; fix the level by an outlet pressure patch (B 412)"<<std::endl;

        mkdir("./REEF3D_CFD_Cavitation",0777);
        const bool restart = p->I40>0;
        out.open("./REEF3D_CFD_Cavitation/REEF3D-CFD-Cavitation.dat",restart ? std::ios::app : std::ios::out);
        if(!restart)
        out<<"time \t vapour_volume \t equivalent_radius \t vapour_production \t p_min"<<std::endl;
    }
}

cavitation::~cavitation()
{
    if(active==this)
    active = nullptr;
}

double cavitation::rho_mix(lexer *p, double al) const
{
    return (1.0-al)*rhol + al*rhov;
}

void cavitation::coefficients(lexer *p, fdm *a, ghostcell *pgc)
{
    const double rate = std::sqrt(2.0/(3.0*rhol));

    LOOP
    {
        const double al = alpha(i,j,k);
        const double dp = a->press(i,j,k) - pv;
        double c = 0.0;

        if(dp<0.0)
        {
            // evaporation, from the nuclei in pure liquid
            const double ae = std::min(al + alnuc,1.0);
            const double RB = std::cbrt(3.0*ae/(4.0*M_PI*n0*std::max(1.0-al,1.0e-6)));
            c = Cv*ae*(1.0-al)*3.0/RB*rate/std::sqrt(std::fabs(dp)+eps);
        }
        else if(al>1.0e-10)
        {
            // condensation; max(1-alpha, alpha_nuc) lets a pure vapour cell condense
            const double RB = std::cbrt(3.0*al/(4.0*M_PI*n0*std::max(1.0-al,1.0e-6)));
            c = Cc*al*std::max(1.0-al,alnuc)*3.0/RB*rate/std::sqrt(std::fabs(dp)+eps);
        }

        // no phase change in a direct-forcing body
        if(p->X10==1)
        c *= std::max(0.0,1.0-a->fbh4(i,j,k));

        Cs(i,j,k) = c;
    }
}

void cavitation::projection(lexer *p, fdm *a, double ark)
{
    // div u = Cs (rho_l - rho_v)/rho (p_v - p - pcorr): implicit in the pressure correction
    n=0;
    LOOP
    {
        const double c = Cs(i,j,k)*(rhol-rhov)/rho_mix(p,alpha(i,j,k));
        if(c>0.0)
        {
        a->M.p[n] += c/(ark*p->dt);
        a->rhsvec.V[n] += c*(pv - a->press(i,j,k))/(ark*p->dt);
        }
        ++n;
    }
}

void cavitation::start(lexer *p, fdm *a, ghostcell *pgc)
{
    // vapour source from the final pressure: the source whose divergence the last projection imposed
    LOOP
    Sv(i,j,k) = Cs(i,j,k)*(pv - a->press(i,j,k));

    // divergence of the convection scheme: the scheme applied to a uniform field
    {
    field4 one(p);
    LOOP
    one(i,j,k) = 1.0;
    pgc->start4(p,one,gcval_alpha);
    LOOP
    a->L(i,j,k) = 0.0;
    pconvec->start(p,a,one,4,a->u,a->v,a->w);
    LOOP
    Dv(i,j,k) = -a->L(i,j,k);
    }

    field4 ark1(p),ark2(p);

    transport_step(p,a,pgc,alpha,ark1,1.0,0.0,1.0);
    transport_step(p,a,pgc,ark1,ark2,0.75,0.25,0.25);
    transport_step(p,a,pgc,ark2,alpha,1.0/3.0,2.0/3.0,2.0/3.0);

    properties(p,a,pgc);
    log(p,a,pgc);
}

void cavitation::transport_step(lexer *p, fdm *a, ghostcell *pgc, field &f, field &g, double wa, double wf, double wl)
{
    // SSP RK3 stage: g = wa alpha + wf f + wl dt (L(f) + f div u + Sv), alpha the value at the start of the step
    LOOP
    a->L(i,j,k) = 0.0;

    pconvec->start(p,a,f,4,a->u,a->v,a->w);

    LOOP
    {
        double val = (wf>0.0 ? wa*alpha(i,j,k) + wf*f(i,j,k) : wa*f(i,j,k))
                   + wl*p->dt*(a->L(i,j,k) + f(i,j,k)*Dv(i,j,k) + Sv(i,j,k));

        if(p->X10==1 && a->fbh4(i,j,k)>0.5)
        val = 0.0;

        g(i,j,k) = std::min(std::max(val,0.0),1.0);
    }
    pgc->start4(p,g,gcval_alpha);
}

void cavitation::properties(lexer *p, fdm *a, ghostcell *pgc)
{
    LOOP
    {
        const double al = alpha(i,j,k);
        const double ro = rho_mix(p,al);
        a->ro(i,j,k) = ro;
        a->visc(i,j,k) = ((1.0-al)*rhol*nul + al*rhov*nuv)/ro;
    }
    pgc->start4(p,a->ro,1);
    pgc->start4(p,a->visc,1);
}

void cavitation::log(lexer *p, fdm *a, ghostcell *pgc)
{
    double V=0.0, S=0.0, pmin=1.0e20;

    LOOP
    {
        const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        V += alpha(i,j,k)*dV;
        S += Cs(i,j,k)*(pv - a->press(i,j,k))*dV;
        pmin = std::min(pmin,a->press(i,j,k));
    }
    V = pgc->globalsum(V);
    S = pgc->globalsum(S);
    pmin = pgc->globalmin(pmin);

    if(p->mpirank==0)
    {
        const double R = std::cbrt(3.0*V/(4.0*M_PI));
        out<<std::setprecision(10)<<p->simtime+p->dt<<" \t "<<V<<" \t "<<R<<" \t "<<S<<" \t "<<pmin<<std::endl;

        if(p->count%p->P12==0)
        std::cout<<"cavitation: vapour volume "<<V<<" R_eq "<<R<<" production "<<S<<" p_min "<<pmin<<std::endl;
    }
}

void cavitation::name_ParaView_parallel(lexer *p, std::ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"alpha_v\"/>\n";
}

void cavitation::name_ParaView(lexer *p, std::ostream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"alpha_v\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void cavitation::offset_ParaView(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
    ++n;
}

void cavitation::print_3D(lexer *p, fdm *a, ghostcell *pgc, std::vector<char> &buffer, size_t &m)
{
    int iin;
    float ffn;

    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPLOOP
    {
    ffn=float(p->ipol4(alpha));
    std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
    }
}
