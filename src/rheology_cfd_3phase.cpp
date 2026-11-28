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

#include"rheology_cfd_3phase.h"
#include"rheology_core.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"
#include"ckpt_registry.h"
#include"heaviside.h"
#include"interface_width.h"
#include<cmath>
#include<cstring>
#include<algorithm>

rheology_cfd_3phase *rheology_cfd_3phase::self = nullptr;

namespace
{
// ln(cosh(x)) without overflow
inline double lncosh(double x)
{
    x = std::fabs(x);
    return x + std::log1p(std::exp(-2.0*x)) - 0.6931471805599453;
}

// point (x,z) in the triangle (x1,z1) (x2,z2) (x3,z3), either orientation
inline bool in_triangle(double x, double z, double x1, double z1, double x2, double z2, double x3, double z3)
{
    const double d1 = (x - x2)*(z1 - z2) - (x1 - x2)*(z - z2);
    const double d2 = (x - x3)*(z2 - z3) - (x2 - x3)*(z - z3);
    const double d3 = (x - x1)*(z3 - z1) - (x3 - x1)*(z - z1);
    const bool neg = (d1<0.0) || (d2<0.0) || (d3<0.0);
    const bool pos = (d1>0.0) || (d2>0.0) || (d3>0.0);
    return !(neg && pos);
}
}

rheology_cfd_3phase::rheology_cfd_3phase(lexer *p, fdm *a, ghostcell *pgc) : strain(p),
    alpha(p), alphas(p), pband(p,p->W131==3), sat(p,p->W136_S0<1.0), dsrc(p,p->W136_S0<1.0), imb(p->W136_S0<1.0), alpha0(p), flux(p), phyd(p),
    pe(p,p->W143==2), phis(p,p->W143==2 && p->W145_phi0>0.0), tpsi(p,p->W143==2 && p->W145_phi0>0.0),
    qpe(p,p->W143==2), qphi(p,p->W143==2 && p->W145_phi0>0.0),
    fpe(p,p->W143==2), fphi(p,p->W143==2 && p->W145_phi0>0.0),
    qsat(p,p->W136_S0<1.0), fsat(p,p->W136_S0<1.0),
    pore(p->W143), dil(p->W143==2 && p->W145_phi0>0.0), beta_thinc(3.5)
{
    self = this;

    rho_s = p->W141;
    rho_w = p->W1;
    rho_a = p->W3;
    eta_w = p->W1*p->W2;

    lambda0 = p->W144_lambda;
    kperm = p->W144_k;
    beta = MAX(p->W144_beta,1.0e-30);
    phi0 = p->W145_phi0;
    K3 = p->W145_K3;
    phic = p->W145_phic;
    dgrain = p->W124_d;

    pe_ini = (pore!=2 || lambda0==0.0);

    S0 = MIN(MAX(p->W136_S0,0.0),1.0);
    vimb = p->W136_v;
    nporo = MAX(1.0 - p->W124_phis,0.0);
    imbibed = 0.0;
    sweeps = 0;
    volume = volume0 = 0.0;

    eta.resize(size_t(MAX(p->gknox,1))*size_t(MAX(p->gknoy,1)));

    ckstate[0] = ckstate[1] = ckstate[2] = ckstate[3] = 0.0;
    ckid = ckpt_registry::add_raw(p,ckstate,4,false);

    if(p->mpirank==0)
    {
    cout<<"CFD three-phase slide (W 140): rho_s "<<rho_s<<", grains "<<p->W142<<" kg/m3, water "<<rho_w<<", air "<<rho_a
        <<"; mixture viscosity "<<(p->W149==2 ? "harmonic" : "arithmetic")<<"; pore pressure "
        <<(pore==0 ? "none" : (pore==1 ? "drained hydrostatic" : "hydrostatic + excess p_e"))<<endl;

        if(imb)
        cout<<"CFD three-phase slide: imbibition from saturation "<<S0<<" (porosity "<<nporo<<", saturated bulk density "
            <<rho_s + (1.0-S0)*nporo*rho_w<<"), front "<<(vimb>0.0 ? vimb : 0.0)<<(vimb>0.0 ? " m/s" : " (at once)")<<endl;

        if(pore==2)
        {
        cout<<"CFD three-phase slide: p_e lambda_0 "<<lambda0<<", k "<<(kperm>0.0 ? kperm : -1.0)<<(kperm>0.0 ? " m2" : " (Carman-Kozeny)")
            <<", beta "<<beta<<" 1/Pa";
        if(dil)
        cout<<"; dilatancy phi_0 "<<phi0<<", K3 "<<K3<<", phi_c "<<phic;
        cout<<endl;
        }
    }
}

rheology_cfd_3phase::~rheology_cfd_3phase()
{
    ckpt_registry::remove(ckid);

    if(self==this)
    self = nullptr;
}

void rheology_cfd_3phase::bc(lexer *p, ghostcell *pgc, field &f)
{
    pgc->start4(p,f,1);   // zero gradient at all boundaries, MPI and periodic exchange
}

void rheology_cfd_3phase::smooth(lexer *p, ghostcell *pgc)
{
    // two passes of (1/4 1/2 1/4) in each direction; a solid neighbour counts with the value of the cell
    const int size = p->imax*p->jmax*p->kmax;
    std::memcpy(alphas.V,alpha.V,sizeof(double)*size);

    for(int pass=0; pass<2; ++pass)
    for(int dir=0; dir<3; ++dir)
    {
        if(dir==1 && p->j_dir==0)
        continue;

        LOOP
        {
            const int di = dir==0, dj = dir==1, dk = dir==2;
            const double c = alphas(i,j,k);

            const int nm = (i-di-p->imin)*p->jmax*p->kmax + (j-dj-p->jmin)*p->kmax + k-dk-p->kmin;
            const int np = (i+di-p->imin)*p->jmax*p->kmax + (j+dj-p->jmin)*p->kmax + k+dk-p->kmin;

            const double am = p->flag4[nm]>0 ? alphas(i-di,j-dj,k-dk) : c;
            const double ap = p->flag4[np]>0 ? alphas(i+di,j+dj,k+dk) : c;

            flux(i,j,k) = 0.25*am + 0.5*c + 0.25*ap;
        }

        LOOP
        alphas(i,j,k) = flux(i,j,k);

        bc(p,pgc,alphas);
    }
}

void rheology_cfd_3phase::band(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->W131!=3)
    return;

    // the buoyant weight of the smoothed band above the cell: from the top of the column of the rank down,
    // sum of (rho_s - rho_f) alpha_s dz over the contiguous band cells (1e-3 < alpha_s < 1 - 1e-3), half of the
    // own cell; none below the band (the solved pressure is larger there), as rho_1 (g.n) int H ds of the
    // two-phase W 131 3. The band is 4 - 5 cells thick: a rank
    // boundary inside it starts the sum again (the floor is then lower, not wrong in sign)
    const double g = std::fabs(p->W22);
    const double eps = 1.0e-3;

    for(i=0; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    {
        double sum = 0.0;

        for(k=p->knoz-1; k>=0; --k)
        {
            pband(i,j,k) = 0.0;

            if(p->flag4[IJK]<=0)
            {
            sum = 0.0;
            continue;
            }

            const double al = alphas(i,j,k);

            // only the band: below it (alpha_s = 1) the solved pressure carries the overburden, and a floor there
            // would hold p' at the lithostatic value where the pore pressure fluidises the slide
            if(al<=eps || al>=1.0-eps)
            {
            sum = 0.0;
            continue;
            }

            const double H = heaviside(a->phi(i,j,k),interface_width(p,a->phi,i,j,k));
            const double drho = MAX(rho_slide(i,j,k) - (H*rho_w + (1.0-H)*rho_a),0.0);
            const double w = drho*g*al*p->DZN[KP];

            pband(i,j,k) = sum + 0.5*w;
            sum += w;
        }
    }

    bc(p,pgc,pband);
}

void rheology_cfd_3phase::imbibe(lexer *p, fdm *a, ghostcell *pgc)
{
    // a cell of the rank, water (outside the slide, below the level set surface)
    auto inside = [&](int ii, int jj, int kk) -> bool
    {
        return ii>=0 && ii<p->knox && jj>=0 && jj<p->knoy && kk>=0 && kk<p->knoz;
    };

    auto fl4 = [&](int ii, int jj, int kk) -> int
    {
        return p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin];
    };

    auto water = [&](int ii, int jj, int kk) -> bool
    {
        return fl4(ii,jj,kk)>0 && alpha(ii,jj,kk)<0.5 && a->phi(ii,jj,kk)>0.0;
    };

    const int nb = p->j_dir==1 ? 6 : 4;
    const int di[6] = {-1,1,0,0,0,0}, dk[6] = {0,0,-1,1,0,0}, dj[6] = {0,0,0,0,-1,1};

    LOOP
    dsrc(i,j,k) = 0.0;

    // the saturation of the step: the cells that take up water and the water they take, given to their water
    // neighbours on the rank (as a sink rate per cell volume), the rest into the pool of all contact cells
    double Qtot = 0.0, Qpool = 0.0;
    std::vector<double> dS(size_t(p->knox)*p->knoy*p->knoz,0.0);
    int q = 0;

    LOOP
    {
        q = (i*p->knoy + j)*p->knoz + k;

        if(alpha(i,j,k)<=0.5 || sat(i,j,k)>=1.0 || phyd(i,j,k)<=0.0)
        continue;

        bool wet = false;

        for(int n=0; n<nb && !wet; ++n)
        {
            const int ii=i+di[n], jj=j+dj[n], kk=k+dk[n];

            if(fl4(ii,jj,kk)<=0)
            continue;

            // a full neighbour, S >= 0.99 (with S > 0.5 the front ran at twice the speed v; S = 1 exactly stalled it,
            // the transport at the surface of the slide mixes a little S into the cells)
            wet = water(ii,jj,kk) || (alpha(ii,jj,kk)>0.5 && sat(ii,jj,kk)>=0.99);
        }

        if(!wet)
        continue;

        double h = MIN(p->DXN[IP],p->DZN[KP]);
        if(p->j_dir==1)
        h = MIN(h,p->DYN[JP]);

        dS[q] = vimb>0.0 ? MIN(1.0-sat(i,j,k),p->dt*vimb/h) : 1.0-sat(i,j,k);
    }

    LOOP
    {
        q = (i*p->knoy + j)*p->knoz + k;

        if(dS[q]<=0.0)
        continue;

        const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        const double Q = nporo*alpha(i,j,k)*dS[q]*dV/p->dt;   // water taken up [m3/s]

        sat(i,j,k) += dS[q];
        Qtot += Q;

        int nw = 0;
        for(int n=0; n<nb; ++n)
        if(inside(i+di[n],j+dj[n],k+dk[n]) && water(i+di[n],j+dj[n],k+dk[n]))
        ++nw;

        if(nw==0)
        {
        Qpool += Q;
        continue;
        }

        for(int n=0; n<nb; ++n)
        {
            const int ii=i+di[n], jj=j+dj[n], kk=k+dk[n];

            if(inside(ii,jj,kk) && water(ii,jj,kk))
            dsrc(ii,jj,kk) -= Q/double(nw)/(p->DXN[ii+marge]*p->DYN[jj+marge]*p->DZN[kk+marge]);
        }
    }

    Qtot = pgc->globalsum(Qtot);
    Qpool = pgc->globalsum(Qpool);

    // the pool: into all water cells next to the slide
    if(Qpool>0.0)
    {
        double Vc = 0.0;

        LOOP
        if(water(i,j,k))
        for(int n=0; n<nb; ++n)
        if(fl4(i+di[n],j+dj[n],k+dk[n])>0 && alpha(i+di[n],j+dj[n],k+dk[n])>0.5)
        {
        Vc += p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        break;
        }

        Vc = pgc->globalsum(Vc);

        if(Vc>0.0)
        LOOP
        if(water(i,j,k))
        for(int n=0; n<nb; ++n)
        if(fl4(i+di[n],j+dj[n],k+dk[n])>0 && alpha(i+di[n],j+dj[n],k+dk[n])>0.5)
        {
        dsrc(i,j,k) -= Qpool/Vc;
        break;
        }
    }

    // the pore air that escapes: a source in all air cells
    if(Qtot>0.0)
    {
        double Va = 0.0;

        LOOP
        if(alpha(i,j,k)<0.5 && a->phi(i,j,k)<0.0)
        Va += p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

        Va = pgc->globalsum(Va);

        if(Va>0.0)
        LOOP
        if(alpha(i,j,k)<0.5 && a->phi(i,j,k)<0.0)
        dsrc(i,j,k) += Qtot/Va;
    }

    imbibed += Qtot*p->dt;

    bc(p,pgc,sat);
    bc(p,pgc,dsrc);
}

void rheology_cfd_3phase::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    // sub-cell sampling of the boxes and triangles: the initial volume to about 1/ns^2 of a cell at the edges
    const int ns = 10;
    const int nsy = p->j_dir==1 ? ns : 1;

    LOOP
    {
        int in=0, tot=0;

        for(int qi=0; qi<ns; ++qi)
        for(int qj=0; qj<nsy; ++qj)
        for(int qk=0; qk<ns; ++qk)
        {
            const double x = p->XN[IP] + (qi+0.5)/double(ns)*p->DXN[IP];
            const double y = p->YN[JP] + (qj+0.5)/double(nsy)*p->DYN[JP];
            const double z = p->ZN[KP] + (qk+0.5)/double(ns)*p->DZN[KP];

            bool inside = false;

            for(int q=0; q<p->W147 && !inside; ++q)
            if(x>=p->W147_xs[q] && x<=p->W147_xe[q] && z>=p->W147_zs[q] && z<=p->W147_ze[q]
               && (p->j_dir==0 || (y>=p->W147_ys[q] && y<=p->W147_ye[q])))
            inside = true;

            for(int q=0; q<p->W148 && !inside; ++q)
            if(in_triangle(x,z,p->W148_x1[q],p->W148_z1[q],p->W148_x2[q],p->W148_z2[q],p->W148_x3[q],p->W148_z3[q]))
            inside = true;

            in += inside ? 1 : 0;
            ++tot;
        }

        alpha(i,j,k) = double(in)/double(tot);
    }

    bc(p,pgc,alpha);

    if(imb)
    {
        LOOP
        {
        sat(i,j,k) = S0;
        dsrc(i,j,k) = 0.0;
        }
        bc(p,pgc,sat);
        bc(p,pgc,dsrc);
    }

    smooth(p,pgc);
    band(p,a,pgc);

    if(pore==2)
    {
        LOOP
        pe(i,j,k) = 0.0;
        bc(p,pgc,pe);
    }

    if(dil)
    {
        LOOP
        {
        phis(i,j,k) = phi0;
        tpsi(i,j,k) = 0.0;
        }
        bc(p,pgc,phis);
        bc(p,pgc,tpsi);
    }

    volume = 0.0;
    LOOP
    volume += alpha(i,j,k)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
    volume = pgc->globalsum(volume);
    volume0 = volume;

    hydrostatic(p,a,pgc);

    if(p->mpirank==0)
    cout<<"CFD three-phase slide: initial volume "<<setprecision(10)<<volume0<<setprecision(6)
        <<" m3 ("<<p->W147<<" boxes W 147, "<<p->W148<<" triangles W 148)"<<endl;
}

bool rheology_cfd_3phase::face_open(lexer *p, int dir, int ii, int jj, int kk) const
{
    // the face between (ii,jj,kk) and its neighbour in direction dir: both cells fluid
    const int n0 = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
    const int n1 = n0 + (dir==0 ? p->jmax*p->kmax : (dir==1 ? p->kmax : 1));

    return p->flag4[n0]>0 && p->flag4[n1]>0;
}

double rheology_cfd_3phase::face_vel(lexer *p, fdm *a, int dir, int ii, int jj, int kk) const
{
    if(!face_open(p,dir,ii,jj,kk))
    return 0.0;

    if(dir==0)
    return a->u(ii,jj,kk);

    if(dir==1)
    return a->v(ii,jj,kk);

    return a->w(ii,jj,kk);
}

// fraction (times the cell length) of the upwind cell that crosses its downstream face in a step of
// Courant number c, from the THINC profile 0.5 (1 + g tanh(beta (xi - s))), xi = 0 at the upstream face
double rheology_cfd_3phase::thinc(double aback, double aup, double aahead, double c) const
{
    const double eps = 1.0e-8;

    if(c<=0.0)
    return 0.0;

    if(aup<eps || aup>1.0-eps || (aahead-aup)*(aup-aback)<=0.0)
    return aup*c;

    const double g = aahead>=aback ? 1.0 : -1.0;
    const double B = beta_thinc;
    const double Q = (2.0*aup - 1.0)*B*g;

    // the jump position s from the cell mean: e^{2 B s} = (e^B - e^Q)/(e^Q - e^-B)
    const double X2 = (std::exp(B) - std::exp(Q))/(std::exp(Q) - std::exp(-B));

    if(!(X2>0.0))
    return aup*c;

    const double s = std::log(X2)/(2.0*B);

    double amount = 0.5*c + g/(2.0*B)*(lncosh(B*(1.0-s)) - lncosh(B*(1.0-c-s)));

    return std::min(std::max(amount,0.0),std::min(c,aup));
}

void rheology_cfd_3phase::sweep(lexer *p, fdm *a, ghostcell *pgc, int dir, double dt)
{
    const bool q1 = (pore==2);
    const bool q2 = dil;
    const bool q3 = imb;

    // fluxes through the faces +dir of the cells -1 .. kno-1 in dir (the faces of the rank's cells)
    const int i0 = dir==0 ? -1 : 0;
    const int j0 = dir==1 ? -1 : 0;
    const int k0 = dir==2 ? -1 : 0;

    for(i=i0; i<p->knox; ++i)
    for(j=j0; j<p->knoy; ++j)
    for(k=k0; k<p->knoz; ++k)
    {
        const double uf = face_vel(p,a,dir,i,j,k);

        flux(i,j,k) = 0.0;

        if(q1)
        fpe(i,j,k) = 0.0;

        if(q2)
        fphi(i,j,k) = 0.0;

        if(q3)
        fsat(i,j,k) = 0.0;

        if(uf==0.0)
        continue;

        const int di = dir==0, dj = dir==1, dk = dir==2;

        int ui, uj, uk;
        double ab, aup, aah, dx;

        if(uf>0.0)
        {
        ui=i; uj=j; uk=k;
        ab  = alpha(i-di,j-dj,k-dk);
        aup = alpha(i,j,k);
        aah = alpha(i+di,j+dj,k+dk);
        dx = dir==0 ? p->DXN[IP] : (dir==1 ? p->DYN[JP] : p->DZN[KP]);
        }
        else
        {
        ui=i+di; uj=j+dj; uk=k+dk;
        ab  = alpha(i+2*di,j+2*dj,k+2*dk);
        aup = alpha(ui,uj,uk);
        aah = alpha(i,j,k);
        dx = dir==0 ? p->DXN[IP1] : (dir==1 ? p->DYN[JP1] : p->DZN[KP1]);
        }

        const double c = std::fabs(uf)*dt/dx;
        const double F = (uf>0.0 ? 1.0 : -1.0)*thinc(ab,aup,aah,c)*dx;

        flux(i,j,k) = F;

        // per unit volume of slide: p_e, phi and the saturation go with the slide of the upwind cell
        if(q1 || q2 || q3)
        {
            const double au = alpha(ui,uj,uk);

            if(au>1.0e-12)
            {
                if(q1)
                fpe(i,j,k) = F*qpe(ui,uj,uk)/au;

                if(q2)
                fphi(i,j,k) = F*qphi(ui,uj,uk)/au;

                if(q3)
                fsat(i,j,k) = F*qsat(ui,uj,uk)/au;
            }
        }
    }

    // update: -div(F) + c div(u) (Weymouth & Yue 2010), c = 1 in cells with alpha > 1/2 at the start of the step
    double clipped = 0.0;

    LOOP
    {
        const int di = dir==0, dj = dir==1, dk = dir==2;
        const double dx = dir==0 ? p->DXN[IP] : (dir==1 ? p->DYN[JP] : p->DZN[KP]);
        const double du = (face_vel(p,a,dir,i,j,k) - face_vel(p,a,dir,i-di,j-dj,k-dk))*dt/dx;
        const double cc = alpha0(i,j,k)>0.5 ? 1.0 : 0.0;

        const double aold = alpha(i,j,k);
        double an = aold - (flux(i,j,k) - flux(i-di,j-dj,k-dk))/dx + cc*du;

        if(q1 || q2 || q3)
        {
            const double r1 = q1 && aold>1.0e-12 ? qpe(i,j,k)/aold : 0.0;
            const double r2 = q2 && aold>1.0e-12 ? qphi(i,j,k)/aold : 0.0;
            const double r3 = q3 && aold>1.0e-12 ? qsat(i,j,k)/aold : 0.0;

            if(q3)
            qsat(i,j,k) += -(fsat(i,j,k) - fsat(i-di,j-dj,k-dk))/dx + cc*r3*du;

            if(q1)
            qpe(i,j,k) += -(fpe(i,j,k) - fpe(i-di,j-dj,k-dk))/dx + cc*r1*du;

            if(q2)
            qphi(i,j,k) += -(fphi(i,j,k) - fphi(i-di,j-dj,k-dk))/dx + cc*r2*du;
        }

        if(an<0.0 || an>1.0)
        {
        const double ac = std::min(std::max(an,0.0),1.0);
        clipped += std::fabs(an-ac)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        an = ac;
        }

        alpha(i,j,k) = an;

        if(an<1.0e-12)
        {
            if(q1)
            qpe(i,j,k) = 0.0;

            if(q2)
            qphi(i,j,k) = 0.0;

            if(q3)
            qsat(i,j,k) = 0.0;
        }
    }

    bc(p,pgc,alpha);

    if(q1)
    bc(p,pgc,qpe);

    if(q2)
    bc(p,pgc,qphi);

    if(q3)
    bc(p,pgc,qsat);

    clipsum += clipped;
}

void rheology_cfd_3phase::start(lexer *p, fdm *a, ghostcell *pgc, solver *psolv)
{
    // the scalars as the last step (or the checkpoint of a restart) left them
    sweeps = int(ckstate[0]);
    pe_ini = pe_ini || ckstate[1]>0.5;
    if(ckstate[2]>0.0)
    volume0 = ckstate[2];
    imbibed = ckstate[3];

    const int size = p->imax*p->jmax*p->kmax;
    const int ndir = p->j_dir==1 ? 3 : 2;
    const int order3[3] = {0,1,2};
    const int order2[2] = {0,2};

    // Courant number per direction: sub-steps above 0.4 (the split scheme is bounded up to 1/2)
    double cmax = 0.0;
    LOOP
    {
    cmax = MAX(cmax,std::fabs(a->u(i,j,k))*p->dt/p->DXN[IP]);
    cmax = MAX(cmax,std::fabs(a->u(i-1,j,k))*p->dt/p->DXN[IP]);
    if(p->j_dir==1)
    {
    cmax = MAX(cmax,std::fabs(a->v(i,j,k))*p->dt/p->DYN[JP]);
    cmax = MAX(cmax,std::fabs(a->v(i,j-1,k))*p->dt/p->DYN[JP]);
    }
    cmax = MAX(cmax,std::fabs(a->w(i,j,k))*p->dt/p->DZN[KP]);
    cmax = MAX(cmax,std::fabs(a->w(i,j,k-1))*p->dt/p->DZN[KP]);
    }
    cmax = pgc->globalmax(cmax);

    const int nsub = MAX(1,int(std::ceil(cmax/0.4)));
    const double dts = p->dt/double(nsub);

    // p_e and phi per unit volume of the mixture during the transport
    if(pore==2)
    for(int n=0; n<size; ++n)
    qpe.V[n] = alpha.V[n]*pe.V[n];

    if(dil)
    for(int n=0; n<size; ++n)
    qphi.V[n] = alpha.V[n]*phis.V[n];

    if(imb)
    for(int n=0; n<size; ++n)
    qsat.V[n] = alpha.V[n]*sat.V[n];

    clipsum = 0.0;

    // the indicator c of the divergence correction from the start of the step, for all sub-steps: with the sources of
    // the imbibition (div(u) != 0 in the water and the air) a cell that crossed alpha = 1/2 in a sub-step took a
    // share of the source and the slide volume changed (1e-3 in case 15)
    std::memcpy(alpha0.V,alpha.V,sizeof(double)*size);

    for(int s=0; s<nsub; ++s)
    {

        for(int d=0; d<ndir; ++d)
        {
            int dir = ndir==3 ? order3[d] : order2[d];

            // the order of the directions alternates from step to step (x y z, z y x)
            if(sweeps%2==1)
            dir = ndir==3 ? order3[2-d] : order2[1-d];

            sweep(p,a,pgc,dir,dts);
        }

        ++sweeps;
    }

    if(pore==2)
    {
        LOOP
        pe(i,j,k) = alpha(i,j,k)>1.0e-12 ? qpe(i,j,k)/alpha(i,j,k) : 0.0;
        bc(p,pgc,pe);
    }

    if(dil)
    {
        LOOP
        phis(i,j,k) = alpha(i,j,k)>1.0e-12 ? qphi(i,j,k)/alpha(i,j,k) : phi0;
        bc(p,pgc,phis);
    }

    if(imb)
    {
        LOOP
        sat(i,j,k) = alpha(i,j,k)>1.0e-12 ? MIN(MAX(qsat(i,j,k)/alpha(i,j,k),0.0),1.0) : S0;
        bc(p,pgc,sat);
    }

    double clipped = pgc->globalsum(clipsum);

    smooth(p,pgc);
    band(p,a,pgc);

    hydrostatic(p,a,pgc);

    // the pores below the water surface take up water (W 136); the band and the density of the step with it
    if(imb)
    {
    imbibe(p,a,pgc);
    band(p,a,pgc);
    }

    if(pore==2)
    excess(p,a,pgc,psolv);

    diagnostics(p,a,pgc,clipped,nsub);

    ckstate[0] = double(sweeps);
    ckstate[1] = pe_ini ? 1.0 : 0.0;
    ckstate[2] = volume0;
    ckstate[3] = imbibed;
}

void rheology_cfd_3phase::hydrostatic(lexer *p, fdm *a, ghostcell *pgc)
{
    if(pore==0)
    {
        LOOP
        phyd(i,j,k) = 0.0;
        bc(p,pgc,phyd);
        return;
    }

    // water surface of the column: the highest z + phi of a water cell outside the slide and above the slide of
    // the column (phi <= the distance to the surface above, = it at the surface). Water under the slide (a film
    // the slide rides on) does not count: its phi is the distance to the slide above, and the slide over it
    // counted as dry
    const int gny = MAX(p->gknoy,1);
    std::fill(eta.begin(),eta.end(),-1.0e30);

    std::vector<double> zs(eta.size(),-1.0e30);

    LOOP
    if(alpha(i,j,k)>0.5)
    {
    const size_t q = size_t(i+p->origin_i)*size_t(gny) + size_t(p->j_dir==1 ? j+p->origin_j : 0);
    zs[q] = MAX(zs[q], p->ZP[KP]);
    }

    pgc->globalmax(zs.data(),int(zs.size()));

    LOOP
    if(alpha(i,j,k)<0.5 && a->phi(i,j,k)>0.0)
    {
    const size_t q = size_t(i+p->origin_i)*size_t(gny) + size_t(p->j_dir==1 ? j+p->origin_j : 0);

    if(p->ZP[KP]>zs[q])
    eta[q] = MAX(eta[q], p->ZP[KP] + a->phi(i,j,k));
    }

    pgc->globalmax(eta.data(),int(eta.size()));

    // columns the slide fills (no open water left in them: a slide that has pushed the water aside) take the water
    // surface of the next column with water, spreading through the columns with slide only: the slide below the
    // surrounding water level is submerged (its drained pore pressure, its imbibition W 136); before, these columns
    // had no water surface and the slide in them counted as dry
    {
        const int gnx = MAX(p->gknox,1);
        std::vector<double> sc(eta.size(),0.0);

        for(size_t q=0; q<sc.size(); ++q)
        sc[q] = zs[q]>-1.0e29 ? 1.0 : 0.0;

        std::vector<double> en(eta);
        bool changed = true;

        for(int it=0; it<gnx+gny && changed; ++it)
        {
            changed = false;

            for(int ii=0; ii<gnx; ++ii)
            for(int jj=0; jj<gny; ++jj)
            {
                const size_t q = size_t(ii)*size_t(gny) + size_t(jj);

                if(eta[q]>-1.0e29 || sc[q]<0.5)
                continue;

                double e = -1.0e30;

                if(ii>0)
                e = MAX(e,eta[q-gny]);
                if(ii<gnx-1)
                e = MAX(e,eta[q+gny]);
                if(jj>0)
                e = MAX(e,eta[q-1]);
                if(jj<gny-1)
                e = MAX(e,eta[q+1]);

                if(e>-1.0e29)
                {
                en[q] = e;
                changed = true;
                }
            }

            eta = en;
        }
    }

    const double g = std::fabs(p->W22);

    LOOP
    {
    const size_t q = size_t(i+p->origin_i)*size_t(gny) + size_t(p->j_dir==1 ? j+p->origin_j : 0);
    const double e = eta[q];
    phyd(i,j,k) = (e>-1.0e29 && e>p->ZP[KP]) ? rho_w*g*(e - p->ZP[KP]) : 0.0;
    }

    bc(p,pgc,phyd);
}

double rheology_cfd_3phase::permeability(double phi) const
{
    if(kperm>0.0)
    return kperm;

    // Carman-Kozeny
    phi = MIN(MAX(phi,0.3),0.7);
    return dgrain*dgrain*(1.0-phi)*(1.0-phi)*(1.0-phi)/(180.0*phi*phi);
}

void rheology_cfd_3phase::excess(lexer *p, fdm *a, ghostcell *pgc, solver *psolv)
{
    // initial excess pore pressure after the first step: lambda_0 (p - p_hyd) of the slide at rest
    if(!pe_ini && p->count>=2)
    {
        LOOP
        pe(i,j,k) = alpha(i,j,k)>0.5 ? lambda0*MAX(a->press(i,j,k) - p->pressgage - ppore(i,j,k),0.0) : 0.0;
        bc(p,pgc,pe);
        pe_ini = true;

        if(p->mpirank==0)
        cout<<"CFD three-phase slide: initial excess pore pressure lambda_0 "<<lambda0<<" (p - p_hyd)"<<endl;
    }

    // dilatancy: tan(psi) of the step and the source gamma tan(psi) [1/s], from the second step on (the first
    // step starts from the initial pressure, not the solved one). The viscous number takes p' of at least the
    // buoyant weight of one cell: at the surface of the slide p' -> 0, J -> infinity and phi_eq -> 0, and with a
    // small compressibility beta the dilation of the surface cells sucked the whole loose column to p_e < 0
    if(dil)
    {
        const double drho = std::fabs(rho_s - rho_w);
        const double g = std::fabs(p->W22);

        LOOP
        {
            double tp = 0.0;

            if(alpha(i,j,k)>0.5 && p->count>=2)
            {
                const double pp = MAX(peff(p,a),drho*g*p->DZN[KP]);
                const double gam = strainterm(p,a->u,a->v,a->w);
                const double J = rheo::viscous_number(gam,pp,eta_w);
                tp = MIN(MAX(rheo::tan_psi(K3,phis(i,j,k),phic,J),-0.6),0.6);
            }

            tpsi(i,j,k) = tp;
        }
        bc(p,pgc,tpsi);
    }

    // beta (p_e - p_e*)/dt - div(K grad p_e) = - gamma tan(psi),  K = k/eta_w;  p_e = 0 outside the slide
    a->M.reset();

    const double idt = beta/p->dt;
    int nn=0;

    LOOP
    a->rhsvec.V[nn++] = 0.0;

    auto Kc = [&](int ii, int jj, int kk) -> double
    {
        const double ph = dil ? phis(ii,jj,kk) : (phi0>0.0 ? phi0 : 0.6);
        return permeability(ph)/eta_w;
    };

    auto Kf = [&](double K1, double K2) -> double
    {
        return (K1>0.0 && K2>0.0) ? 2.0*K1*K2/(K1+K2) : 0.0;
    };

    nn=0;
    LOOP
    {
        if(alpha(i,j,k)<=0.5)
        {
        a->M.p[nn] = 1.0;
        a->M.s[nn] = a->M.n[nn] = a->M.e[nn] = a->M.w[nn] = a->M.b[nn] = a->M.t[nn] = 0.0;
        a->rhsvec.V[nn] = 0.0;
        ++nn;
        continue;
        }

        const double K0 = Kc(i,j,k);

        // face coefficients K_f/(dx_cell dx_centres); walls: no flux; outside the slide: p_e = 0
        double cs = Kf(K0,Kc(i-1,j,k))/(p->DXN[IP]*p->DXP[IM1]);
        double cn = Kf(K0,Kc(i+1,j,k))/(p->DXN[IP]*p->DXP[IP]);
        double ce = p->j_dir==1 ? Kf(K0,Kc(i,j-1,k))/(p->DYN[JP]*p->DYP[JM1]) : 0.0;
        double cw = p->j_dir==1 ? Kf(K0,Kc(i,j+1,k))/(p->DYN[JP]*p->DYP[JP]) : 0.0;
        double cb = Kf(K0,Kc(i,j,k-1))/(p->DZN[KP]*p->DZP[KM1]);
        double ct = Kf(K0,Kc(i,j,k+1))/(p->DZN[KP]*p->DZP[KP]);

        if(p->flag4[Im1JK]<0)
        cs = 0.0;
        if(p->flag4[Ip1JK]<0)
        cn = 0.0;
        if(p->flag4[IJm1K]<0)
        ce = 0.0;
        if(p->flag4[IJp1K]<0)
        cw = 0.0;
        if(p->flag4[IJKm1]<0)
        cb = 0.0;
        if(p->flag4[IJKp1]<0)
        ct = 0.0;

        a->M.p[nn] = idt + cs + cn + ce + cw + cb + ct;

        // the neighbours outside the slide hold p_e = 0: their coupling drops out of the row (symmetric matrix)
        a->M.s[nn] = alpha(i-1,j,k)>0.5 ? -cs : 0.0;
        a->M.n[nn] = alpha(i+1,j,k)>0.5 ? -cn : 0.0;
        a->M.e[nn] = alpha(i,j-1,k)>0.5 ? -ce : 0.0;
        a->M.w[nn] = alpha(i,j+1,k)>0.5 ? -cw : 0.0;
        a->M.b[nn] = alpha(i,j,k-1)>0.5 ? -cb : 0.0;
        a->M.t[nn] = alpha(i,j,k+1)>0.5 ? -ct : 0.0;

        const double src = dil ? strainterm(p,a->u,a->v,a->w)*tpsi(i,j,k) : 0.0;

        a->rhsvec.V[nn] = idt*pe(i,j,k) - src;

        ++nn;
    }

    // the solver needs the old values for the dilation below
    if(dil)
    std::memcpy(flux.V,pe.V,sizeof(double)*p->imax*p->jmax*p->kmax);

    psolv->start(p,a,pgc,pe,a->rhsvec,4);

    // no tension in the skeleton: p_e at most the effective stress p - p_hyd (a contracting loose slide is fully
    // fluidised there, p' = 0); with a small beta the contraction alone drove p_e to 5 times the overburden
    LOOP
    if(alpha(i,j,k)>0.5)
    pe(i,j,k) = MIN(pe(i,j,k),MAX(a->press(i,j,k) - p->pressgage - ppore(i,j,k),0.0));

    bc(p,pgc,pe);

    if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
    cout<<"p_e iter: "<<p->solveriter<<endl;

    // dilation D_v = beta dp_e/dt + gamma tan(psi) (the pore fluid that flowed in), dphi/dt = -phi D_v
    if(dil)
    {
        LOOP
        if(alpha(i,j,k)>0.5)
        {
            const double Dv = idt*(pe(i,j,k) - flux(i,j,k)) + strainterm(p,a->u,a->v,a->w)*tpsi(i,j,k);
            phis(i,j,k) = MIN(MAX(phis(i,j,k)*(1.0 - p->dt*Dv),0.3),0.74);
        }
        bc(p,pgc,phis);
    }
}

double rheology_cfd_3phase::peff(lexer *p, fdm *a)
{
    double pp = a->press(i,j,k) - p->pressgage;

    // the pore water carries the hydrostatic pressure; with imbibition (W 136) only the saturated share S of the pores:
    // the air in the dry pores of a slide that stands out of the water is open to the atmosphere, so dry grains under
    // the water surface carry their full weight (no buoyancy until the water reaches them)
    if(pore>=1)
    pp -= ppore(i,j,k);

    if(pore==2)
    pp -= pe(i,j,k);

    // W 131 3: at least the weight of the interface band above the cell
    if(p->W131==3)
    pp = MAX(pp,pband(i,j,k));

    return pp;
}

double rheology_cfd_3phase::rho_face(lexer *p, int aa, int bb, int cc, double H)
{
    const double af = 0.5*(alphas(i,j,k) + alphas(i+aa,j+bb,k+cc));

    if(imb)
    return 0.5*(alphas(i,j,k)*rho_slide(i,j,k) + alphas(i+aa,j+bb,k+cc)*rho_slide(i+aa,j+bb,k+cc))
           + (1.0-af)*(H*rho_w + (1.0-H)*rho_a);

    return af*rho_s + (1.0-af)*(H*rho_w + (1.0-H)*rho_a);
}

void rheology_cfd_3phase::diagnostics(lexer *p, fdm *a, ghostcell *pgc, double clipped, int nsub)
{
    double amin=1.0e10, amax=-1.0e10;
    volume = 0.0;

    LOOP
    {
    volume += alpha(i,j,k)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
    amin = MIN(amin,alpha(i,j,k));
    amax = MAX(amax,alpha(i,j,k));
    }

    volume = pgc->globalsum(volume);
    amin = pgc->globalmin(amin);
    amax = pgc->globalmax(amax);

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"Slide volume: "<<setprecision(10)<<volume<<setprecision(6)<<" (rel. change "<<(volume0>0.0 ? volume/volume0-1.0 : 0.0)
        <<"), alpha_s "<<amin<<" - "<<amax<<", clipped "<<clipped<<" m3, sub-steps "<<nsub<<endl;

    if(imb)
    cout<<"Slide pore water taken up: "<<imbibed<<" m3"<<endl;
    }
}

int rheology_cfd_3phase::nvtu(lexer *p) const
{
    return 2 + (pore>=1 ? 1 : 0) + (pore==2 ? 1 : 0) + (dil ? 1 : 0) + (imb ? 1 : 0);
}

void rheology_cfd_3phase::offset_vtu(lexer *p, int *offset, int &n)
{
    for(int q=0; q<nvtu(p); ++q)
    {
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
    ++n;
    }
}

void rheology_cfd_3phase::name_vtu(lexer *p, std::ostream &result, int *offset, int &n)
{
    const char *names[6] = {"alpha_s","p_eff","p_hyd","p_e","phi_s","saturation"};
    const bool on[6] = {true,true,pore>=1,pore==2,dil,imb};

    for(int q=0; q<6; ++q)
    if(on[q])
    {
    result<<"<DataArray type=\"Float32\" Name=\""<<names[q]<<"\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
}

void rheology_cfd_3phase::name_pvtu(lexer *p, std::ostream &result)
{
    const char *names[6] = {"alpha_s","p_eff","p_hyd","p_e","phi_s","saturation"};
    const bool on[6] = {true,true,pore>=1,pore==2,dil,imb};

    for(int q=0; q<6; ++q)
    if(on[q])
    result<<"<PDataArray type=\"Float32\" Name=\""<<names[q]<<"\"/>\n";
}

void rheology_cfd_3phase::print_vtu(lexer *p, fdm *a, ghostcell *pgc, std::vector<char> &buffer, size_t &m)
{
    // p_eff into the scratch field (the flux of the last sweep is no longer needed)
    LOOP
    flux(i,j,k) = peff(p,a);
    bc(p,pgc,flux);

    field *f[6] = {&alpha,&flux,&phyd,(pore==2 ? &pe : nullptr),(dil ? &phis : nullptr),(imb ? &sat : nullptr)};
    const bool on[6] = {true,true,pore>=1,pore==2,dil,imb};

    for(int q=0; q<6; ++q)
    if(on[q])
    {
        int iin = sizeof(float)*p->pointnum;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m += sizeof(int);

        TPLOOP
        {
        float ffn = float(p->ipol4(*f[q]));
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m += sizeof(float);
        }
    }
}
