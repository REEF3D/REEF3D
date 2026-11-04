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

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include<vector>

// Phase 2 (F 81 2): the momentum of the water from the particles, the air Eulerian.
// step_begin, before the momentum step: the particle velocities go to the velocity faces in
// the water (phi at the face >= 0), weighted with the trilinear weights of the interpolation
// (the particles of the neighbours near the subdomain boundary come as a halo); faces in the water
// that no particle reaches and the faces in the air keep the CFD velocity. The momentum step then
// leaves out the convection on the faces in the water (mask_convection, from momentum_rk): the
// particles carry it. After the step (pic_f::start) the particles take the change of the water
// velocity (FLIP) and, by F 82, the velocity itself (PIC):
//     v_p = F82 (v_p + I(u_new - u_old)) + (1 - F82) I(u_new)
// with I the interpolation of the water velocity extended into the air (pic_f_extend.cpp)
// (Brackbill & Ruppel 1986, Zhu & Bridson 2005; Chen et al. 2018: F 82 = 0.96).

namespace
{
inline int fid(lexer *p, int i, int j, int k)
{
    return (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin;
}
}

void pic_f::step_begin(lexer *p, fdm *a, ghostcell *pgc)
{
    if(!flipmode)
    return;

    std::vector<double> H;
    halo(p,pgc,H,true);

    const int size = p->imax*p->jmax*p->kmax;
    std::vector<double> sw(size), sv(size);

    for(int c=0; c<3; ++c)
    {
        if(c==1 && p->j_dir==0)
        continue;

        field &f = c==0 ? static_cast<field&>(a->u) : (c==1 ? static_cast<field&>(a->v) : static_cast<field&>(a->w));
        const int *flag = c==0 ? p->flag1 : (c==1 ? p->flag2 : p->flag3);
        const int di = c==0, dj = c==1, dk = c==2;

        std::fill(sw.begin(),sw.end(),0.0);
        std::fill(sv.begin(),sv.end(),0.0);

        auto scatter = [&](double x, double y, double z, double val)
        {
            int ii[8], jj[8], kk[8];
            double ww[8];
            const int n = stencil(p,c+1,x,y,z,ii,jj,kk,ww);

            for(int m=0; m<n; ++m)
            {
                const int g = fid(p,ii[m],jj[m],kk[m]);
                sw[g] += ww[m];
                sv[g] += ww[m]*val;
            }
        };

        for(const auto &q : P)
        scatter(q.x,q.y,q.z, c==0 ? q.u : (c==1 ? q.v : q.w));

        for(size_t m=0; m+6<=H.size(); m+=6)
        scatter(H[m],H[m+1],H[m+2],H[m+3+c]);

        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        for(int k=0; k<p->knoz; ++k)
        {
            const int g = fid(p,i,j,k);

            if(flag[g]<=0 || sw[g]<=1.0e-12)
            continue;

            if(0.5*(a->phi(i,j,k) + a->phi(i+di,j+dj,k+dk)) < 0.0)
            continue;

            f(i,j,k) = sv[g]/sw[g];
        }

        if(c==0)
        pgc->start1(p,f,10);
        if(c==1)
        pgc->start2(p,f,11);
        if(c==2)
        pgc->start3(p,f,12);
    }

    // the velocity of the start of the step, extended into the air, for the FLIP update
    water_velocity(p,a,pgc);

    for(int n=0; n<size; ++n)
    {
        uo->V[n] = uw->V[n];
        vo->V[n] = vw->V[n];
        wo->V[n] = ww->V[n];
    }
}

// after the momentum step, before the particles move: the FLIP / PIC update of their velocity
void pic_f::g2p(lexer *p, fdm *a)
{

    for(auto &q : P)
    {
        double un,vn,wn;
        velocity(p,a,q.x,q.y,q.z,un,vn,wn);

        // in the relaxation zones the particles take the relaxed velocity of the water (PIC)
        const double c = (relaxzones && inzone(p,q.x,q.y)) ? 0.0 : p->F82;

        const double du = un - ipol(p,*uo,1,q.x,q.y,q.z);
        const double dv = p->j_dir==1 ? vn - ipol(p,*vo,2,q.x,q.y,q.z) : 0.0;
        const double dw = wn - ipol(p,*wo,3,q.x,q.y,q.z);

        q.u = c*(q.u + du) + (1.0-c)*un;
        q.v = p->j_dir==1 ? c*(q.v + dv) + (1.0-c)*vn : 0.0;
        q.w = c*(q.w + dw) + (1.0-c)*wn;
    }
}

// momentum_rk: no convection on the faces in the water (the particles carry the momentum there)
void pic_f::mask_convection(lexer *p, fdm *a, int c, const double *before, field &F)
{
    const int di = c==0, dj = c==1, dk = c==2;

    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    if(0.5*(a->phi(i,j,k) + a->phi(i+di,j+dj,k+dk)) >= 0.0)
    {
        const int g = fid(p,i,j,k);
        F.V[g] = before[g];
    }
}
