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

#include"fnpf_ice.h"
#include"lexer.h"
#include"ice_contact.h"

namespace
{
    // 3x3 solve by Cramer's rule, A x = b
    void solve3(const double (*A)[3], const double *b, double *x)
    {
        const double det = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
                         - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                         + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

        if(fabs(det)<1.0e-300)
        {
        x[0]=x[1]=x[2]=0.0;
        return;
        }

        for(int c=0; c<3; ++c)
        {
            double M[3][3];
            for(int r=0;r<3;++r)
            for(int s=0;s<3;++s)
            M[r][s] = (s==c) ? b[r] : A[r][s];

            x[c] = ( M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                   - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                   + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]) )/det;
        }
    }
}

void fnpf_ice::derivative(const fnpf_ice_floe &fl, double *D) const
{
    // dx/dt = v,  dq/dt = 0.5*(0,w)*q,  m dv/dt = F - m g e_z,  I dw/dt = M - w x I w   (world frame)
    const double *q = fl.q;
    const double *w = fl.w;

    for(int a=0;a<3;++a)
    D[a] = fl.v[a];

    D[3] = 0.5*(-w[0]*q[1] - w[1]*q[2] - w[2]*q[3]);
    D[4] = 0.5*( w[0]*q[0] + w[1]*q[3] - w[2]*q[2]);
    D[5] = 0.5*(-w[0]*q[3] + w[1]*q[0] + w[2]*q[1]);
    D[6] = 0.5*( w[0]*q[2] - w[1]*q[1] + w[2]*q[0]);

    D[7] = fl.F[0]/fl.mass;
    D[8] = fl.F[1]/fl.mass;
    D[9] = fl.F[2]/fl.mass - g;

    double I[3][3];
    for(int a=0;a<3;++a)
    for(int b=0;b<3;++b)
    {
    I[a][b]=0.0;
    for(int k1=0;k1<3;++k1)
    for(int k2=0;k2<3;++k2)
    I[a][b] += fl.R[a][k1]*fl.Ib[k1][k2]*fl.R[b][k2];
    }

    double Iw[3], rhs[3];
    for(int a=0;a<3;++a)
    Iw[a] = I[a][0]*w[0] + I[a][1]*w[1] + I[a][2]*w[2];

    rhs[0] = fl.M[0] - (w[1]*Iw[2] - w[2]*Iw[1]);
    rhs[1] = fl.M[1] - (w[2]*Iw[0] - w[0]*Iw[2]);
    rhs[2] = fl.M[2] - (w[0]*Iw[1] - w[1]*Iw[0]);

    solve3(I,rhs,&D[10]);

    // 2D flume: surge, heave, pitch only
    if(is2D)
    {
    D[8]  = 0.0;
    D[10] = 0.0;
    D[12] = 0.0;
    }
}

void fnpf_ice::rk_stage(lexer *p)
{
    // Same stage combination as the free surface (fnpf_RK3: SSP-RK3, fnpf_RK4: classical RK4),
    // called once per stage after the stage loads F, M are known. On exit the floes hold the state
    // of the next stage, after the last stage the state at t^n+1.
    const double dt = p->dt;
    double D[NY], Y[NY];

    for(size_t n=0; n<floe.size(); ++n)
    {
        fnpf_ice_floe &fl = floe[n];

        if(fl.type!=0)
        continue;

        derivative(fl,D);
        get_state(fl,Y);
        
        // stage accelerations, the last stage feeds the bending moments
        for(int a=0;a<3;++a)
        {
        fl.acc[a] = D[7+a];
        fl.alp[a] = D[10+a];
        }

        const double *Y0 = Yn[n].data();

        if(nstage==3)
        {
            if(stage==1)
            for(int a=0;a<NY;++a) Y[a] = Y0[a] + dt*D[a];

            if(stage==2)
            for(int a=0;a<NY;++a) Y[a] = 0.75*Y0[a] + 0.25*(Y[a] + dt*D[a]);

            if(stage==3)
            for(int a=0;a<NY;++a) Y[a] = (1.0/3.0)*Y0[a] + (2.0/3.0)*(Y[a] + dt*D[a]);
        }

        if(nstage==4)
        {
            if(stage==1)
            for(int a=0;a<NY;++a) {D1[n][a]=D[a]; Y[a] = Y0[a] + 0.5*dt*D[a];}

            if(stage==2)
            for(int a=0;a<NY;++a) {D2[n][a]=D[a]; Y[a] = Y0[a] + 0.5*dt*D[a];}

            if(stage==3)
            for(int a=0;a<NY;++a) {D3[n][a]=D[a]; Y[a] = Y0[a] + dt*D[a];}

            if(stage==4)
            for(int a=0;a<NY;++a) Y[a] = Y0[a] + (dt/6.0)*(D1[n][a] + 2.0*D2[n][a] + 2.0*D3[n][a] + D[a]);
        }

        set_state(fl,Y);
        quat_normalize(fl.q);
    }
}

void fnpf_ice::contact(lexer *p)
{
    // Rank 0, after the RK step. Planar non-smooth contact on (u, v, w_z) of floes and obstacles,
    // applied as a velocity jump; positions follow the corrected velocities over the step
    // (x += dt*du), as in Moreau-Jean time stepping.
    const double dt = p->dt;

    // floes that left the domain are taken out (type 2)
    for(auto &fl : floe)
    if(fl.type==0 && fl.Awet<=0.0)
    if(fl.x[0]<p->global_xmin || fl.x[0]>p->global_xmax || (!is2D && (fl.x[1]<p->global_ymin || fl.x[1]>p->global_ymax)))
    {
    fl.type=2;
    cout<<"FNPF ice: floe "<<fl.id<<" left the domain at t = "<<p->simtime+dt<<" and is removed"<<endl;
    }

    cmap.clear();
    for(size_t n=0; n<floe.size(); ++n)
    if(floe[n].type==0 || floe[n].type==1)
    cmap.push_back(int(n));
    
    const vector<int> &map = cmap;

    vector<ice_body2D> body(map.size());

    for(size_t nb=0; nb<map.size(); ++nb)
    {
        fnpf_ice_floe &fl = floe[map[nb]];
        ice_body2D &b = body[nb];

        quat_to_matrix(fl.q,fl.R);

        b.x  = fl.x[0];
        b.y  = fl.x[1];
        b.vx = fl.v[0];
        b.vy = fl.v[1];
        b.w  = fl.w[2];
        b.rb = fl.rbound;
        b.id = fl.id;

        if(fl.type==0)
        {
        b.im = 1.0/fl.mass;
        double Izz=0.0;
        for(int k1=0;k1<3;++k1)
        for(int k2=0;k2<3;++k2)
        Izz += fl.R[2][k1]*fl.Ib[k1][k2]*fl.R[2][k2];
        b.iI = (is2D || Izz<=0.0) ? 0.0 : 1.0/Izz;
        }

        // planform at the level of the centre of mass (obstacles: as given)
        const int nv = int(fl.bx.size());
        b.px.resize(nv);
        b.py.resize(nv);
        for(int q=0; q<nv; ++q)
        {
        b.px[q] = fl.x[0] + fl.R[0][0]*fl.bx[q] + fl.R[0][1]*fl.by[q];
        b.py[q] = fl.x[1] + fl.R[1][0]*fl.bx[q] + fl.R[1][1]*fl.by[q];
        }
    }

    pcontact->solve(body,dt,is2D);

    for(size_t nb=0; nb<map.size(); ++nb)
    {
        fnpf_ice_floe &fl = floe[map[nb]];
        const ice_body2D &b = body[nb];

        fl.Fc[0] = b.fx;
        fl.Fc[1] = b.fy;
        fl.ncontact = b.ncontact;

        if(fl.type!=0 || b.ncontact==0)
        continue;

        const double du = b.vx - fl.v[0];
        const double dv = b.vy - fl.v[1];
        const double dw = b.w  - fl.w[2];

        fl.v[0] = b.vx;
        fl.v[1] = b.vy;
        fl.w[2] = b.w;

        fl.x[0] += dt*du;
        fl.x[1] += dt*dv;

        // yaw correction: rotate about the vertical by dt*dw
        const double ha = 0.5*dt*dw;
        const double rq[4] = {cos(ha), 0.0, 0.0, sin(ha)};
        const double q0=fl.q[0], q1=fl.q[1], q2=fl.q[2], q3=fl.q[3];
        fl.q[0] = rq[0]*q0 - rq[3]*q3;
        fl.q[1] = rq[0]*q1 - rq[3]*q2;
        fl.q[2] = rq[0]*q2 + rq[3]*q1;
        fl.q[3] = rq[0]*q3 + rq[3]*q0;
        quat_normalize(fl.q);
    }
}
