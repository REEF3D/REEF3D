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

#include"net_membrane.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<mpi.h>
#include<sys/stat.h>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<cmath>

net_membrane::net_membrane(int num, const membrane_param &mp) : nMem(num), prm(mp), Afloor(0.0),
                           delta(0.0), Kn(0.0), Kt(0.0), Fx(0.0), Fy(0.0), Fz(0.0), Fzfloor(0.0), Qleak(0.0), urelmax(0.0), dh(0.0),
                           printtime(0.0), printcount(0), outdir("./REEF3D_NHFLOW_Membrane")
{
}

net_membrane::~net_membrane()
{
}

void net_membrane::initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->A520!=1)
    {
        // The incremental pressure correction scheme (A 520 2) accumulates the membrane pressure jump in
        // the old pressure of the predictor, which the wide collocated gradients see differently from the
        // compact Poisson operator: in tests a spurious circulation below the bag floor keeps growing.
        if(p->mpirank==0)
        cout<<"\n!!! X 330 membrane: requires the non-hydrostatic pressure projection A 520 1 !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    if(prm.shape==2 && p->j_dir==0)
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330 membrane: cylinder bag in a 2D simulation, use a box bag !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    // default layer half width: 1.5 times the largest cell size in any direction,
    // so that the layer is at least ~3 cells thick across walls (dx) and floor (sigma layers)
    double dmax=0.0, hmin=1.0e20;

    for(i=0; i<p->knox; ++i)
    {
        dmax = MAX(dmax, p->DXN[IP]);
        hmin = MIN(hmin, p->DXN[IP]);
    }

    if(p->j_dir==1)
    for(j=0; j<p->knoy; ++j)
    {
        dmax = MAX(dmax, p->DYN[JP]);
        hmin = MIN(hmin, p->DYN[JP]);
    }

    LOOP
    if(p->wet[IJ]==1)
    dmax = MAX(dmax, p->DZN[KP]*d->WL(i,j));

    dmax = pgc->globalmax(dmax);
    hmin = pgc->globalmin(hmin);

    delta = prm.delta>0.0 ? prm.delta : 1.5*dmax;

    if(prm.h<=0.0)
    prm.h = hmin;

    // the integral of the indicator across the layer is 1.5 delta (plateau delta, two tapers delta/4)
    Kn = prm.Rn/(1.5*delta);
    Kt = prm.Rt/(1.5*delta);

    mesh(p);

    // cell map storage
    slot_.assign(p->imax*p->jmax*(p->kmax+2),-1);

    xc_.resize(p->knox);
    for(i=0; i<p->knox; ++i)
    xc_[i] = p->XP[IP];

    yc_.resize(p->knoy);
    for(j=0; j<p->knoy; ++j)
    yc_[j] = p->YP[JP];

    tf_.assign(3*tri_.size(),0.0);

    if(p->mpirank==0)
    {
        mkdir(outdir.c_str(),0777);

        cout<<"Membrane "<<nMem<<" ("<<prm.name<<"): "<<(prm.shape==1?"box":"cylinder")<<", "<<tri_.size()<<" triangles, "
            <<"delta = "<<delta<<" m, R_n = "<<prm.Rn<<" m/s (K_n = "<<Kn<<" 1/s), R_t = "<<prm.Rt<<" m/s, "
            <<"floor area = "<<Afloor<<" m^2, fill = "<<prm.fill<<" m"<<(prm.poisson==1?"":", Poisson mobility OFF")<<endl;

        if(delta < 1.5*dmax)
        cout<<"Membrane "<<nMem<<": delta < 1.5 max(dx,dy,dz) - the smeared layer may leave gaps between cells"<<endl;

        ofstream ts((outdir+"/REEF3D_NHFLOW_Membrane_"+to_string(nMem)+".dat").c_str());
        ts<<"# membrane "<<nMem<<" "<<prm.name<<"  delta "<<delta<<"  Rn "<<prm.Rn<<"  Rt "<<prm.Rt<<"  Afloor "<<Afloor<<"\n";
        ts<<"# time  eta_in  eta_out  dh  Q_leak[m3/s]  Fx  Fy  Fz  Fz_floor  Fz_floor_hydrostatic(-rho g dh A)  max|u_n,rel|_layer  max|U|  water_volume\n";
        ts.close();
    }
}

// ---------------------------------------------------------------------------------------------
// geometry
// ---------------------------------------------------------------------------------------------

void net_membrane::mesh(lexer *p)
{
    x_.clear();
    tri_.clear();
    tn_.clear();
    tc_.clear();
    ta_.clear();
    ttag_.clear();

    const double h = prm.h;
    const double lz = prm.zt - prm.zb;
    const int nz = MAX(1,(int)ceil(lz/h));

    if(prm.shape==1)
    {
        // in 2D the y-extent is irrelevant: the panels are extended far beyond the single cell row
        const double ylo = p->j_dir==1 ? prm.y0 : -1.0e3;
        const double yhi = p->j_dir==1 ? prm.y1 : 1.0e3;
        const double lx = prm.x1 - prm.x0;
        const double ly = yhi - ylo;
        const int nx = MAX(1,(int)ceil(lx/h));
        const int ny = p->j_dir==1 ? MAX(1,(int)ceil(ly/h)) : 1;

        using V3 = Eigen::Vector3d;

        // walls x = x0, x = x1
        add_panel(V3(prm.x0,ylo,prm.zb), V3(0.0,ly,0.0), V3(0.0,0.0,lz), ny, nz, V3(-1.0,0.0,0.0), 0);
        add_panel(V3(prm.x1,ylo,prm.zb), V3(0.0,ly,0.0), V3(0.0,0.0,lz), ny, nz, V3( 1.0,0.0,0.0), 0);

        // walls y = y0, y = y1
        if(p->j_dir==1)
        {
        add_panel(V3(prm.x0,prm.y0,prm.zb), V3(lx,0.0,0.0), V3(0.0,0.0,lz), nx, nz, V3(0.0,-1.0,0.0), 0);
        add_panel(V3(prm.x0,prm.y1,prm.zb), V3(lx,0.0,0.0), V3(0.0,0.0,lz), nx, nz, V3(0.0, 1.0,0.0), 0);
        }

        // floor
        add_panel(V3(prm.x0,ylo,prm.zb), V3(lx,0.0,0.0), V3(0.0,ly,0.0), nx, ny, V3(0.0,0.0,-1.0), 1);
    }

    if(prm.shape==2)
    {
        const int nt = MAX(16,(int)ceil(2.0*PI*prm.R/h));
        const int nr = MAX(1,(int)ceil(prm.R/h));

        add_cylinder_wall(nt,nz);
        add_disk(nt,nr);
    }

    xdot_.assign(x_.size(),Eigen::Vector3d::Zero());

    Afloor=0.0;
    for(size_t t=0; t<tri_.size(); ++t)
    if(ttag_[t]==1)
    Afloor += ta_[t];

    // 2D: loads act on the width of the single cell row
    if(p->j_dir==0)
    {
    j=0;
    Afloor = (prm.x1 - prm.x0)*p->DYN[JP];
    }
}

void net_membrane::add_panel(const Eigen::Vector3d &o, const Eigen::Vector3d &a, const Eigen::Vector3d &b, int na, int nb,
                             const Eigen::Vector3d &nout, int tag)
{
    const int n0 = x_.size();

    for(int ib=0; ib<=nb; ++ib)
    for(int ia=0; ia<=na; ++ia)
    x_.push_back(o + a*(double(ia)/double(na)) + b*(double(ib)/double(nb)));

    for(int ib=0; ib<nb; ++ib)
    for(int ia=0; ia<na; ++ia)
    {
        const int q0 = n0 + ib*(na+1) + ia;
        const int q1 = q0 + 1;
        const int q2 = q0 + (na+1);
        const int q3 = q2 + 1;

        add_tri(q0,q1,q3,nout,tag);
        add_tri(q0,q3,q2,nout,tag);
    }
}

void net_membrane::add_cylinder_wall(int nt, int nz)
{
    const int n0 = x_.size();

    for(int iz=0; iz<=nz; ++iz)
    for(int it=0; it<nt; ++it)
    {
        const double th = 2.0*PI*double(it)/double(nt);
        x_.push_back(Eigen::Vector3d(prm.xc + prm.R*cos(th), prm.yc + prm.R*sin(th), prm.zb + (prm.zt-prm.zb)*double(iz)/double(nz)));
    }

    for(int iz=0; iz<nz; ++iz)
    for(int it=0; it<nt; ++it)
    {
        const int q0 = n0 + iz*nt + it;
        const int q1 = n0 + iz*nt + (it+1)%nt;
        const int q2 = q0 + nt;
        const int q3 = q1 + nt;

        const double th = 2.0*PI*(double(it)+0.5)/double(nt);
        const Eigen::Vector3d nout(cos(th),sin(th),0.0);

        add_tri(q0,q1,q3,nout,0);
        add_tri(q0,q3,q2,nout,0);
    }
}

void net_membrane::add_disk(int nt, int nr)
{
    const Eigen::Vector3d nout(0.0,0.0,-1.0);
    const int nc = x_.size();

    x_.push_back(Eigen::Vector3d(prm.xc,prm.yc,prm.zb));

    const int n0 = x_.size();

    for(int ir=1; ir<=nr; ++ir)
    for(int it=0; it<nt; ++it)
    {
        const double th = 2.0*PI*double(it)/double(nt);
        const double r = prm.R*double(ir)/double(nr);
        x_.push_back(Eigen::Vector3d(prm.xc + r*cos(th), prm.yc + r*sin(th), prm.zb));
    }

    // centre fan
    for(int it=0; it<nt; ++it)
    add_tri(nc, n0 + it, n0 + (it+1)%nt, nout, 1);

    // rings
    for(int ir=1; ir<nr; ++ir)
    for(int it=0; it<nt; ++it)
    {
        const int q0 = n0 + (ir-1)*nt + it;
        const int q1 = n0 + (ir-1)*nt + (it+1)%nt;
        const int q2 = q0 + nt;
        const int q3 = q1 + nt;

        add_tri(q0,q2,q3,nout,1);
        add_tri(q0,q3,q1,nout,1);
    }
}

void net_membrane::add_tri(int a, int b, int c, const Eigen::Vector3d &nout, int tag)
{
    Eigen::Vector3d n = (x_[b]-x_[a]).cross(x_[c]-x_[a]);
    const double A2 = n.norm();

    if(A2<1.0e-20)
    return;

    n/=A2;

    // orient outward (towards the outside of the bag)
    if(n.dot(nout)<0.0)
    {
        swap(b,c);
        n = -n;
    }

    tri_.push_back({a,b,c});
    tn_.push_back(n);
    tc_.push_back((x_[a]+x_[b]+x_[c])/3.0);
    ta_.push_back(0.5*A2);
    ttag_.push_back(tag);
}

bool net_membrane::inside_footprint(double xp, double yp, double margin) const
{
    if(prm.shape==1)
    {
        if(xp<=prm.x0+margin || xp>=prm.x1-margin)
        return false;

        if(yc_.size()>1 && (yp<=prm.y0+margin || yp>=prm.y1-margin))
        return false;

        return true;
    }

    const double r = sqrt((xp-prm.xc)*(xp-prm.xc) + (yp-prm.yc)*(yp-prm.yc));

    return r < prm.R - margin;
}

bool net_membrane::outside_footprint(double xp, double yp, double margin) const
{
    if(prm.shape==1)
    {
        if(xp<prm.x0-margin || xp>prm.x1+margin)
        return true;

        if(yc_.size()>1 && (yp<prm.y0-margin || yp>prm.y1+margin))
        return true;

        return false;
    }

    const double r = sqrt((xp-prm.xc)*(xp-prm.xc) + (yp-prm.yc)*(yp-prm.yc));

    return r > prm.R + margin;
}

void net_membrane::fill_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(fabs(prm.fill)<1.0e-20)
    return;

    SLICELOOP4
    if(inside_footprint(p->XP[IP],p->YP[JP],0.0))
    {
        d->eta(i,j) += prm.fill;
        d->eta_n(i,j) = d->eta(i,j);
        d->WL(i,j) = d->eta(i,j) + d->depth(i,j);
    }

    pgc->gcsl_start4(p,d->eta,50);
    pgc->gcsl_start4(p,d->eta_n,50);
    pgc->gcsl_start4(p,d->WL,50);
}

// ---------------------------------------------------------------------------------------------
// cell map: cells within delta of the triangulation
// ---------------------------------------------------------------------------------------------

Eigen::Vector3d net_membrane::closest_point(const Eigen::Vector3d &P, const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                                            const Eigen::Vector3d &c, double &u, double &v, double &w)
{
    // closest point on triangle abc to P, barycentric weights (u,v,w) for (a,b,c); Ericson, RTCD 5.1.5
    const Eigen::Vector3d ab = b-a, ac = c-a, ap = P-a;
    const double d1 = ab.dot(ap), d2 = ac.dot(ap);
    if(d1<=0.0 && d2<=0.0) {u=1.0; v=0.0; w=0.0; return a;}

    const Eigen::Vector3d bp = P-b;
    const double d3 = ab.dot(bp), d4 = ac.dot(bp);
    if(d3>=0.0 && d4<=d3) {u=0.0; v=1.0; w=0.0; return b;}

    const double vc = d1*d4 - d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0)
    {
        const double t = d1/(d1-d3);
        u=1.0-t; v=t; w=0.0;
        return a + t*ab;
    }

    const Eigen::Vector3d cp = P-c;
    const double d5 = ab.dot(cp), d6 = ac.dot(cp);
    if(d6>=0.0 && d5<=d6) {u=0.0; v=0.0; w=1.0; return c;}

    const double vb = d5*d2 - d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0)
    {
        const double t = d2/(d2-d6);
        u=1.0-t; v=0.0; w=t;
        return a + t*ac;
    }

    const double va = d3*d6 - d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0)
    {
        const double t = (d4-d3)/((d4-d3)+(d5-d6));
        u=0.0; v=1.0-t; w=t;
        return b + t*(c-b);
    }

    const double denom = 1.0/(va+vb+vc);
    v = vb*denom;
    w = vc*denom;
    u = 1.0-v-w;
    return a + ab*v + ac*w;
}

double net_membrane::indicator(double dist) const
{
    // plateau of full resistance for dist < delta/2, cosine taper to zero at dist = delta.
    // With a uniform resistance across most of the layer the pressure (and free surface) drop over
    // the membrane is a linear ramp over the layer, which the wide collocated gradients of the
    // pressure correction and the compact hydrostatic gradient see alike. A peaked indicator
    // concentrates the jump in one cell, where the two gradients disagree and drive a spurious
    // circulation under the bag floor.
    const double h = 0.5*delta;
    
    if(dist<=h)
    return 1.0;
    
    if(dist>=delta)
    return 0.0;
    
    return 0.5*(1.0 + cos(PI*(dist-h)/h));
}

Eigen::Vector3d net_membrane::membrane_vel(int t, double w0, double w1, double w2) const
{
    return w0*xdot_[tri_[t][0]] + w1*xdot_[tri_[t][1]] + w2*xdot_[tri_[t][2]];
}

void net_membrane::build_map(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // The sigma grid moves with the free surface, so the map is rebuilt in every stage.
    // Cost ~ (triangles) x (cells in the delta-expanded bounding box of a triangle).
    cells_.clear();

    const double xs = p->XN[0+marge]-delta;
    const double xe = p->XN[p->knox+marge]+delta;
    const double ys = p->YN[0+marge]-delta;
    const double ye = p->YN[p->knoy+marge]+delta;

    // normals closer than this are treated as the same surface
    const double samenormal = 0.9;

    double u,v,w;

    for(size_t t=0; t<tri_.size(); ++t)
    {
        const Eigen::Vector3d &a = x_[tri_[t][0]];
        const Eigen::Vector3d &b = x_[tri_[t][1]];
        const Eigen::Vector3d &c = x_[tri_[t][2]];

        const double bxmin = min(a(0),min(b(0),c(0)))-delta;
        const double bxmax = max(a(0),max(b(0),c(0)))+delta;
        const double bymin = min(a(1),min(b(1),c(1)))-delta;
        const double bymax = max(a(1),max(b(1),c(1)))+delta;
        const double bzmin = min(a(2),min(b(2),c(2)))-delta;
        const double bzmax = max(a(2),max(b(2),c(2)))+delta;

        if(bxmax<xs || bxmin>xe)
        continue;

        if(p->j_dir==1 && (bymax<ys || bymin>ye))
        continue;

        const int ia = lower_bound(xc_.begin(),xc_.end(),bxmin) - xc_.begin();
        const int ib = int(upper_bound(xc_.begin(),xc_.end(),bxmax) - xc_.begin()) - 1;

        int ja=0, jb=0;
        if(p->j_dir==1)
        {
        ja = lower_bound(yc_.begin(),yc_.end(),bymin) - yc_.begin();
        jb = int(upper_bound(yc_.begin(),yc_.end(),bymax) - yc_.begin()) - 1;
        }

        for(i=ia; i<=ib; ++i)
        for(j=ja; j<=jb; ++j)
        {
            if(p->wet[IJ]==0)
            continue;

            for(k=0; k<p->knoz; ++k)
            {
                if(p->flag4[IJK]<=0)
                continue;

                const double zp = p->ZSP[IJK];

                if(zp<bzmin || zp>bzmax)
                continue;

                const Eigen::Vector3d P(p->XP[IP], p->YP[JP], zp);
                const Eigen::Vector3d C = closest_point(P,a,b,c,u,v,w);
                const double dist = (P-C).norm();

                if(dist>=delta)
                continue;

                int &s = slot_[IJK];

                if(s<0)
                {
                    s = cells_.size();
                    cellentry e;
                    e.i=i; e.j=j; e.k=k;
                    e.t1=t; e.d1=dist; e.w0=u; e.w1=v; e.w2=w;
                    e.t2=-1; e.d2=1.0e20;
                    e.H1=e.H2=0.0;
                    cells_.push_back(e);
                    continue;
                }

                cellentry &e = cells_[s];

                const bool differs = fabs(tn_[t].dot(tn_[e.t1])) < samenormal;

                if(dist<e.d1)
                {
                    // the old closest triangle becomes the second surface if its normal differs
                    if(differs)
                    {
                        e.t2 = e.t1;
                        e.d2 = e.d1;
                    }

                    e.t1=t; e.d1=dist; e.w0=u; e.w1=v; e.w2=w;
                }
                else if(differs && dist<e.d2)
                {
                    e.t2=t;
                    e.d2=dist;
                }
            }
        }
    }

    for(auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;
        slot_[IJK] = -1;

        e.H1 = indicator(e.d1);
        e.H2 = e.t2>=0 ? indicator(e.d2) : 0.0;
    }
}

// ---------------------------------------------------------------------------------------------
// coupling
// ---------------------------------------------------------------------------------------------

void net_membrane::mobility_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha)
{
    build_map(p,d,pgc);
    
    if(prm.poisson==0)
    return;
    
    const double a = alpha*p->dt;
    
    // mobility in the pressure Poisson equation and velocity correction
    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;
        d->MBETA[IJK] = MIN(d->MBETA[IJK], 1.0/(1.0 + a*Kn*MAX(e.H1,e.H2)));
    }
}

void net_membrane::forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    const double a = alpha*p->dt;
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        const Eigen::Vector3d &n1 = tn_[e.t1];

        Eigen::Matrix3d A = Kn*e.H1*(n1*n1.transpose()) + Kt*e.H1*(I - n1*n1.transpose());

        if(e.t2>=0)
        {
            const Eigen::Vector3d &n2 = tn_[e.t2];
            A += Kn*e.H2*(n2*n2.transpose());
        }

        const Eigen::Vector3d um = membrane_vel(e.t1,e.w0,e.w1,e.w2);

        Eigen::Vector3d rel(d->U[IJK]-um(0), d->V[IJK]-um(1), d->W[IJK]-um(2));

        if(p->j_dir==0)
        rel(1)=0.0;

        // implicit resistance over the stage: (I + a A)(u - u_m)^new = (u - u_m)^*
        const Eigen::Vector3d relnew = (I + a*A).ldlt().solve(rel);
        const Eigen::Vector3d du = relnew - rel;

        d->U[IJK] += du(0);
        UH[IJK]   += du(0)*WL(i,j);

        if(p->j_dir==1)
        {
        d->V[IJK] += du(1);
        VH[IJK]   += du(1)*WL(i,j);
        }

        d->W[IJK] += du(2);
        WH[IJK]   += du(2)*WL(i,j);
    }
}

void net_membrane::reaction_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, slice &WL, bool finalize)
{
    // load on the membrane from the final (projected) velocity: F = rho H A (u^{n+1} - u_m) dV
    fill(tf_.begin(),tf_.end(),0.0);

    double Q=0.0, urm=0.0;
    const double rho = p->W1;

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        if(p->wet[IJ]==0)
        continue;

        const Eigen::Vector3d um = membrane_vel(e.t1,e.w0,e.w1,e.w2);

        Eigen::Vector3d rel(d->U[IJK]-um(0), d->V[IJK]-um(1), d->W[IJK]-um(2));

        if(p->j_dir==0)
        rel(1)=0.0;

        const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);

        const Eigen::Vector3d &n1 = tn_[e.t1];
        const double un1 = n1.dot(rel);

        const Eigen::Vector3d f1 = rho*dV*e.H1*(Kn*un1*n1 + Kt*(rel - un1*n1));

        tf_[3*e.t1+0] += f1(0);
        tf_[3*e.t1+1] += f1(1);
        tf_[3*e.t1+2] += f1(2);

        double un2=0.0;

        if(e.t2>=0)
        {
            const Eigen::Vector3d &n2 = tn_[e.t2];
            un2 = n2.dot(rel);

            const Eigen::Vector3d f2 = rho*dV*e.H2*Kn*un2*n2;

            tf_[3*e.t2+0] += f2(0);
            tf_[3*e.t2+1] += f2(1);
            tf_[3*e.t2+2] += f2(2);
        }

        // flux through the membrane: int u_n dA = (1/(1.5 delta)) int H u_n dV, outward positive
        Q += (e.H1*un1 + e.H2*un2)*dV/(1.5*delta);

        if(e.H1>0.5)
        urm = MAX(urm,fabs(un1));
    }

    MPI_Allreduce(MPI_IN_PLACE,tf_.data(),(int)tf_.size(),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    Qleak = pgc->globalsum(Q);
    urelmax = pgc->globalmax(urm);

    Fx=Fy=Fz=Fzfloor=0.0;

    for(size_t t=0; t<tri_.size(); ++t)
    {
        Fx += tf_[3*t+0];
        Fy += tf_[3*t+1];
        Fz += tf_[3*t+2];

        if(ttag_[t]==1)
        Fzfloor += tf_[3*t+2];
    }

    if(finalize)
    {
        print_timeseries(p,d,pgc);

        if(prm.printdt>0.0 && (p->simtime>=printtime || p->count==0))
        {
            print_vtp(p);
            printtime += prm.printdt;
        }
    }
}

void net_membrane::netForces(lexer *p, double &Xne, double &Yne, double &Zne, double &Kne, double &Mne, double &Nne)
{
    // total hydrodynamic load on the membrane and its moment about the origin
    Xne=Yne=Zne=Kne=Mne=Nne=0.0;

    for(size_t t=0; t<tri_.size(); ++t)
    {
        const Eigen::Vector3d F(tf_[3*t+0],tf_[3*t+1],tf_[3*t+2]);
        const Eigen::Vector3d M = tc_[t].cross(F);

        Xne += F(0);
        Yne += F(1);
        Zne += F(2);
        Kne += M(0);
        Mne += M(1);
        Nne += M(2);
    }
}

// ---------------------------------------------------------------------------------------------
// static overpressure of the bag below the floor
// ---------------------------------------------------------------------------------------------

double net_membrane::smoothstep(double s) const
{
    // integral of the layer indicator, normalised to go from 0 (s <= -delta) to 1 (s >= delta):
    // the pressure drop of a uniform Darcy flux through the layer has exactly this shape
    const double h = 0.5*delta;
    const double u = fabs(s);
    double G;
    
    if(u<=h)
    G = u;
    else if(u<delta)
    G = h + 0.5*(u-h) + h/(2.0*PI)*sin(PI*(u-h)/h);
    else
    G = 0.75*delta;
    
    return 0.5 + (s>=0.0 ? 1.0 : -1.0)*G/(1.5*delta);
}

double net_membrane::footprint_distance(double xp, double yp, double &gx, double &gy) const
{
    // signed horizontal distance to the footprint boundary, positive inside, and its gradient
    gx=gy=0.0;
    
    if(prm.shape==1)
    {
        double s = xp-prm.x0;   gx = 1.0;
        
        if(prm.x1-xp<s) {s = prm.x1-xp; gx=-1.0;}
        
        if(yc_.size()>1)
        {
            if(yp-prm.y0<s) {s = yp-prm.y0; gx=0.0; gy= 1.0;}
            if(prm.y1-yp<s) {s = prm.y1-yp; gx=0.0; gy=-1.0;}
        }
        
        return s;
    }
    
    const double rx = xp-prm.xc, ry = yp-prm.yc;
    const double r = sqrt(rx*rx + ry*ry);
    
    if(r>1.0e-12)
    {
        gx = -rx/r;
        gy = -ry/r;
    }
    
    return prm.R - r;
}

void net_membrane::static_pressure_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    // NHFLOW carries one free surface per column, so a column below the bag floor sees the inner
    // water level eta_in in its hydrostatic pressure, while the water there is connected to the
    // outside, eta_out. The non-hydrostatic pressure would have to jump by -rho g dh across the floor,
    // and that sharp jump on the sloping sigma levels at the bag edge (where eta ramps from eta_in to
    // eta_out) leaves pressure gradient errors that drive a spurious circulation under the floor.
    // The known static part is therefore applied as a prescribed pressure
    //
    //      p_m = -rho g dh A(s_xy) B(z_b - z),     dh = eta_in - eta_out (footprint means)
    //
    // with A, B the integral of the layer indicator across the walls and the floor. Its gradient is
    // evaluated analytically in physical coordinates and added as a body force in the stage, so the
    // solved non-hydrostatic pressure stays smooth.
    double ein=0.0, ain=0.0, eout=0.0, aout=0.0, gx, gy;
    
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        const double A = p->DXN[IP]*p->DYN[JP];
        const double s = footprint_distance(p->XP[IP],p->YP[JP],gx,gy);
        
        if(s>delta)
        {
            ein += d->eta(i,j)*A;
            ain += A;
        }
        
        if(s<-delta && s>-3.0*delta)
        {
            eout += d->eta(i,j)*A;
            aout += A;
        }
    }
    
    ein  = pgc->globalsum(ein);
    ain  = pgc->globalsum(ain);
    eout = pgc->globalsum(eout);
    aout = pgc->globalsum(aout);
    
    if(ain<=0.0 || aout<=0.0)
    return;
    
    dh = ein/ain - eout/aout;
    
    const double a = alpha*p->dt;
    const double g = fabs(p->W22);
    
    LOOP
    WETDRYDEEP
    {
        const double s = footprint_distance(p->XP[IP],p->YP[JP],gx,gy);
        const double sz = prm.zb - p->ZSP[IJK];
        
        if(s<=-delta || sz<=-delta)
        continue;
        
        const double A  = smoothstep(s);
        const double B  = smoothstep(sz);
        const double dA = s<delta  ? indicator(fabs(s))/(1.5*delta)  : 0.0;
        const double dB = sz<delta ? indicator(fabs(sz))/(1.5*delta) : 0.0;
        
        d->MCHI[IJK] = MAX(d->MCHI[IJK], A*B);

        // f = -grad(p_m)/rho
        const double fx =  g*dh*dA*B*gx;
        const double fy =  g*dh*dA*B*gy;
        const double fz = -g*dh*A*dB;
        
        d->U[IJK] += a*fx;
        UH[IJK]   += a*fx*WL(i,j);
        
        if(p->j_dir==1)
        {
        d->V[IJK] += a*fy;
        VH[IJK]   += a*fy*WL(i,j);
        }
        
        d->W[IJK] += a*fz;
        WH[IJK]   += a*fz*WL(i,j);
    }
}

// ---------------------------------------------------------------------------------------------
// output
// ---------------------------------------------------------------------------------------------

void net_membrane::print_timeseries(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double ein=0.0, ain=0.0, eout=0.0, aout=0.0, umax=0.0, vol=0.0;

    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        const double A = p->DXN[IP]*p->DYN[JP];

        vol += d->WL(i,j)*A;

        if(inside_footprint(p->XP[IP],p->YP[JP],delta))
        {
            ein += d->eta(i,j)*A;
            ain += A;
        }

        if(outside_footprint(p->XP[IP],p->YP[JP],delta))
        {
            eout += d->eta(i,j)*A;
            aout += A;
        }
    }

    LOOP
    if(p->wet[IJ]==1)
    umax = MAX(umax, sqrt(d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK] + d->W[IJK]*d->W[IJK]));

    ein  = pgc->globalsum(ein);
    ain  = pgc->globalsum(ain);
    eout = pgc->globalsum(eout);
    aout = pgc->globalsum(aout);
    umax = pgc->globalmax(umax);
    vol  = pgc->globalsum(vol);

    if(p->mpirank==0)
    {
        const double etain  = ain>0.0  ? ein/ain  : 0.0;
        const double etaout = aout>0.0 ? eout/aout : 0.0;
        const double dhl = etain - etaout;

        ofstream ts((outdir+"/REEF3D_NHFLOW_Membrane_"+to_string(nMem)+".dat").c_str(), ios::app);
        ts<<setprecision(10)<<p->simtime<<" "<<etain<<" "<<etaout<<" "<<dhl<<" "<<Qleak<<" "
          <<Fx<<" "<<Fy<<" "<<Fz<<" "<<Fzfloor<<" "<<-p->W1*fabs(p->W22)*dhl*Afloor<<" "<<urelmax<<" "<<umax<<" "<<vol<<"\n";
    }
}

void net_membrane::print_vtp(lexer *p)
{
    if(p->mpirank!=0)
    return;

    const string name = outdir+"/REEF3D_NHFLOW_Membrane_"+to_string(nMem)+"_"+to_string(printcount)+".vtp";
    ++printcount;

    ofstream out(name.c_str());

    out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n";
    out<<"<Piece NumberOfPoints=\""<<x_.size()<<"\" NumberOfPolys=\""<<tri_.size()<<"\">\n";
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"<<p->simtime<<"</DataArray></FieldData>\n";

    out<<"<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(const auto &x : x_)
    out<<x(0)<<" "<<x(1)<<" "<<x(2)<<"\n";
    out<<"</DataArray>\n</Points>\n";

    out<<"<CellData Scalars=\"dp\">\n";

    // pressure jump across the panel, positive when the inside pressure is higher
    out<<"<DataArray type=\"Float64\" Name=\"dp\" format=\"ascii\">\n";
    for(size_t t=0; t<tri_.size(); ++t)
    {
        const Eigen::Vector3d F(tf_[3*t+0],tf_[3*t+1],tf_[3*t+2]);
        out<<F.dot(tn_[t])/ta_[t]<<"\n";
    }
    out<<"</DataArray>\n";

    out<<"<DataArray type=\"Float64\" Name=\"force\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(size_t t=0; t<tri_.size(); ++t)
    out<<tf_[3*t+0]<<" "<<tf_[3*t+1]<<" "<<tf_[3*t+2]<<"\n";
    out<<"</DataArray>\n";

    out<<"<DataArray type=\"Int32\" Name=\"tag\" format=\"ascii\">\n";
    for(size_t t=0; t<tri_.size(); ++t)
    out<<ttag_[t]<<"\n";
    out<<"</DataArray>\n</CellData>\n";

    out<<"<Polys>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(const auto &tr : tri_)
    out<<tr[0]<<" "<<tr[1]<<" "<<tr[2]<<"\n";
    out<<"</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(size_t t=0; t<tri_.size(); ++t)
    out<<3*(t+1)<<"\n";
    out<<"</DataArray>\n</Polys>\n</Piece>\n</PolyData>\n</VTKFile>\n";
}
