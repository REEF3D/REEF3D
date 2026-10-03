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
#include"iqn_ils.h"
#include<mpi.h>
#include<sys/stat.h>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<cmath>

net_membrane::net_membrane(int num, const membrane_param &mp) : nMem(num), prm(mp), Afloor(0.0),
                           delta(0.0), Kn(0.0), Kt(0.0), Fx(0.0), Fy(0.0), Fz(0.0), Fzfloor(0.0), Qleak(0.0), urelmax(0.0), dh(0.0), etaref(0.0),
                           printtime(0.0), printcount(0), outdir("./REEF3D_NHFLOW_Membrane"), vtpdir("./REEF3D_NHFLOW_Membrane_VTP")
{
}

net_membrane::~net_membrane()
{
    delete pqn_;
    free_factor();
}

void net_membrane::initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->A520!=1 && p->A520!=2)
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330 membrane: requires the non-hydrostatic pressure scheme A 520 1 or A 520 2 !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    if(prm.structure>0 && p->X10==0 && prm.structure==1)
    {
        if(p->mpirank==0)
        cout<<"Membrane "<<nMem<<": structure rigid without floating body (X 10), the membrane stays in place"<<endl;
    }
    
    // default coupling of a flexible membrane: iterated. The staggered scheme gives the fabric an extra inertia
    // ~ rho R_n dt per area, so its dynamics depend on the time step; it stays available as 'coupling staggered'
    if(prm.coupling<0)
    prm.coupling = prm.structure==2 ? 1 : 0;
    
    if(prm.coupling==1 && prm.structure!=2 && p->mpirank==0)
    cout<<"Membrane "<<nMem<<": coupling iterated applies to flexible membranes only, ignored"<<endl;
    
    if(prm.structure==2 && p->X10==1 && prm.coupling==0 && p->mpirank==0)
    cout<<"Membrane "<<nMem<<": WARNING flexible membrane on a freely floating body (X 10 1): the collar coupling is "
        <<"implicit, but the flexible membrane follows the fluid only with an extra inertia ~ rho R_n dt per area "
        <<"(implicit porous coupling); in waves the collar is held by the bag and surge was unstable in the tests. "
        <<"Use coupling iterated, a prescribed collar motion (X 10 2) or structure rigid"<<endl;
    
    // tangential resistance: a moving membrane carries the fluid of its layer along (the Poisson mobility beta is
    // isotropic, so the layer can not be moved tangentially by pressure); a fixed one only blocks the normal flow
    if(prm.Rt<0.0)
    prm.Rt = moving() ? prm.Rn : 0.0;
    
    // floorpressure 3 matches the discrete free-surface ramp in grid coordinates, a moving membrane uses 0
    if(prm.floorp<0)
    prm.floorp = moving() ? 0 : 3;
    
    if(moving() && prm.floorp!=0)
    {
        if(p->mpirank==0)
        cout<<"Membrane "<<nMem<<": floorpressure "<<prm.floorp<<" needs a fixed membrane, using floorpressure 0"<<endl;
        
        prm.floorp=0;
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

    // triangle size: fixed membrane the smallest cell size (geometry), moving membrane the largest, so that every
    // node of the structure receives loads from the cells of the layer
    // strong coupling: elements as large as the coupling layer (1.5 delta); shorter structural modes are not seen by the
    // smeared layer, are nearly free and make the coupled fixed point slow to converge (3D: 7 instead of > 100 iterations)
    if(prm.h<=0.0)
    prm.h = iterated() ? 1.5*delta : (moving() ? dmax : hmin);

    // the integral of the indicator across the layer is 1.5 delta (plateau delta, two tapers delta/4)
    Kn = prm.Rn/(1.5*delta);
    Kt = prm.Rt/(1.5*delta);

    mesh(p);
    
    tn0_ = tn_;
    
    ini_structure(p,pgc);

    // cell map storage
    slot_.assign(p->imax*p->jmax*(p->kmax+2),-1);

    xc_.resize(p->knox);
    for(i=0; i<p->knox; ++i)
    xc_[i] = p->XP[IP];

    yc_.resize(p->knoy);
    for(j=0; j<p->knoy; ++j)
    yc_[j] = p->YP[JP];

    tf_.assign(3*tri_.size(),0.0);
    nf_.assign(3*x_.size(),0.0);
    zf_.assign(p->imax*p->jmax,prm.zb);
    
    floor_geometry(p);

    if(p->mpirank==0)
    {
        mkdir(outdir.c_str(),0777);
        
        if(prm.printdt!=0.0)
        mkdir(vtpdir.c_str(),0777);

        cout<<"Membrane "<<nMem<<" ("<<prm.name<<"): "<<(prm.shape==1?"box":"cylinder")<<", "<<tri_.size()<<" triangles, "
            <<"delta = "<<delta<<" m, R_n = "<<prm.Rn<<" m/s (K_n = "<<Kn<<" 1/s), R_t = "<<prm.Rt<<" m/s, "
            <<"floor area = "<<Afloor<<" m^2, fill = "<<prm.fill<<" m, floorpressure "<<prm.floorp
            <<", projections "<<prm.projections<<(prm.poisson==1?"":", Poisson mobility OFF")<<endl;
        
        const char *sname[3] = {"fixed","rigid","flexible"};
        cout<<"Membrane "<<nMem<<": structure "<<sname[prm.structure]<<", "<<x_.size()<<" nodes";
        
        if(prm.structure==2)
        {
            int na=0;
            for(size_t q=0; q<att_.size(); ++q)
            na += att_[q];
            
            cout<<", "<<edge_.size()<<" edges, "<<na<<" attached nodes, mass "<<prm.mA<<" kg/m^2, density "<<prm.rhom
                <<" kg/m^3, E t = "<<prm.EA<<" N/m, damping "<<prm.zeta<<", sinker "<<prm.sinker<<" N/m";
        }
        cout<<endl;

        if(delta < 1.5*dmax)
        cout<<"Membrane "<<nMem<<": delta < 1.5 max(dx,dy,dz) - the smeared layer may leave gaps between cells"<<endl;

        ofstream ts((outdir+"/REEF3D_NHFLOW_Membrane_"+to_string(nMem)+".dat").c_str());
        ts<<"# membrane "<<nMem<<" "<<prm.name<<"  delta "<<delta<<"  Rn "<<prm.Rn<<"  Rt "<<prm.Rt<<"  Afloor "<<Afloor<<"\n";
        ts<<"# time  eta_in  eta_out  dh  Q_leak[m3/s]  Fx  Fy  Fz  Fz_floor  Fz_floor_hydrostatic(-rho g dh A)  max|u_n,rel|_layer  max|U|  water_volume"
            <<"  Fx_body  Fy_body  Fz_body  floor_z_mean  floor_z_min  max|u_node|  max_tension[N/m]"
            <<(iterated() ? "  coupling_iterations  coupling_residual" : "");
        
        if(collar())
        {
            ts<<"  collar_x  collar_y  collar_z  collar_zmin  collar_zmax  collar_Mbend_max[Nm]  collar_Naxial_max[N]";
            for(size_t m=0; m<prm.moor.size(); ++m)
            ts<<"  T_mooring"<<m;
        }
        ts<<"\n";
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
        // in 2D the panels span the single cell row, so masses and loads refer to the same width
        j=0;
        const double ylo = p->j_dir==1 ? prm.y0 : p->YN[JP];
        const double yhi = p->j_dir==1 ? prm.y1 : p->YN[JP1];
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

    merge_nodes();
    
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

void net_membrane::merge_nodes()
{
    // the panels are meshed separately: merge coincident nodes along the edges, so the structure is connected
    double lmin=1.0e20;
    for(const auto &t : tri_)
    for(int q=0; q<3; ++q)
    lmin = MIN(lmin, (x_[t[q]]-x_[t[(q+1)%3]]).norm());
    
    const double tol = 1.0e-3*lmin;
    
    vector<int> order(x_.size());
    for(size_t q=0; q<order.size(); ++q)
    order[q]=q;
    
    sort(order.begin(),order.end(),[&](int a, int b){return x_[a](0)<x_[b](0);});
    
    vector<int> newid(x_.size(),-1);
    vector<Eigen::Vector3d> xn;
    
    for(size_t qa=0; qa<order.size(); ++qa)
    {
        const int a = order[qa];
        
        if(newid[a]>=0)
        continue;
        
        newid[a] = xn.size();
        xn.push_back(x_[a]);
        
        for(size_t qb=qa+1; qb<order.size() && x_[order[qb]](0)-x_[a](0)<tol; ++qb)
        {
            const int b = order[qb];
            
            if(newid[b]<0 && (x_[b]-x_[a]).norm()<tol)
            newid[b] = newid[a];
        }
    }
    
    x_ = xn;
    
    for(auto &t : tri_)
    for(int q=0; q<3; ++q)
    t[q] = newid[t[q]];
}

void net_membrane::update_geometry()
{
    // normals, centroids and areas of the (moved) triangles; the vertex order keeps the outward orientation
    for(size_t t=0; t<tri_.size(); ++t)
    {
        const Eigen::Vector3d &a = x_[tri_[t][0]];
        const Eigen::Vector3d &b = x_[tri_[t][1]];
        const Eigen::Vector3d &c = x_[tri_[t][2]];
        
        Eigen::Vector3d n = (b-a).cross(c-a);
        const double A2 = n.norm();
        
        if(A2>1.0e-20)
        {
        tn_[t] = n/A2;
        ta_[t] = 0.5*A2;
        }
        
        tc_[t] = (a+b+c)/3.0;
    }
}

bool net_membrane::inside_footprint(double xp, double yp, double margin) const
{
    xp -= offx_;
    yp -= offy_;

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
    xp -= offx_;
    yp -= offy_;

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

Eigen::Vector3d net_membrane::membrane_vel(const vector<Eigen::Vector3d> &v, int t, const double *w) const
{
    return w[0]*v[tri_[t][0]] + w[1]*v[tri_[t][1]] + w[2]*v[tri_[t][2]];
}

void net_membrane::layer_matrix(const cellentry &e, const vector<Eigen::Vector3d> &v, Eigen::Matrix3d &A, Eigen::Vector3d &b) const
{
    // resistance of a layer cell, f = -(A u - b): each surface orientation with the velocity of its own closest
    // triangle (at the floor edge the floor and the wall move differently), tangential with the closest triangle
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d &nc = tn_[e.tc];
    const Eigen::Matrix3d Pt = Kt*e.Hc*(I - nc*nc.transpose());
    const double wc[3] = {e.w0,e.w1,e.w2};

    A = Pt;
    b = Pt*membrane_vel(v,e.tc,wc);

    for(int r=0; r<e.ns; ++r)
    {
        const Eigen::Vector3d &nr = tn_[e.t[r]];
        const Eigen::Matrix3d Pr = Kn*e.H[r]*(nr*nr.transpose());
        
        A += Pr;
        b += Pr*membrane_vel(v,e.t[r],e.bw[r]);
    }
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

                // inside a floating body (X 10, level set FB < 0): the direct forcing of the body sets the velocity
                // there and the part of the membrane inside the body (the edge clamped in the collar) is carried by
                // it. A membrane resistance in these cells acts against the body forcing: the projected velocity
                // differs from the body velocity before the reforcing, the layer turns that into a load K (u - u_m)
                // on the attached nodes, which goes to the body (3D ring collar: diverged in the first steps)
                if(p->X10>0 && d->FB[IJK]<0.0)
                continue;

                const double zp = p->ZSP[IJK];

                if(zp<bzmin || zp>bzmax)
                continue;

                const Eigen::Vector3d P(p->XP[IP], p->YP[JP], zp);
                const Eigen::Vector3d C = closest_point(P,a,b,c,u,v,w);
                const double dist = (P-C).norm();

                if(dist>=delta)
                continue;
                
                // moving membrane, beyond a corner edge (the cell projects outside the triangle across an edge
                // to a panel of different orientation): the panel ends there. The layers of the two panels meet
                // on the inner side of the corner and seal it; on the outer side the water has to pass freely
                // below the edge of a moving floor, the rounded outer layer would have to be dragged along.
                // A fixed membrane keeps the rounded corner (it steadies the wall columns below the floor edge).
                if(moving() && !cornerE_.empty())
                {
                    const Eigen::Vector3d dP = P-C;
                    const Eigen::Vector3d dtan = dP - dP.dot(tn_[t])*tn_[t];
                    
                    if(dtan.norm()>1.0e-6*delta)
                    {
                        const double eps = 1.0e-9;
                        int q=-1, vtx=-1;
                        
                        if(u<eps && v<eps)      vtx = tri_[t][2];
                        else if(v<eps && w<eps) vtx = tri_[t][0];
                        else if(u<eps && w<eps) vtx = tri_[t][1];
                        else if(w<eps)          q = 0;
                        else if(u<eps)          q = 1;
                        else if(v<eps)          q = 2;
                        
                        if((vtx>=0 && cornerN_[vtx]) || (q>=0 && cornerE_[tedge_[t][q]]))
                        continue;
                    }
                }

                int &s = slot_[IJK];

                if(s<0)
                {
                    s = cells_.size();
                    cellentry e;
                    e.i=i; e.j=j; e.k=k;
                    e.ns=1;
                    e.t[0]=t; e.dd[0]=dist; e.H[0]=0.0;
                    e.bw[0][0]=u; e.bw[0][1]=v; e.bw[0][2]=w;
                    e.tc=t; e.w0=u; e.w1=v; e.w2=w;
                    e.dc=dist;
                    cells_.push_back(e);
                    continue;
                }

                cellentry &e = cells_[s];

                // closest triangle overall: membrane velocity and tangential resistance
                if(dist<e.dc)
                {
                    e.dc=dist; e.tc=t; e.w0=u; e.w1=v; e.w2=w;
                }

                // closest triangle of each distinct surface orientation (up to three: walls meeting
                // the floor at the bottom corners of a box need all three normals)
                // (orientation of the panel in the undeformed mesh: wrinkles of a flexible membrane must not split
                // a wall into several orientations, which would switch the resistance of its cells)
                int q=-1;
                for(int r=0; r<e.ns; ++r)
                if(fabs(tn0_[t].dot(tn0_[e.t[r]]))>=samenormal)
                q=r;

                if(q>=0)
                {
                    if(dist<e.dd[q])
                    {
                        e.t[q]=t;
                        e.dd[q]=dist;
                        e.bw[q][0]=u; e.bw[q][1]=v; e.bw[q][2]=w;
                    }
                }
                else if(e.ns<3)
                {
                    e.t[e.ns]=t;
                    e.dd[e.ns]=dist;
                    e.bw[e.ns][0]=u; e.bw[e.ns][1]=v; e.bw[e.ns][2]=w;
                    ++e.ns;
                }
                else
                {
                    int f=0;
                    for(int r=1; r<3; ++r)
                    if(e.dd[r]>e.dd[f])
                    f=r;

                    if(dist<e.dd[f])
                    {
                        e.t[f]=t;
                        e.dd[f]=dist;
                        e.bw[f][0]=u; e.bw[f][1]=v; e.bw[f][2]=w;
                    }
                }
            }
        }
    }

    for(auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;
        slot_[IJK] = -1;

        e.Hmax=0.0;
        for(int r=0; r<e.ns; ++r)
        {
            e.H[r] = indicator(e.dd[r]);
            e.Hmax = MAX(e.Hmax,e.H[r]);
        }
        
        e.Hc = indicator(e.dc);
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
        d->MBETA[IJK] = MIN(d->MBETA[IJK], 1.0/(1.0 + a*Kn*e.Hmax));
    }
}

void net_membrane::forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    const double a = alpha*p->dt;
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    // node velocities of this forcing (strong coupling: reforce_nhflow adds the response to their change)
    if(iterated())
    xdotf_ = xdot_;

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        // f = -sum_r K_n H_r n_r n_r^T (u - u_m,r) - K_t H_c (I - n_c n_c^T)(u - u_m,c)
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        layer_matrix(e,xdot_,A,b);

        Eigen::Vector3d u(d->U[IJK], d->V[IJK], d->W[IJK]);

        if(p->j_dir==0)
        u(1)=0.0;

        // implicit resistance over the stage: (I + a A) u^new = u^* + a b
        const Eigen::Vector3d unew = (I + a*A).ldlt().solve(u + a*b);
        Eigen::Vector3d du = unew - u;
        
        if(p->j_dir==0)
        du(1)=0.0;

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
    compute_loads(p,d,pgc,WL);

    // structure: the flexible membrane is advanced once per time step, in the final stage, with the loads
    // of the final velocity; the load on the floating body is updated with it and held over the next step.
    // The rigid membrane passes its load to the body in every stage.
    // Strong coupling: the structure was solved together with the fluid in this stage (couple_nhflow), its
    // result is taken over here.
    if(iterated())
    {
        if(pending_)
        stage_commit(p);
    }
    else
    if(prm.structure==2 && finalize)
    {
        sample_collar(p,d,pgc);
        advance_structure(p,p->dt);
        update_geometry();
        update_body_load(p);
        
        if(body_)
        attach_response(p,p->dt);
    }
    
    if(prm.structure==1)
    update_body_load(p);
    
    if(prm.structure==2)
    update_body_fluid_load(p);

    if(finalize)
    {
        collar_step_end(p);
        print_timeseries(p,d,pgc);

        if(print_now(p))
        print_vtp(p);
    }
}

void net_membrane::compute_loads(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    // load on the membrane from the final (projected) velocity: F = rho H A (u^{n+1} - u_m) dV,
    // per triangle and per node (barycentric weights of the closest point)
    fill(tf_.begin(),tf_.end(),0.0);
    fill(nf_.begin(),nf_.end(),0.0);
    
    auto addnode = [&](int t, const double *w, const Eigen::Vector3d &f)
    {
        for(int q=0; q<3; ++q)
        {
            const int nd = tri_[t][q];
            nf_[3*nd+0] += w[q]*f(0);
            nf_[3*nd+1] += w[q]*f(1);
            nf_[3*nd+2] += w[q]*f(2);
        }
    };

    double Q=0.0, urm=0.0;
    const double rho = p->W1;

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        if(p->wet[IJ]==0)
        continue;

        const Eigen::Vector3d um = membrane_vel(e.tc,e.w0,e.w1,e.w2);

        Eigen::Vector3d rel(d->U[IJK]-um(0), d->V[IJK]-um(1), d->W[IJK]-um(2));

        if(p->j_dir==0)
        rel(1)=0.0;

        const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);

        // normal resistance, per surface orientation
        double hun=0.0;
        
        for(int r=0; r<e.ns; ++r)
        {
            const Eigen::Vector3d &nr = tn_[e.t[r]];
            const Eigen::Vector3d umr = moving() ? membrane_vel(e.t[r],e.bw[r][0],e.bw[r][1],e.bw[r][2]) : um;
            const double un = nr.dot(Eigen::Vector3d(d->U[IJK]-umr(0), p->j_dir==1 ? d->V[IJK]-umr(1) : 0.0, d->W[IJK]-umr(2)));
            const Eigen::Vector3d f = rho*dV*e.H[r]*Kn*un*nr;

            tf_[3*e.t[r]+0] += f(0);
            tf_[3*e.t[r]+1] += f(1);
            tf_[3*e.t[r]+2] += f(2);
            
            addnode(e.t[r],e.bw[r],f);
            
            hun += e.H[r]*un;
            
            if(e.H[r]>0.5)
            urm = MAX(urm,fabs(un));
        }
        
        // tangential resistance, closest triangle
        if(Kt>0.0)
        {
            const Eigen::Vector3d &nc = tn_[e.tc];
            const Eigen::Vector3d f = rho*dV*e.Hc*Kt*(rel - nc.dot(rel)*nc);

            tf_[3*e.tc+0] += f(0);
            tf_[3*e.tc+1] += f(1);
            tf_[3*e.tc+2] += f(2);
            
            const double w[3] = {e.w0,e.w1,e.w2};
            addnode(e.tc,w,f);
        }

        // flux through the membrane: int u_n dA = (1/(1.5 delta)) int H u_n dV, outward positive
        Q += hun*dV/(1.5*delta);
    }

    MPI_Allreduce(MPI_IN_PLACE,tf_.data(),(int)tf_.size(),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,nf_.data(),(int)nf_.size(),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

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
}

void net_membrane::netForces(lexer *p, double &Xne, double &Yne, double &Zne, double &Kne, double &Mne, double &Nne)
{
    // load of the membrane on the floating body, moment about the centre of gravity
    body_load(p,Eigen::Vector3d(p->xg,p->yg,p->zg),Rb_,Xne,Yne,Zne,Kne,Mne,Nne);
}

void net_membrane::kinematics_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // node positions and velocities for this stage: with a floating body the rigid membrane moves with it
    // entirely, the flexible one at its attached nodes (the free nodes move in advance_structure)
    if(!moving())
    return;
    
    if(body_)
    for(size_t q=0; q<x_.size(); ++q)
    if(prm.structure==1 || att_[q])
    {
        x_[q] = attached_position(q);
        xdot_[q] = body_velocity(x_[q]);
    }
    
    update_geometry();
    floor_geometry(p);
}

void net_membrane::floor_geometry(lexer *p)
{
    // footprint offset from the floor edge and floor height per column for a moving membrane
    if(!moving() || ring_.empty())
    return;
    
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    for(int q : ring_)
    c += x_[q];
    c /= double(ring_.size());
    
    offx_ = c(0) - ringc0_(0);
    offy_ = p->j_dir==1 ? c(1) - ringc0_(1) : 0.0;
    zring_ = c(2);
    
    fill(zf_.begin(),zf_.end(),zring_);
    
    // floor triangles: height at the column centres inside their horizontal projection
    for(size_t t=0; t<tri_.size(); ++t)
    {
        if(ttag_[t]!=1)
        continue;
        
        const Eigen::Vector3d &a = x_[tri_[t][0]];
        const Eigen::Vector3d &b = x_[tri_[t][1]];
        const Eigen::Vector3d &e = x_[tri_[t][2]];
        
        const double bxmin = min(a(0),min(b(0),e(0))), bxmax = max(a(0),max(b(0),e(0)));
        const double bymin = min(a(1),min(b(1),e(1))), bymax = max(a(1),max(b(1),e(1)));
        
        const int ia = lower_bound(xc_.begin(),xc_.end(),bxmin) - xc_.begin();
        const int ib = int(upper_bound(xc_.begin(),xc_.end(),bxmax) - xc_.begin()) - 1;
        
        int ja=0, jb=0;
        if(p->j_dir==1)
        {
        ja = lower_bound(yc_.begin(),yc_.end(),bymin) - yc_.begin();
        jb = int(upper_bound(yc_.begin(),yc_.end(),bymax) - yc_.begin()) - 1;
        }
        
        const double det = (b(0)-a(0))*(e(1)-a(1)) - (e(0)-a(0))*(b(1)-a(1));
        
        if(fabs(det)<1.0e-20)
        continue;
        
        for(i=ia; i<=ib; ++i)
        for(j=ja; j<=jb; ++j)
        {
            const double xp = p->XP[IP];
            const double yp = p->YP[JP];
            
            const double l1 = ((xp-a(0))*(e(1)-a(1)) - (e(0)-a(0))*(yp-a(1)))/det;
            const double l2 = ((b(0)-a(0))*(yp-a(1)) - (xp-a(0))*(b(1)-a(1)))/det;
            
            if(l1<-1.0e-10 || l2<-1.0e-10 || l1+l2>1.0+1.0e-10)
            continue;
            
            zf_[(i-p->imin)*p->jmax + (j-p->jmin)] = a(2) + l1*(b(2)-a(2)) + l2*(e(2)-a(2));
        }
    }
}

double net_membrane::zfloor(lexer *p, int ii, int jj) const
{
    return moving() ? zf_[(ii-p->imin)*p->jmax + (jj-p->jmin)] : prm.zb;
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
    xp -= offx_;
    yp -= offy_;
    
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

void net_membrane::footprint_weight(double xp, double yp, double &A, double &Ax, double &Ay) const
{
    // Smoothed footprint indicator, 1 inside, 0 outside, with the Darcy profile across the walls.
    // Box: tensor product of the 1D steps across the four walls - smooth at the vertical corners,
    // where a distance function min(x-x0, x1-x, ...) has a kink along the diagonal that leaves a
    // spurious torque below the floor. Cylinder: step in the radial distance.
    Ax=Ay=0.0;
    xp -= offx_;
    yp -= offy_;
    
    if(prm.shape==1)
    {
        const double sx0 = xp-prm.x0, sx1 = prm.x1-xp;
        const double Sx = smoothstep(sx0)*smoothstep(sx1);
        const double dSx = dstep(sx0)*smoothstep(sx1) - smoothstep(sx0)*dstep(sx1);
        
        double Sy=1.0, dSy=0.0;
        
        if(yc_.size()>1)
        {
            const double sy0 = yp-prm.y0, sy1 = prm.y1-yp;
            Sy  = smoothstep(sy0)*smoothstep(sy1);
            dSy = dstep(sy0)*smoothstep(sy1) - smoothstep(sy0)*dstep(sy1);
        }
        
        A  = Sx*Sy;
        Ax = dSx*Sy;
        Ay = Sx*dSy;
        return;
    }
    
    const double rx = xp-prm.xc, ry = yp-prm.yc;
    const double r = sqrt(rx*rx + ry*ry);
    const double s = prm.R - r;
    
    A = smoothstep(s);
    
    if(r>1.0e-12)
    {
        Ax = -dstep(s)*rx/r;
        Ay = -dstep(s)*ry/r;
    }
}

void net_membrane::footprint_weight_ext(double xp, double yp, double &A, double &Ax, double &Ay) const
{
    // footprint weight widened across the whole wall layer: 1 for s > -delta, 0 for s < -2 delta
    auto W  = [this](double s) {return smoothstep(2.0*(s + 1.5*delta));};
    auto dW = [this](double s) {return 2.0*dstep(2.0*(s + 1.5*delta));};
    
    Ax=Ay=0.0;
    xp -= offx_;
    yp -= offy_;
    
    if(prm.shape==1)
    {
        const double sx0 = xp-prm.x0, sx1 = prm.x1-xp;
        const double Sx  = W(sx0)*W(sx1);
        const double dSx = dW(sx0)*W(sx1) - W(sx0)*dW(sx1);
        
        double Sy=1.0, dSy=0.0;
        
        if(yc_.size()>1)
        {
            const double sy0 = yp-prm.y0, sy1 = prm.y1-yp;
            Sy  = W(sy0)*W(sy1);
            dSy = dW(sy0)*W(sy1) - W(sy0)*dW(sy1);
        }
        
        A  = Sx*Sy;
        Ax = dSx*Sy;
        Ay = Sx*dSy;
        return;
    }
    
    const double rx = xp-prm.xc, ry = yp-prm.yc;
    const double r = sqrt(rx*rx + ry*ry);
    const double s = prm.R - r;
    
    A = W(s);
    
    if(r>1.0e-12)
    {
        Ax = -dW(s)*rx/r;
        Ay = -dW(s)*ry/r;
    }
}

void net_membrane::footprint_weight_int(double xp, double yp, double &A, double &Ax, double &Ay) const
{
    // interior weight: 0 for s < delta, 1 for s > 3 delta
    auto W  = [this](double s) {return smoothstep(s - 2.0*delta);};
    auto dW = [this](double s) {return dstep(s - 2.0*delta);};
    
    Ax=Ay=0.0;
    xp -= offx_;
    yp -= offy_;
    
    if(prm.shape==1)
    {
        const double sx0 = xp-prm.x0, sx1 = prm.x1-xp;
        const double Sx  = W(sx0)*W(sx1);
        const double dSx = dW(sx0)*W(sx1) - W(sx0)*dW(sx1);
        
        double Sy=1.0, dSy=0.0;
        
        if(yc_.size()>1)
        {
            const double sy0 = yp-prm.y0, sy1 = prm.y1-yp;
            Sy  = W(sy0)*W(sy1);
            dSy = dW(sy0)*W(sy1) - W(sy0)*dW(sy1);
        }
        
        A  = Sx*Sy;
        Ax = dSx*Sy;
        Ay = Sx*dSy;
        return;
    }
    
    const double rx = xp-prm.xc, ry = yp-prm.yc;
    const double r = sqrt(rx*rx + ry*ry);
    const double s = prm.R - r;
    
    A = W(s);
    
    if(r>1.0e-12)
    {
        Ax = -dW(s)*rx/r;
        Ay = -dW(s)*ry/r;
    }
}

double net_membrane::dindicator(double dist) const
{
    // d H / d dist
    const double h = 0.5*delta;
    
    if(dist<=h || dist>=delta)
    return 0.0;
    
    return -0.5*PI/h*sin(PI*(dist-h)/h);
}

double net_membrane::dstep(double s) const
{
    return fabs(s)<delta ? indicator(fabs(s))/(1.5*delta) : 0.0;
}

void net_membrane::static_pressure_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    // NHFLOW carries one free surface per column, so a column below the bag floor sees the inner
    // water level in its hydrostatic pressure, while the water there is connected to the outside.
    // The non-hydrostatic pressure would have to jump by -rho g dh across the floor, and that sharp
    // jump on the sloping sigma levels at the bag edge leaves pressure gradient errors that drive a
    // spurious circulation under the floor. The static part is therefore prescribed, below the floor
    // (B) and over the footprint (A or W), as
    //
    //   floorpressure 0:  p_m = -rho g dh A B                A = Darcy profile across the walls
    //   floorpressure 1:  p_m = -rho g (eta - eta_out) W B   local excess head, W = 1 for s > -delta
    //   floorpressure 3:  p_m = -rho g dh S W B              S = (eta_avg - eta_out,avg)/dh_avg, the shape
    //                                                         of the free-surface ramp, running average
    //                                                         (tau) frozen at t = 2 tau
    //
    // On a Cartesian grid the discrete free-surface ramp across the wall layer of a curved bag differs
    // between the axis and the diagonal directions; with the analytic profile A (0) the mismatch drives a
    // slowly growing, grid-aligned circulation below the bag edge. 3 matches the actual ramp; it is frozen
    // after the start, since a shape that keeps following the free surface also absorbs slow drifts of the
    // levels in the wall columns, which then lose their restoring pressure below the floor (unstable with
    // A 520 2). 3 is the default for a fixed membrane, 0 for a moving one (the ramp moves with it).
    // 1 follows the instantaneous level and is kept for comparison only.
    // Deep inside the footprint (s > 3 delta) the uniform dh is used. eta_out is the mean level in a
    // ring 1 - 3 delta outside the footprint. The gradient of p_m is added as a body force in the stage.
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
    etaref = eout/aout;
    
    // floorpressure 3: shape of the discrete free-surface ramp, running average (relaxation time tau) until
    // t = 2 tau, frozen afterwards. A shape that keeps following the free surface also absorbs slow drifts of the
    // levels in the wall columns, which then have no restoring pressure below the floor; with A 520 2 these
    // drifts grow.
    if(prm.floorp==3 && p->simtime<=2.0*prm.tau)
    {
        const double r = etab_.empty() ? 1.0 : MIN(1.0, alpha*p->dt/prm.tau);
        
        if(etab_.empty())
        etab_.assign(p->imax*p->jmax,0.0);
        
        for(int ii=-1; ii<=p->knox; ++ii)
        for(int jj=(p->j_dir==1?-1:0); jj<=(p->j_dir==1?p->knoy:0); ++jj)
        {
            double &e = etab_[(ii-p->imin)*p->jmax + (jj-p->jmin)];
            e += r*(d->eta(ii,jj) - e);
        }
        
        erefb += r*(etaref - erefb);
        dhb   += r*(dh - dhb);
    }
    
    const double a = alpha*p->dt;
    const double g = fabs(p->W22);
    
    LOOP
    WETDRYDEEP
    {
        const double sz = zfloor(p,i,j) - p->ZSP[IJK];
        
        if(sz<=-delta)
        continue;
        
        double A, Ax, Ay, fx, fy, fz;
        
        if(prm.floorp==1 || prm.floorp==3)
        {
            // local excess head of the column over the outside level, applied over the footprint and
            // across the whole wall ramp (W = 1 for s > -delta): below the floor the hydrostatic head is
            // then eta_out everywhere, whatever the shape of the discrete free-surface ramp
            footprint_weight_ext(p->XP[IP],p->YP[JP],A,Ax,Ay);
            
            if(A<=0.0)
            continue;
            
            // excess head: local (eta - eta_out) across the wall ramp, blended into the uniform dh
            // deeper inside (I = 1 for s > 3 delta), where the local value would decouple the columns
            // below the floor from the inner free surface and let slow modes grow
            double I, Ix, Iy;
            footprint_weight_int(p->XP[IP],p->YP[JP],I,Ix,Iy);
            
            double E, Ex, Ey;
            
            if(prm.floorp==3)
            {
                // floorpressure 3: shape of the ramp from the averaged level, normalised with the averaged
                // head, S = (eta_avg - eta_ref,avg)/dh_avg; amplitude from the current head dh.
                // The discrete (grid-dependent) ramp across the wall layer is matched as in 1, but neither the
                // amplitude nor the shape feeds back on the local free surface.
                auto EV = [&](int ii, int jj) {return etab_[(ii-p->imin)*p->jmax + (jj-p->jmin)];};
                const double sc = fabs(dhb)>1.0e-8 ? 1.0/dhb : 0.0;
                
                double S  = (EV(i,j) - erefb)*sc;
                double Sx = (EV(i+1,j) - EV(i-1,j))/(p->XP[IP+1] - p->XP[IM1])*sc;
                double Sy = p->j_dir==1 ? (EV(i,j+1) - EV(i,j-1))/(p->YP[JP+1] - p->YP[JM1])*sc : 0.0;
                
                E  = dh*((1.0-I)*S + I);
                Ex = dh*((1.0-I)*Sx + Ix*(1.0 - S));
                Ey = dh*((1.0-I)*Sy + Iy*(1.0 - S));
            }
            else
            {
                const double El = d->eta(i,j) - etaref;
                const double Elx = (d->eta(i+1,j) - d->eta(i-1,j))/(p->XP[IP+1] - p->XP[IM1]);
                const double Ely = p->j_dir==1 ? (d->eta(i,j+1) - d->eta(i,j-1))/(p->YP[JP+1] - p->YP[JM1]) : 0.0;
                
                E  = (1.0-I)*El + I*dh;
                Ex = (1.0-I)*Elx + Ix*(dh - El);
                Ey = (1.0-I)*Ely + Iy*(dh - El);
            }
            const double B  = smoothstep(sz);
            const double dB = sz<delta ? indicator(fabs(sz))/(1.5*delta) : 0.0;
            
            // the continuity dissipation is switched off below the floor inside the wall only (as for 0):
            // the widened weight W reaches into outside columns, whose free surface needs it
            {
                double A0, A0x, A0y;
                footprint_weight(p->XP[IP],p->YP[JP],A0,A0x,A0y);
                d->MCHI[IJK] = MAX(d->MCHI[IJK], A0*B);
            }
            
            // f = -grad(p_m)/rho,  p_m = -rho g (eta - eta_ref) W B
            fx =  g*B*(A*Ex + E*Ax);
            fy =  g*B*(A*Ey + E*Ay);
            fz = -g*E*A*dB;
        }
        else
        {
            // uniform head difference dh over the smoothed footprint
            footprint_weight(p->XP[IP],p->YP[JP],A,Ax,Ay);
            
            if(A<=0.0)
            continue;
            
            const double B  = smoothstep(sz);
            const double dB = sz<delta ? indicator(fabs(sz))/(1.5*delta) : 0.0;
            
            d->MCHI[IJK] = MAX(d->MCHI[IJK], A*B);
            
            // f = -grad(p_m)/rho
            fx =  g*dh*Ax*B;
            fy =  g*dh*Ay*B;
            fz = -g*dh*A*dB;
        }
        
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
