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
#include"ghostcell.h"
#include<map>
#include<Eigen/Sparse>
#include<algorithm>
#include<cmath>

// Structure of the membrane (membrane.dat: structure fixed | rigid | flexible).
//
// flexible: mass-spring membrane on the triangulation
//   nodes     mass m_A A_i (A_i: a third of the adjacent triangle areas), fabric weight with buoyancy
//             below the still water level, m_A g (1 - rho_w/rho_m); sinker weight w along the floor edge
//   edges     springs k_e = E t sum(A_T)/L0^2 (A_T: the one or two adjacent triangles), tension only with
//             1 % stiffness in compression (wrinkling fabric); dampers c_e = 2 zeta sqrt(k_e m_e)
//   attached  the nodes at the top edge (or above 'attach z') follow the floating body; without floating
//             body they stay in place
//   fluid     node loads of the porous jump (net_membrane::reaction_nhflow)
//
// Coupling. The fluid of the smeared layer is tied to the membrane velocity by the resistance K = R/(1.5 delta),
// i.e. per node the load responds to the membrane velocity as a damper C = rho sum_T A_T/3 (R_n n n^T +
// R_t (I - n n^T)) wherever the fluid can not follow the membrane at once (in-plane stretching of the layer,
// the floor edge, grid-scale undulations). With a fabric of ~1 kg/m^2 an explicit (lagged) use of this load is
// unstable by orders of magnitude, and an added-mass correction does not cure the damper part. The load is
// therefore linearised in the node velocity, F(v) = F^n - C (v - v_f), with v_f the velocity the fluid step
// used, and integrated implicitly with the structure:
//
//   (m I + h C)(v^{n+1}) + h (h S + D)(v_a - v_b) = m v^n + h (F_spring + F^n + C v_f + weight + sinker)
//
// linearised backward Euler with the spring stiffness S and damper matrix D of the edges, sub-steps h of at
// most 0.01 s, the fluid load held over the fluid time step. This is unconditionally stable. The equilibrium
// (F_structure + F_fluid = 0) is not affected; the dynamics carry an additional inertia of about rho R dt per
// area (the membrane can move relative to the fluid of the last step only by the Darcy slip of the net force),
// so fast deformations are damped: the flexible membrane is meant for the deformed shape under current, fill
// and slow motions of the collar; for R dt >> the size of the bag (fluid added mass) the response to waves
// is too stiff. The in-plane modes of the light, stiff fabric are damped by backward Euler.
//
// The fluid sees the step-averaged node velocity (x^{n+1} - x^n)/dt.
// 2D: the nodes on both sides of the single cell row move together, without motion in y.

void net_membrane::ini_structure(lexer *p, ghostcell *pgc)
{
    const size_t nn = x_.size();
    const double rho = p->W1;
    
    x0_ = x_;
    vs_.assign(nn,Eigen::Vector3d::Zero());
    nn_.assign(nn,Eigen::Vector3d(0.0,0.0,1.0));
    
    // node areas
    an_.assign(nn,0.0);
    for(size_t t=0; t<tri_.size(); ++t)
    for(int q=0; q<3; ++q)
    an_[tri_[t][q]] += ta_[t]/3.0;
    
    // edges, adjacent area, wall/floor adjacency
    map<pair<int,int>,int> eid;
    vector<double> eA;
    vector<array<int,2> > etag;
    
    edge_.clear();
    tedge_.assign(tri_.size(),{-1,-1,-1});
    
    for(size_t t=0; t<tri_.size(); ++t)
    for(int q=0; q<3; ++q)
    {
        int a = tri_[t][q], b = tri_[t][(q+1)%3];
        if(a>b)
        swap(a,b);
        
        auto it = eid.find({a,b});
        int e;
        
        if(it==eid.end())
        {
            e = edge_.size();
            eid[{a,b}] = e;
            edge_.push_back({a,b});
            eA.push_back(0.0);
            etag.push_back({0,0});
        }
        else
        e = it->second;
        
        eA[e] += ta_[t];
        etag[e][ttag_[t]] += 1;
        tedge_[t][q] = e;
    }
    
    // corner edges: between two triangles of clearly different orientation (floor edge, vertical box corners)
    {
        vector<vector<int> > et(edge_.size());
        for(size_t t=0; t<tri_.size(); ++t)
        for(int q=0; q<3; ++q)
        et[tedge_[t][q]].push_back(t);
        
        cornerE_.assign(edge_.size(),0);
        cornerN_.assign(nn,0);
        
        for(size_t e=0; e<edge_.size(); ++e)
        if(et[e].size()==2 && fabs(tn_[et[e][0]].dot(tn_[et[e][1]]))<0.9)
        {
            cornerE_[e] = 1;
            cornerN_[edge_[e][0]] = cornerN_[edge_[e][1]] = 1;
        }
    }
    
    // floor edge (edges between a floor and a wall triangle) and floor nodes
    vector<char> onring(nn,0), onfloor(nn,0);
    ws_.assign(nn,0.0);
    
    for(size_t e=0; e<edge_.size(); ++e)
    if(etag[e][0]>0 && etag[e][1]>0)
    {
        const int a = edge_[e][0], b = edge_[e][1];
        const double L = (x_[b]-x_[a]).norm();
        
        onring[a] = onring[b] = 1;
        ws_[a] += 0.5*prm.sinker*L;
        ws_[b] += 0.5*prm.sinker*L;
    }
    
    for(size_t t=0; t<tri_.size(); ++t)
    if(ttag_[t]==1)
    for(int q=0; q<3; ++q)
    onfloor[tri_[t][q]] = 1;
    
    ring_.clear();
    flo_.clear();
    ringc0_ = Eigen::Vector3d::Zero();
    
    for(size_t q=0; q<nn; ++q)
    {
        if(onring[q])
        {
        ring_.push_back(q);
        ringc0_ += x_[q];
        }
        
        if(onfloor[q])
        flo_.push_back(q);
    }
    
    if(!ring_.empty())
    ringc0_ /= double(ring_.size());
    
    // attached nodes: top edge
    double lmin = 1.0e20;
    for(size_t e=0; e<edge_.size(); ++e)
    lmin = MIN(lmin, (x_[edge_[e][1]]-x_[edge_[e][0]]).norm());
    
    const double zatt = prm.zattach>-1.0e19 ? prm.zattach : prm.zt;
    
    att_.assign(nn,0);
    for(size_t q=0; q<nn; ++q)
    if(x_[q](2) >= zatt - 1.0e-3*lmin)
    att_[q] = 1;
    
    // groups: in 2D the two nodes across the cell row
    grp_.resize(nn);
    for(size_t q=0; q<nn; ++q)
    grp_[q] = q;
    
    if(p->j_dir==0)
    {
        map<pair<long,long>,int> key;
        
        for(size_t q=0; q<nn; ++q)
        {
            const pair<long,long> k(lround(x_[q](0)/(1.0e-3*lmin)), lround(x_[q](2)/(1.0e-3*lmin)));
            auto it = key.find(k);
            
            if(it==key.end())
            key[k] = q;
            else
            grp_[q] = it->second;
        }
    }
    
    {
        map<int,int> gi;
        gm_.clear();
        
        for(size_t q=0; q<nn; ++q)
        {
            if(att_[q])
            continue;
            
            auto it = gi.find(grp_[q]);
            
            if(it==gi.end())
            {
                gi[grp_[q]] = gm_.size();
                gm_.push_back({(int)q});
            }
            else
            gm_[it->second].push_back(q);
        }
    }
    
    // masses, springs, dampers
    mn_.resize(nn);
    for(size_t q=0; q<nn; ++q)
    mn_[q] = prm.mA*an_[q];
    
    L0_.resize(edge_.size());
    ke_.resize(edge_.size());
    ce_.resize(edge_.size());
    
    for(size_t e=0; e<edge_.size(); ++e)
    {
        const int a = edge_[e][0], b = edge_[e][1];
        
        L0_[e] = (x_[b]-x_[a]).norm();
        ke_[e] = prm.EA*eA[e]/(L0_[e]*L0_[e]);
        ce_[e] = 2.0*prm.zeta*sqrt(ke_[e]*0.5*(mn_[a]+mn_[b]));
    }
    
    xan_ = x_;
}

Eigen::Vector3d net_membrane::attached_position(int q) const
{
    return body_ ? Eigen::Vector3d(cb_ + Rb_*rb_[q]) : x0_[q];
}

Eigen::Vector3d net_membrane::body_velocity(const Eigen::Vector3d &x) const
{
    return body_ ? Eigen::Vector3d(vb_ + wb_.cross(x - cb_)) : Eigen::Vector3d::Zero();
}

void net_membrane::attach_body(lexer *p, const Eigen::Vector3d &c, const Eigen::Matrix3d &R)
{
    // node positions in the body frame at the start
    rb_.resize(x0_.size());
    
    for(size_t q=0; q<x0_.size(); ++q)
    rb_[q] = R.transpose()*(x0_[q] - c);
    
    body_ = true;
    cb_ = c;
    Rb_ = R;
}

void net_membrane::set_body(const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Vector3d &v, const Eigen::Vector3d &w)
{
    cb_ = c;
    Rb_ = R;
    vb_ = v;
    wb_ = w;
}

void net_membrane::internal_forces(lexer *p, const vector<Eigen::Vector3d> &x, const vector<Eigen::Vector3d> &v, vector<Eigen::Vector3d> &F) const
{
    for(auto &f : F)
    f.setZero();
    
    for(size_t e=0; e<edge_.size(); ++e)
    {
        const int a = edge_[e][0], b = edge_[e][1];
        
        if(grp_[a]==grp_[b])
        continue;
        
        const Eigen::Vector3d dx = x[b]-x[a];
        const double L = dx.norm();
        
        if(L<1.0e-20)
        continue;
        
        const Eigen::Vector3d ev = dx/L;
        
        double T = ke_[e]*(L - L0_[e]);
        
        // fabric: compression wrinkles the membrane
        if(T<0.0)
        T *= 0.01;
        
        T += ce_[e]*(v[b]-v[a]).dot(ev);
        
        F[a] += T*ev;
        F[b] -= T*ev;
    }
}

void net_membrane::external_forces(lexer *p, vector<Eigen::Vector3d> &F) const
{
    // fluid load, weight (buoyancy below the still water level), sinker
    const double g = fabs(p->W22);
    const double rho = p->W1;
    
    for(size_t q=0; q<x_.size(); ++q)
    {
        const double wet = x_[q](2) < p->wd ? 1.0 - rho/prm.rhom : 1.0;
        
        F[q] = Eigen::Vector3d(nf_[3*q+0], nf_[3*q+1], nf_[3*q+2]);
        F[q](2) -= mn_[q]*g*wet + ws_[q];
    }
}

void net_membrane::advance_structure(lexer *p, double dt)
{
    if(dt<=0.0 || gm_.empty())
    return;
    
    const size_t nn = x_.size();
    const int ng = gm_.size();
    const double rho = p->W1;
    
    // linearised backward Euler with the implicit porous-jump coupling, sub-steps of at most dtmax
    const double dtmax = 0.01;
    nsub_ = MAX(1, (int)ceil(dt/dtmax - 1.0e-10));
    const double h = dt/double(nsub_);
    
    // node normals (area weighted) and added-mass tensors
    vector<Eigen::Matrix3d> Cq(nn,Eigen::Matrix3d::Zero());
    
    for(auto &n : nn_)
    n.setZero();
    
    // porous-jump coupling per node, dF/dv = -C: C = rho sum_T A_T/3 (R_n n n^T + R_t (I - n n^T))
    for(size_t t=0; t<tri_.size(); ++t)
    for(int q=0; q<3; ++q)
    {
        const Eigen::Matrix3d P = tn_[t]*tn_[t].transpose();
        
        nn_[tri_[t][q]] += ta_[t]*tn_[t];
        Cq[tri_[t][q]] += rho*ta_[t]/3.0*(prm.Rn*P + prm.Rt*(Eigen::Matrix3d::Identity() - P));
    }
    
    for(auto &n : nn_)
    {
        const double l = n.norm();
        n = l>1.0e-20 ? Eigen::Vector3d(n/l) : Eigen::Vector3d(0.0,0.0,1.0);
    }
    
    // attached nodes: from their positions at the start of the step to the current body position
    vector<Eigen::Vector3d> xa1(nn);
    
    for(size_t q=0; q<nn; ++q)
    if(att_[q])
    xa1[q] = attached_position(q);
    
    const vector<Eigen::Vector3d> xs = x_;
    
    // implicit coupling: F_fluid(v) = F^n - C (v - v_seen), v_seen the velocity the fluid step used
    structure_solve(p,dt,nsub_,Cq,xdot_,xan_,xa1);
    
    // end of the step
    vmax_ = 0.0;
    
    for(size_t q=0; q<nn; ++q)
    {
        if(att_[q])
        {
            xdot_[q] = body_velocity(x_[q]);
            xan_[q] = x_[q];
        }
        else
        {
            // the fluid sees the step-averaged velocity
            xdot_[q] = (x_[q] - xs[q])/dt;
            vmax_ = MAX(vmax_, xdot_[q].norm());
        }
    }
    
    // largest membrane tension, E t strain
    Tmax_ = 0.0;
    for(size_t e=0; e<edge_.size(); ++e)
    if(grp_[edge_[e][0]]!=grp_[edge_[e][1]])
    Tmax_ = MAX(Tmax_, prm.EA*((x_[edge_[e][1]]-x_[edge_[e][0]]).norm() - L0_[e])/L0_[e]);
}

struct structure_factor
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
};

void net_membrane::free_factor()
{
    delete sfac_;
    sfac_ = nullptr;
    sfac_ok_ = false;
}

void net_membrane::structure_solve(lexer *p, double dt, int nsub, const vector<Eigen::Matrix3d> &Cq, const vector<Eigen::Vector3d> &vseen,
                                   const vector<Eigen::Vector3d> &xa0, const vector<Eigen::Vector3d> &xa1, bool reuse)
{
    // reuse (single step only): the matrix depends on the base geometry, dt and Cq, not on the loads or vseen; it is
    // factorised at the first call after sfac_ok_ was reset and reused until the next reset (coupling iterations)
    // linearised backward Euler of the mass-spring membrane from x_, vs_ over dt in nsub sub-steps, fluid load
    // nf_ held, linearised in the node velocity with the matrices Cq about vseen; the attached nodes move from
    // xa0 to xa1. Result in x_, vs_ (attached nodes: xa1 and their mean velocity).
    const size_t nn = x_.size();
    const int ng = gm_.size();
    const double h = dt/double(nsub);
    
    vector<int> gid(nn,-1);
    vector<Eigen::Matrix3d> Mg(ng,Eigen::Matrix3d::Zero());
    
    for(int g=0; g<ng; ++g)
    for(int q : gm_[g])
    {
        gid[q] = g;
        Mg[g] += mn_[q]*Eigen::Matrix3d::Identity() + h*Cq[q];
    }
    
    // attached nodes: from xa0 to xa1 over the step
    vector<Eigen::Vector3d> va(nn,Eigen::Vector3d::Zero());
    
    for(size_t q=0; q<nn; ++q)
    if(att_[q])
    va[q] = (xa1[q] - xa0[q])/dt;
    
    // forces held over the step: fluid, weight, sinker, and the added-mass force of the last step
    vector<Eigen::Vector3d> Fc(nn);
    external_forces(p,Fc);
    
    // linearised fluid load: F_fluid(v) = F^n - C (v - v_seen)
    for(size_t q=0; q<nn; ++q)
    if(!att_[q])
    Fc[q] += Cq[q]*vseen[q];
    
    const bool twod = p->j_dir==0;
    
    typedef Eigen::Triplet<double> Trip;
    vector<Trip> trip;
    Eigen::VectorXd rhs(3*ng), sol(3*ng);
    
    auto addblock = [&](int ga, int gb, const Eigen::Matrix3d &B)
    {
        for(int r=0; r<3; ++r)
        for(int c=0; c<3; ++c)
        if(!(twod && (r==1 || c==1)))
        trip.push_back(Trip(3*ga+r,3*gb+c,B(r,c)));
    };
    
    for(int s=0; s<nsub; ++s)
    {
        const double r0 = double(s)/double(nsub);
        
        for(size_t q=0; q<nn; ++q)
        if(att_[q])
        {
            x_[q] = xa0[q] + r0*(xa1[q] - xa0[q]);
            vs_[q] = va[q];
        }
        
        trip.clear();
        rhs.setZero();
        
        for(int g=0; g<ng; ++g)
        {
            addblock(g,g,Mg[g]);
            
            Eigen::Vector3d b = Eigen::Vector3d::Zero();
            for(int q : gm_[g])
            b += mn_[q]*vs_[q] + h*Fc[q];
            
            rhs.segment<3>(3*g) += b;
            
            if(twod)
            trip.push_back(Trip(3*g+1,3*g+1,1.0));
        }
        
        for(size_t e=0; e<edge_.size(); ++e)
        {
            const int a = edge_[e][0], b = edge_[e][1];
            
            if(grp_[a]==grp_[b])
            continue;
            
            const int ga = gid[a], gb = gid[b];
            
            if(ga<0 && gb<0)
            continue;
            
            const Eigen::Vector3d dx = x_[b]-x_[a];
            const double L = dx.norm();
            
            if(L<1.0e-20)
            continue;
            
            const Eigen::Vector3d ev = dx/L;
            const Eigen::Matrix3d P = ev*ev.transpose();
            
            // spring force on a at x^n, tension only (1 % stiffness in compression)
            const double kef = L>L0_[e] ? ke_[e] : 0.01*ke_[e];
            const double geo = L>L0_[e] ? 1.0 - L0_[e]/L : 0.0;
            const Eigen::Vector3d fa = kef*(L - L0_[e])*ev;
            
            const Eigen::Matrix3d B = h*h*kef*(P + geo*(Eigen::Matrix3d::Identity() - P)) + h*ce_[e]*P;
            
            if(ga>=0)
            {
                addblock(ga,ga,B);
                rhs.segment<3>(3*ga) += h*fa;
            }
            
            if(gb>=0)
            {
                addblock(gb,gb,B);
                rhs.segment<3>(3*gb) -= h*fa;
            }
            
            if(ga>=0 && gb>=0)
            {
                addblock(ga,gb,-B);
                addblock(gb,ga,-B);
            }
            else if(ga>=0)
            rhs.segment<3>(3*ga) += B*vs_[b];
            
            else
            rhs.segment<3>(3*gb) += B*vs_[a];
        }
        
        if(twod)
        for(int g=0; g<ng; ++g)
        rhs(3*g+1) = 0.0;
        
        if(reuse && nsub==1)
        {
            if(!sfac_ok_)
            {
                Eigen::SparseMatrix<double> A(3*ng,3*ng);
                A.setFromTriplets(trip.begin(),trip.end());
                
                if(sfac_==nullptr)
                sfac_ = new structure_factor;
                
                sfac_->solver.compute(A);
                sfac_ok_ = true;
            }
            
            sol = sfac_->solver.solve(rhs);
        }
        else
        {
        Eigen::SparseMatrix<double> A(3*ng,3*ng);
        A.setFromTriplets(trip.begin(),trip.end());
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        solver.compute(A);
        sol = solver.solve(rhs);
        }
        
        for(int g=0; g<ng; ++g)
        {
            Eigen::Vector3d v = sol.segment<3>(3*g);
            
            if(twod)
            v(1) = 0.0;
            
            for(int q : gm_[g])
            {
                vs_[q] = v;
                x_[q] += h*v;
            }
        }
    }
    
    for(size_t q=0; q<nn; ++q)
    if(att_[q])
    {
        x_[q] = xa1[q];
        vs_[q] = va[q];
    }
}

void net_membrane::update_body_load(lexer *p)
{
    // load of the membrane on the floating body: the rigid membrane passes its whole load (fluid and
    // weight), the flexible one the forces on its attached nodes (edge forces and weight here, at the membrane
    // step; the fluid load on the attached nodes in every stage, update_body_fluid_load)
    const size_t nn = x_.size();
    vector<Eigen::Vector3d> Fe(nn), Fi(nn,Eigen::Vector3d::Zero());
    
    external_forces(p,Fe);
    
    if(prm.structure==2)
    {
        internal_forces(p,x_,vs_,Fi);
        
        for(size_t q=0; q<nn; ++q)
        Fe[q] -= Eigen::Vector3d(nf_[3*q+0],nf_[3*q+1],nf_[3*q+2]);
    }
    
    Fb_.setZero();
    Mb_.setZero();
    
    for(size_t q=0; q<nn; ++q)
    if(prm.structure==1 || att_[q])
    {
        const Eigen::Vector3d F = Fe[q] + Fi[q];
        
        Fb_ += F;
        Mb_ += x_[q].cross(F);
    }
}

void net_membrane::update_body_fluid_load(lexer *p)
{
    // flexible membrane: fluid load on the attached nodes of this stage. The porous jump ties the fluid at the
    // attached edge to the body velocity in every stage; the fluid passes that push on to the hull as pressure
    // within the stage, so its reaction on the attached edge has to reach the body in the same stage, or the
    // exchange is one-sided (it drove the collar away in the tests).
    Ffl_.setZero();
    Mfl_.setZero();
    
    if(prm.structure!=2)
    return;
    
    for(size_t q=0; q<x_.size(); ++q)
    if(att_[q])
    {
        const Eigen::Vector3d F(nf_[3*q+0],nf_[3*q+1],nf_[3*q+2]);
        
        Ffl_ += F;
        Mfl_ += x_[q].cross(F);
    }
}

void net_membrane::body_load(lexer *p, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, double &X, double &Y, double &Z, double &K, double &M, double &N) const
{
    // load on the body at position c, orientation R. The flexible membrane is advanced once per time step; within
    // the step its load responds to the body motion q = (c - c^n, theta) through the impedance of the attached
    // edge (attach_response): F = F^n + J q. Held constant instead, the stiff top edge acts on the light collar
    // like a spring integrated with forward Euler, which grows in every step.
    X=Y=Z=K=M=N=0.0;
    
    if(!moving() || !body_)
    return;
    
    Eigen::Matrix<double,6,1> q;
    q.head<3>() = c - cbn_;
    
    const Eigen::Matrix3d dR = R*Rbn_.transpose();
    q(3) = 0.5*(dR(2,1) - dR(1,2));
    q(4) = 0.5*(dR(0,2) - dR(2,0));
    q(5) = 0.5*(dR(1,0) - dR(0,1));
    
    const Eigen::Matrix<double,6,1> dF = Jb_*q;
    
    const Eigen::Vector3d F  = Fb_ + Ffl_ + dF.head<3>();
    const Eigen::Vector3d Mo = Mb_ + Mfl_ + dF.tail<3>();
    const Eigen::Vector3d Mc = Mo - c.cross(F);
    
    X = F(0);
    Y = F(1);
    Z = F(2);
    K = Mc(0);
    M = Mc(1);
    N = Mc(2);
}

void net_membrane::attach_response(lexer *p, double dt, const vector<Eigen::Matrix3d> *Cext)
{
    // Impedance of the flexible membrane at its attached edge: response of the load on the body to a body motion
    // over the next time step, from the same linearised backward-Euler system as advance_structure (one step dt,
    // porous-jump coupling included). For each generalised displacement q_k (3 translations, 3 rotations about
    // the centre of the body) the attached nodes move by d_a, the free nodes respond with A dv = B d_a/dt, and
    //   dF_a = (dt S + D)(dv_b - d_a/dt)        per edge at an attached node a
    // gives the column k of J = d(F, M_origin)/dq. Used by body_load in the RK stages of the 6DOF.
    Jb_.setZero();
    cbn_ = cb_;
    Rbn_ = Rb_;
    
    if(prm.structure!=2 || gm_.empty() || dt<=0.0)
    return;
    
    const size_t nn = x_.size();
    const int ng = gm_.size();
    const double rho = p->W1;
    const double h = dt;
    const bool twod = p->j_dir==0;
    
    vector<Eigen::Matrix3d> Cq(nn,Eigen::Matrix3d::Zero());
    
    // fluid response of the free nodes: given (strong coupling: local response of the layer fluid),
    // otherwise the implicit porous damper of advance_structure
    if(Cext!=nullptr)
    Cq = *Cext;
    
    else
    for(size_t t=0; t<tri_.size(); ++t)
    for(int q=0; q<3; ++q)
    {
        const Eigen::Matrix3d P = tn_[t]*tn_[t].transpose();
        Cq[tri_[t][q]] += rho*ta_[t]/3.0*(prm.Rn*P + prm.Rt*(Eigen::Matrix3d::Identity() - P));
    }
    
    vector<int> gid(nn,-1);
    vector<Eigen::Matrix3d> Mg(ng,Eigen::Matrix3d::Zero());
    
    for(int g=0; g<ng; ++g)
    for(int q : gm_[g])
    {
        gid[q] = g;
        Mg[g] += mn_[q]*Eigen::Matrix3d::Identity() + h*Cq[q];
    }
    
    typedef Eigen::Triplet<double> Trip;
    vector<Trip> trip;
    
    auto addblock = [&](int ga, int gb, const Eigen::Matrix3d &B)
    {
        for(int r=0; r<3; ++r)
        for(int c=0; c<3; ++c)
        if(!(twod && (r==1 || c==1)))
        trip.push_back(Trip(3*ga+r,3*gb+c,B(r,c)));
    };
    
    for(int g=0; g<ng; ++g)
    {
        addblock(g,g,Mg[g]);
        
        if(twod)
        trip.push_back(Trip(3*g+1,3*g+1,1.0));
    }
    
    // edge matrices B = dt (dt S + D)
    vector<Eigen::Matrix3d> Be(edge_.size(),Eigen::Matrix3d::Zero());
    
    for(size_t e=0; e<edge_.size(); ++e)
    {
        const int a = edge_[e][0], b = edge_[e][1];
        
        if(grp_[a]==grp_[b])
        continue;
        
        const Eigen::Vector3d dx = x_[b]-x_[a];
        const double L = dx.norm();
        
        if(L<1.0e-20)
        continue;
        
        const Eigen::Vector3d ev = dx/L;
        const Eigen::Matrix3d P = ev*ev.transpose();
        const double kef = L>L0_[e] ? ke_[e] : 0.01*ke_[e];
        const double geo = L>L0_[e] ? 1.0 - L0_[e]/L : 0.0;
        
        Be[e] = h*h*kef*(P + geo*(Eigen::Matrix3d::Identity() - P)) + h*ce_[e]*P;
        
        const int ga = gid[a], gb = gid[b];
        
        if(ga>=0)
        addblock(ga,ga,Be[e]);
        
        if(gb>=0)
        addblock(gb,gb,Be[e]);
        
        if(ga>=0 && gb>=0)
        {
            addblock(ga,gb,-Be[e]);
            addblock(gb,ga,-Be[e]);
        }
    }
    
    Eigen::SparseMatrix<double> A(3*ng,3*ng);
    A.setFromTriplets(trip.begin(),trip.end());
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    
    if(solver.info()!=Eigen::Success)
    return;
    
    for(int k=0; k<6; ++k)
    {
        if(twod && (k==1 || k==3 || k==5))
        continue;
        
        // velocity of the attached nodes for a unit generalised displacement over the step
        vector<Eigen::Vector3d> va(nn,Eigen::Vector3d::Zero());
        
        for(size_t q=0; q<nn; ++q)
        if(att_[q])
        {
            Eigen::Vector3d d = Eigen::Vector3d::Zero();
            
            if(k<3)
            d(k) = 1.0;
            else
            {
                Eigen::Vector3d w = Eigen::Vector3d::Zero();
                w(k-3) = 1.0;
                d = w.cross(x_[q] - cb_);
            }
            
            va[q] = d/h;
        }
        
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(3*ng);
        
        for(size_t e=0; e<edge_.size(); ++e)
        {
            const int a = edge_[e][0], b = edge_[e][1];
            const int ga = gid[a], gb = gid[b];
            
            if(ga>=0 && gb<0)
            rhs.segment<3>(3*ga) += Be[e]*va[b];
            
            if(gb>=0 && ga<0)
            rhs.segment<3>(3*gb) += Be[e]*va[a];
        }
        
        if(twod)
        for(int g=0; g<ng; ++g)
        rhs(3*g+1) = 0.0;
        
        const Eigen::VectorXd sol = solver.solve(rhs);
        
        // node velocities of the response
        vector<Eigen::Vector3d> dv(va);
        for(int g=0; g<ng; ++g)
        for(int q : gm_[g])
        dv[q] = sol.segment<3>(3*g);
        
        // load change on the attached nodes
        Eigen::Vector3d dF = Eigen::Vector3d::Zero(), dM = Eigen::Vector3d::Zero();
        
        for(size_t e=0; e<edge_.size(); ++e)
        {
            const int a = edge_[e][0], b = edge_[e][1];
            
            if(grp_[a]==grp_[b] || (!att_[a] && !att_[b]))
            continue;
            
            const Eigen::Vector3d fa = Be[e]/h*(dv[b] - dv[a]);
            
            if(att_[a])
            {
                dF += fa;
                dM += x_[a].cross(fa);
            }
            
            if(att_[b])
            {
                dF -= fa;
                dM -= x_[b].cross(fa);
            }
        }
        
        Jb_.block<3,1>(0,k) = dF;
        Jb_.block<3,1>(3,k) = dM;
    }
}

double net_membrane::body_addedmass(lexer *p) const
{
    // added mass for the stabilised coupling to the floating body: 2 rho V_bag, the water of the bag below the still
    // water level. It moves with a rigid bag; a flexible bag carried along by its top edge passes its inertia to the
    // collar through the pressure under the pontoons (and the attached edge), a comparable outer mass is displaced.
    if(!moving() || !body_)
    return 0.0;
    
    if(prm.Mbody>=0.0)
    return prm.Mbody;
    
    // flexible: the inertia of the bag reaches the collar through its top edge, which is coupled implicitly
    // (attach_response); no added mass by default
    if(prm.structure==2)
    return 0.0;
    

    double A = prm.shape==2 ? PI*prm.R*prm.R : (prm.x1-prm.x0)*(p->j_dir==1 ? prm.y1-prm.y0 : p->DYN[0+marge]);
    
    return 2.0*p->W1*A*MAX(0.0, p->wd - prm.zb);
}
