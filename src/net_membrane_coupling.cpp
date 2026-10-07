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
#include"iqn_ils.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<mpi.h>
#include<iostream>
#include<cmath>

// Strong (iterated) coupling of a flexible membrane with NHFLOW, membrane.dat 'coupling iterated'.
//
// A light membrane in water is the textbook case of the added-mass instability of partitioned FSI: the fluid
// load responds to the membrane acceleration with a mass far larger than the fabric (the water of the bag), so
// any lagged exchange of velocities and loads diverges (Causin, Gerbeau, Nobile 2005). The staggered scheme
// (coupling staggered) is stabilised by the implicit porous damper, which costs an extra inertia ~ rho R_n dt.
// Here the interface condition is converged instead, in every Runge-Kutta stage:
//
//   unknown    v: velocities of the free nodes the fluid uses in the porous-jump forcing of the stage
//   fluid      forcing with v, projection, velocity  ->  node loads F(v)            (nhflow_forcing::projection)
//   structure  backward Euler of the stage from the base state with its real inertia,
//                (m + h C)(v~) + h (h S + D)(..) = m v_b + h (F_spring + F(v) + C v + weight + sinker)
//              with h = alpha dt and the base (x_b, v_b) = (1 - alpha) (x, v)^n + alpha (x, v)^{k-1}, the same
//              convex combination as the fluid stage. C, the local response of the layer fluid to the membrane
//              velocity (robin_matrix), is a Robin preconditioner: C (v~ - v) vanishes at convergence.
//   update     IQN-ILS on r = v~ - v (iqn_ils.h), columns of earlier stages of the same RK index reused
//
// The geometry (cell map, normals) is that of the start of the stage and the floating body (6DOF) is moved
// outside the iteration; its load comes from the converged stage. The fluid side is repeated without changes to
// the pressure solver: nhflow_forcing::projection restores the predicted velocity, and the forcing being linear in
// v, reforce_nhflow adds the response to the change of v, (I + a A)^-1 a b(v - v_f).

void net_membrane::reforce_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    if(!iterated())
    return;

    const double a = alpha*p->dt;
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    vector<Eigen::Vector3d> dv(xdot_.size());
    for(size_t q=0; q<xdot_.size(); ++q)
    dv[q] = xdot_[q] - xdotf_[q];

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        layer_matrix(e,dv,A,b);

        Eigen::Vector3d du = (I + a*A).ldlt().solve(a*b);

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

bool net_membrane::couple_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double alpha, slice &WL, int it)
{
    if(!iterated() || gm_.empty())
    return true;

    const size_t nn = x_.size();
    const double h = alpha*p->dt;

    if(it==0)
    {
    sample_collar(p,d,pgc);
    stage_begin(p,pgc,iter,alpha,WL);
    }

    // loads of this iteration
    compute_loads(p,d,pgc,WL);

    // structure: backward Euler of the stage from the base state, attached nodes to their current position
    const vector<Eigen::Vector3d> xg = x_, vg = vs_;

    x_ = xbase_;
    vs_ = vbase_;

    structure_solve(p,h,1,Cr_,xdot_,xbase_,xg,true);

    xr_ = x_;
    vr_ = vs_;
    x_ = xg;
    vs_ = vg;

    // fixed point on the free node velocities
    Eigen::VectorXd xin, xout;
    group_vector(xdot_,xin,p->j_dir==0);
    group_vector(vr_,xout,p->j_dir==0);

    const double rn = (xout - xin).norm();
    const double on = xout.norm();
    const double sn = sqrt(double(MAX(1,(int)xin.size())));
    // the unconverged Robin term C (v~ - v) acts on the structure as a spurious force; with the Robin matrix scaled
    // by f the velocity tolerance is scaled by 1/f, so that this force error does not depend on f
    const double rtol = prm.crtol/prm.crobin, atol = prm.catol/prm.crobin;
    bool conv = rn <= rtol*on + atol*sn;
    
    // the structure is replicated on all ranks; agree on the decision, so that rounding differences of the reductions
    // can never let one rank leave the loop while the others call the pressure solver again
    {
        int c = conv ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE,&c,1,MPI_INT,MPI_MIN,pgc->mpi_comm);
        conv = c==1;
    }

    ++citstep_;

    if(prm.clog==1 && p->mpirank==0)
    cout<<"Membrane "<<nMem<<" coupling: t "<<p->simtime<<" stage "<<iter<<" it "<<it<<"  |r|/|v| "<<rn/MAX(on,1.0e-20)
        <<"  |v|rms "<<on/sn<<"  columns "<<pqn_->columns()<<endl;

    if(conv || it+1>=prm.citer)
    {
        cres_ = MAX(cres_, rn/MAX(on, atol*sn));

        if(!conv)
        {
            ++cwarn_;

            if(p->mpirank==0 && (cwarn_<=10 || cwarn_%100==0))
            cout<<"Membrane "<<nMem<<": coupling not converged in "<<prm.citer<<" iterations, |r|/|v| = "<<rn/MAX(on,1.0e-20)
                <<" (t = "<<p->simtime<<", "<<cwarn_<<" times)"<<endl;
        }

        // keep the secant information of this stage for the next ones
        pqn_->next(xin,xout);
        pqn_->end();

        // the structure result is taken over in reaction_nhflow (stage_commit), after the loads of the final velocity
        // were evaluated with the node velocities this projection used
        pending_ = true;
        return true;
    }

    // safeguard: a residual far above the best one of this stage means the secant model (columns of earlier stages,
    // noise of the fluid solution) misleads the update: drop the history and continue with relaxed steps
    if(it==0 || rn<rbest_)
    rbest_ = rn;
    
    if(it>0 && rn>10.0*rbest_ && rn>atol*sn)
    pqn_->restart();
    
    const Eigen::VectorXd xnext = pqn_->next(xin,xout);
    group_scatter(xnext,xdot_);

    return false;
}

void net_membrane::stage_begin(lexer *p, ghostcell *pgc, int iter, double alpha, slice &WL)
{
    const size_t nn = x_.size();

    if(pqn_==nullptr)
    pqn_ = new iqn_ils(4,prm.crelax,prm.creuse,prm.ccols,prm.cfilt,prm.cqn==1);

    // structure at the start: initial positions at rest
    if(xk_.size()!=nn)
    {
        xk_ = x0_;
        vsk_.assign(nn,Eigen::Vector3d::Zero());
    }

    // first stage: start of the time step
    if(iter==0)
    {
        xn_ = xk_;
        vsn_ = vsk_;
        citstep_ = 0;
        cres_ = 0.0;
    }

    // base state of the stage, the convex combination of the Runge-Kutta stage
    xbase_.resize(nn);
    vbase_.resize(nn);

    for(size_t q=0; q<nn; ++q)
    {
        xbase_[q] = (1.0-alpha)*xn_[q] + alpha*xk_[q];
        vbase_[q] = (1.0-alpha)*vsn_[q] + alpha*vsk_[q];
    }

    robin_matrix(p,pgc,alpha,WL);
    
    // new base geometry, step and Robin matrix: the structure matrix is factorised again in the first iteration
    sfac_ok_ = false;

    pqn_->begin(iter);
}

void net_membrane::stage_commit(lexer *p)
{
    // converged structure of the stage: new positions and velocities, geometry, load on the floating body
    pending_ = false;

    x_ = xr_;
    vs_ = vr_;
    xk_ = xr_;
    vsk_ = vr_;
    
    // the fluid of the next stage starts from the converged node velocities
    for(size_t q=0; q<x_.size(); ++q)
    if(!att_[q])
    xdot_[q] = vr_[q];

    update_geometry();
    update_body_load(p);

    if(body_)
    attach_response(p,p->dt,&Cr_);

    vmax_ = 0.0;
    for(size_t q=0; q<x_.size(); ++q)
    if(!att_[q])
    vmax_ = MAX(vmax_, xdot_[q].norm());

    Tmax_ = 0.0;
    for(size_t e=0; e<edge_.size(); ++e)
    if(grp_[edge_[e][0]]!=grp_[edge_[e][1]])
    Tmax_ = MAX(Tmax_, prm.EA*((x_[edge_[e][1]]-x_[edge_[e][0]]).norm() - L0_[e])/L0_[e]);
}

void net_membrane::robin_matrix(lexer *p, ghostcell *pgc, double alpha, slice &WL)
{
    // local response of the layer fluid to the node velocities, lumped per node: in a layer cell the forcing gives
    // u = (I + a A)^-1 (u* + a b(v)), so the load F = rho dV (A u - b) changes with the membrane velocity by
    // dF/dv_m = -rho dV (I + a A)^-1 A (before the projection), distributed with the barycentric weights
    const size_t nn = x_.size();
    const double a = alpha*p->dt;
    const double rho = p->W1;
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

    vector<double> buf(9*nn,0.0);

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        if(p->wet[IJ]==0)
        continue;

        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        layer_matrix(e,xdot_,A,b);

        const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
        const Eigen::Matrix3d C = rho*dV*(I + a*A).ldlt().solve(A);
        const double w[3] = {e.w0,e.w1,e.w2};

        for(int q=0; q<3; ++q)
        {
            const int nd = tri_[e.tc][q];

            for(int r=0; r<3; ++r)
            for(int c=0; c<3; ++c)
            buf[9*nd+3*r+c] += w[q]*C(r,c);
        }
    }

    MPI_Allreduce(MPI_IN_PLACE,buf.data(),(int)buf.size(),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    Cr_.assign(nn,Eigen::Matrix3d::Zero());

    for(size_t q=0; q<nn; ++q)
    for(int r=0; r<3; ++r)
    for(int c=0; c<3; ++c)
    Cr_[q](r,c) = 0.5*prm.crobin*(buf[9*q+3*r+c] + buf[9*q+3*c+r]);
}

void net_membrane::group_vector(const vector<Eigen::Vector3d> &v, Eigen::VectorXd &out, bool twod) const
{
    // velocity of the free node groups (first node of each group), 2D without the y component
    const int dim = twod ? 2 : 3;
    out.resize(dim*gm_.size());

    for(size_t g=0; g<gm_.size(); ++g)
    {
        const Eigen::Vector3d &u = v[gm_[g][0]];

        if(twod)
        {
            out(2*g+0) = u(0);
            out(2*g+1) = u(2);
        }
        else
        out.segment<3>(3*g) = u;
    }
}

void net_membrane::group_scatter(const Eigen::VectorXd &x, vector<Eigen::Vector3d> &v) const
{
    const bool twod = x.size()==2*(int)gm_.size();

    for(size_t g=0; g<gm_.size(); ++g)
    {
        const Eigen::Vector3d u = twod ? Eigen::Vector3d(x(2*g+0),0.0,x(2*g+1)) : Eigen::Vector3d(x.segment<3>(3*g));

        for(int q : gm_[g])
        v[q] = u;
    }
}
