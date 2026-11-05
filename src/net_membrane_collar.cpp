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
#include<mpi.h>
#include<cmath>
#include<iostream>
#include<algorithm>

// Flexible collar of a flexible bag (membrane.dat 'collar D m EA EI [Cd [Ca]]', 'mooring xa ya za k T0').
//
// The top edge of the bag is a floating pipe ring instead of being held by a floating body: its nodes are free and
// carry, in addition to the fabric, the collar loads below, solved with the fabric in the same linearised
// backward-Euler step (structure_solve) and so, with 'coupling iterated', strongly coupled to the fluid together
// with the bag.
//
//   pipe       mass m l, weight; axial springs EA/L0 between neighbouring ring nodes (tension and compression)
//   bending    corotational: the ring is compared with its reference shape after the best-fit rigid motion
//              (Kabsch), F = -k_b D^T D (x - c - R x_ref), D the second difference along the ring, k_b = EI/l^3;
//              the stiffness k_b D^T D is constant, so large rigid motions are exact and the bending is linear
//              (torsion neglected)
//   fluid      Morison per node, normal to the pipe axis (P_n = I - t t^T), with the submerged part of the
//              circular section below the local free surface eta (segment area A_s, waterline breadth b):
//                buoyancy   rho g A_s l e_z, implicit hydrostatic stiffness rho g b l
//                inertia    rho (1 + C_a) A_s l P_n a_f (Froude-Krylov and added mass of the fluid acceleration),
//                           added mass rho C_a A_s l P_n on the node (implicit)
//                drag       1/2 rho C_d s l |u_r| u_r, u_r = P_n (u_f - v), s the submerged depth (linearised)
//              u_f and eta are sampled from the flow at the collar nodes (in every stage), a_f from the change of
//              u_f over the last time step. The collar does not act back on the flow (slender pipe, D/L << 1):
//              its node sits on the edge of the membrane layer, whose fluid the porous jump ties to the membrane
//   mooring    linear springs from fixed anchors to the nearest collar node, T = max(T0 + k (L - L0), 0)

void net_membrane::ini_collar(lexer *p)
{
    cn_.clear();
    
    if(!collar())
    return;
    
    if(p->j_dir==0 || p->X10>0)
    {
        if(p->mpirank==0)
        cout<<"\n!!! membrane "<<nMem<<": 'collar' needs a 3D simulation and no floating body (X 10 0) !!!\n"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    // the top edge nodes become the collar; they are free
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    
    for(size_t q=0; q<att_.size(); ++q)
    if(att_[q])
    {
        cn_.push_back(q);
        c += x_[q];
        att_[q] = 0;
    }
    
    if(cn_.size()<4)
    {
        if(p->mpirank==0)
        cout<<"\n!!! membrane "<<nMem<<": collar with "<<cn_.size()<<" nodes !!!\n"<<endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    c /= double(cn_.size());
    
    // ring order: angle about the centroid
    sort(cn_.begin(),cn_.end(),[&](int a, int b)
    {
        return atan2(x_[a](1)-c(1),x_[a](0)-c(0)) < atan2(x_[b](1)-c(1),x_[b](0)-c(0));
    });
    
    const int n = cn_.size();
    cL0_.resize(n);
    cl_.resize(n);
    cxr_.resize(n);
    
    for(int i=0; i<n; ++i)
    {
        cL0_[i] = (x_[cn_[(i+1)%n]] - x_[cn_[i]]).norm();
        cxr_[i] = x_[cn_[i]] - c;
    }
    
    for(int i=0; i<n; ++i)
    cl_[i] = 0.5*(cL0_[(i+n-1)%n] + cL0_[i]);
    
    cuf_.assign(n,Eigen::Vector3d::Zero());
    cuf0_.assign(n,Eigen::Vector3d::Zero());
    caf_.assign(n,Eigen::Vector3d::Zero());
    ceta_.assign(n,p->wd);
    ctn_ = -1.0;
    
    // mooring: fairlead at the collar node nearest to the anchor (horizontal distance)
    mfl_.clear();
    mL0_.clear();
    
    for(auto &m : prm.moor)
    {
        int best=0;
        double dmin=1.0e20;
        
        for(int i=0; i<n; ++i)
        {
            const double dd = hypot(x_[cn_[i]](0)-m[0], x_[cn_[i]](1)-m[1]);
            
            if(dd<dmin)
            {
                dmin = dd;
                best = i;
            }
        }
        
        mfl_.push_back(best);
        mL0_.push_back((Eigen::Vector3d(m[0],m[1],m[2]) - x_[cn_[best]]).norm());
    }
    
    mT_.assign(prm.moor.size(),0.0);
    
    if(p->mpirank==0)
    {
        // floating equilibrium of the pipe alone: submerged area m / rho
        const double r = 0.5*prm.cD, Asub = prm.cm/p->W1;
        double s = 0.0;
        
        if(Asub >= PI*r*r)
        s = prm.cD;
        else
        {
            double lo=0.0, hi=prm.cD;
            
            for(int it=0; it<60; ++it)
            {
                s = 0.5*(lo+hi);
                const double th = 2.0*acos((r-s)/r);
                (r*r*0.5*(th - sin(th)) < Asub ? lo : hi) = s;
            }
        }
        
        cout<<"Membrane "<<nMem<<": collar ring of "<<n<<" nodes, D = "<<prm.cD<<" m, m = "<<prm.cm<<" kg/m, EA = "<<prm.cEA
            <<" N, EI = "<<prm.cEI<<" N m^2, Cd = "<<prm.cCd<<", Ca = "<<prm.cCa<<"; pipe alone floats with its centre at "
            <<p->wd + r - s<<" m (draft "<<s<<" m), collar centre line at "<<prm.zt<<" m; "<<prm.moor.size()<<" mooring springs"<<endl;
    }
}

void net_membrane::sample_collar(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // fluid velocity and free surface at the collar nodes; every rank samples the nodes in its subdomain
    if(cn_.empty())
    return;
    
    const int n = cn_.size();
    vector<double> buf(5*n,0.0);
    
    const double xs = p->XN[0+marge], xe = p->XN[p->knox+marge];
    const double ys = p->YN[0+marge], ye = p->YN[p->knoy+marge];
    
    for(int i=0; i<n; ++i)
    {
        const Eigen::Vector3d &x = x_[cn_[i]];
        
        if(x(0)<xs || x(0)>=xe || x(1)<ys || x(1)>=ye)
        continue;
        
        const double bed = p->ccslipol4(d->bed,x(0),x(1));
        const double eta = p->ccslipol4(d->WL,x(0),x(1)) + bed;
        const double z = MAX(bed + 1.0e-3, MIN(x(2), eta - 1.0e-3));
        
        buf[5*i+0] = p->ccipol4V(d->U,d->WL,d->bed,x(0),x(1),z);
        buf[5*i+1] = p->ccipol4V(d->V,d->WL,d->bed,x(0),x(1),z);
        buf[5*i+2] = p->ccipol4V(d->W,d->WL,d->bed,x(0),x(1),z);
        buf[5*i+3] = eta;
        buf[5*i+4] = 1.0;
    }
    
    MPI_Allreduce(MPI_IN_PLACE,buf.data(),5*n,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    
    for(int i=0; i<n; ++i)
    if(buf[5*i+4]>0.0)
    {
        const double w = 1.0/buf[5*i+4];
        cuf_[i] = w*Eigen::Vector3d(buf[5*i+0],buf[5*i+1],buf[5*i+2]);
        ceta_[i] = w*buf[5*i+3];
    }
    
    if(ctn_<0.0)
    {
        cuf0_ = cuf_;
        ctn_ = p->simtime;
    }
}

void net_membrane::collar_step_end(lexer *p)
{
    // fluid acceleration at the collar nodes from the velocity change over the time step
    if(cn_.empty())
    return;
    
    if(p->simtime > ctn_ + 1.0e-12)
    {
        const double dt = p->simtime - ctn_;
        
        for(size_t i=0; i<cn_.size(); ++i)
        caf_[i] = (cuf_[i] - cuf0_[i])/dt;
    }
    
    cuf0_ = cuf_;
    ctn_ = p->simtime;
}

void net_membrane::collar_assemble(lexer *p, double h, const vector<int> &gid,
                                   const function<void(int,int,const Eigen::Matrix3d&)> &addblock, Eigen::VectorXd &rhs)
{
    // collar terms of one backward-Euler sub-step from x_, vs_: block (m_a + h C + h^2 K), right-hand side
    // h F(x_, vs_) + h C vs_ + m_a vs_ (the node mass and fabric terms are assembled by structure_solve)
    if(cn_.empty())
    return;
    
    const int n = cn_.size();
    const double rho = p->W1;
    const double g = fabs(p->W22);
    const double r = 0.5*prm.cD;
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const Eigen::Vector3d ez(0.0,0.0,1.0);
    
    // pipe: weight, buoyancy, inertia, drag
    for(int i=0; i<n; ++i)
    {
        const int q = cn_[i];
        const int gq = gid[q];
        
        if(gq<0)
        continue;
        
        const double l = cl_[i];
        const Eigen::Vector3d &x = x_[q];
        Eigen::Vector3d t = x_[cn_[(i+1)%n]] - x_[cn_[(i+n-1)%n]];
        t /= MAX(t.norm(),1.0e-20);
        const Eigen::Matrix3d Pn = I - t*t.transpose();
        
        // submerged part of the section
        const double s = MAX(0.0, MIN(prm.cD, ceta_[i] - (x(2) - r)));
        double As=0.0, b=0.0;
        
        if(s>=prm.cD)
        As = PI*r*r;
        else if(s>0.0)
        {
            const double th = 2.0*acos((r-s)/r);
            As = 0.5*r*r*(th - sin(th));
            b = 2.0*sqrt(MAX(0.0, r*r - (r-s)*(r-s)));
        }
        
        Eigen::Vector3d F = (rho*As - prm.cm)*g*l*ez;
        
        // added mass, Froude-Krylov and fluid inertia
        const Eigen::Matrix3d Ma = rho*prm.cCa*As*l*Pn;
        F += rho*(1.0 + prm.cCa)*As*l*(Pn*caf_[i]);
        
        // drag, linearised about the current velocity
        const Eigen::Vector3d ur = Pn*(cuf_[i] - vs_[q]);
        const double um = ur.norm();
        const double cd = 0.5*rho*prm.cCd*s*l;
        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        
        if(um>1.0e-12)
        {
            const Eigen::Vector3d e = ur/um;
            F += cd*um*ur;
            C = cd*um*(Pn + e*e.transpose());
        }
        
        // hydrostatic stiffness
        const Eigen::Matrix3d K = rho*g*b*l*(ez*ez.transpose());
        
        addblock(gq,gq,Ma + h*C + h*h*K);
        rhs.segment<3>(3*gq) += h*F + (Ma + h*C)*vs_[q];
    }
    
    // axial springs between neighbouring ring nodes
    cNmax_ = 0.0;
    
    for(int i=0; i<n; ++i)
    {
        const int a = cn_[i], bb = cn_[(i+1)%n];
        const int ga = gid[a], gb = gid[bb];
        const Eigen::Vector3d dx = x_[bb] - x_[a];
        const double L = dx.norm();
        
        if(L<1.0e-20)
        continue;
        
        const Eigen::Vector3d e = dx/L;
        const Eigen::Matrix3d P = e*e.transpose();
        const double k = prm.cEA/cL0_[i];
        const Eigen::Vector3d fa = k*(L - cL0_[i])*e;
        const Eigen::Matrix3d B = h*h*k*(P + MAX(0.0, 1.0 - cL0_[i]/L)*(I - P));
        
        cNmax_ = MAX(cNmax_, fabs(k*(L - cL0_[i])));
        
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
    }
    
    // corotational bending: best-fit rigid motion of the ring (Kabsch), then F = -k_b D^T D (x - c - R x_ref)
    if(prm.cEI>0.0)
    {
        Eigen::Vector3d c = Eigen::Vector3d::Zero();
        for(int i=0; i<n; ++i)
        c += x_[cn_[i]];
        c /= double(n);
        
        Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
        for(int i=0; i<n; ++i)
        H += cxr_[i]*(x_[cn_[i]] - c).transpose();
        
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d R = svd.matrixV()*svd.matrixU().transpose();
        
        if(R.determinant()<0.0)
        {
            Eigen::Matrix3d V = svd.matrixV();
            V.col(2) *= -1.0;
            R = V*svd.matrixU().transpose();
        }
        
        vector<Eigen::Vector3d> u(n), Du(n);
        for(int i=0; i<n; ++i)
        u[i] = x_[cn_[i]] - c - R*cxr_[i];
        
        for(int i=0; i<n; ++i)
        Du[i] = u[(i+n-1)%n] - 2.0*u[i] + u[(i+1)%n];
        
        double lm = 0.0;
        for(int i=0; i<n; ++i)
        lm += cl_[i];
        lm /= double(n);
        
        const double kb = prm.cEI/(lm*lm*lm);
        const double w[5] = {1.0, -4.0, 6.0, -4.0, 1.0};
        
        cMmax_ = 0.0;
        
        for(int i=0; i<n; ++i)
        {
            const int gi = gid[cn_[i]];
            
            cMmax_ = MAX(cMmax_, prm.cEI*Du[i].norm()/(lm*lm));
            
            if(gi<0)
            continue;
            
            // force: -k_b (D^T D u)_i, D symmetric on the ring
            const Eigen::Vector3d DDu = Du[(i+n-1)%n] - 2.0*Du[i] + Du[(i+1)%n];
            rhs.segment<3>(3*gi) -= h*kb*DDu;
            
            for(int o=-2; o<=2; ++o)
            {
                const int gj = gid[cn_[(i+o+n)%n]];
                
                if(gj>=0)
                addblock(gi,gj,h*h*kb*w[o+2]*I);
            }
        }
    }
    
    // mooring springs
    for(size_t m=0; m<prm.moor.size(); ++m)
    {
        const int q = cn_[mfl_[m]];
        const int gq = gid[q];
        const Eigen::Vector3d dx = Eigen::Vector3d(prm.moor[m][0],prm.moor[m][1],prm.moor[m][2]) - x_[q];
        const double L = dx.norm();
        
        if(L<1.0e-20 || gq<0)
        continue;
        
        const Eigen::Vector3d e = dx/L;
        const double T = MAX(0.0, prm.moor[m][4] + prm.moor[m][3]*(L - mL0_[m]));
        const Eigen::Matrix3d K = prm.moor[m][3]*(e*e.transpose()) + (T/L)*(I - e*e.transpose());
        
        mT_[m] = T;
        addblock(gq,gq,h*h*K);
        rhs.segment<3>(3*gq) += h*T*e;
    }
}
