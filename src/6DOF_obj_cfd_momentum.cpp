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

#include"6DOF_obj_cfd.h"
#include"6DOF_deform.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<vector>

//  X 198: fluid-momentum coupling of a CFD body (rigidity constraint of the fictitious domain,
//  Patankar et al. 2000, Sharma & Patankar 2005, Gazzola et al. 2011). The direct forcing makes the
//  fluid in the hull move with the body; the predicted velocity of the next stage u* (convection,
//  diffusion and the pressure of the last projection on all of the domain) carries in the hull
//  region the momentum the fluid has given the body in that stage. The rigid-body motion fitted to
//  u* there (least squares with the weights rho (2H - 1) dV over the faces inside the hull, H the
//  forcing weight, 1/2 at the hull: the band outside it is fluid and left out, with it the swimmer of
//  validation 03 swam at a fifth of its speed; the deformation velocity of X 195 at the times of the
//  RK combination removed) is the velocity of the body after the stage when it is neutrally buoyant;
//  the loads of the stage follow from it with the RK stage of the body:
//      p^s = c1 p^n + c2 (p^(s-1) + dt F),   F = rho_f V (U_fit - (c1 p^n + c2 p^(s-1))/M)/(c2 dt)
//                                                  + (M - rho_f V) g
//  (rotations alike with the angular momentum in the body frame). The fluid in the hull region is
//  part of the solution of the projection with its added mass, so the coupling needs no added-mass
//  stabilisation (X 197): a neutrally buoyant or heavier body. For a lighter body (rho_b < rho_f)
//  the fluid momentum enters with the factor rho_f/rho_b > 1, which is explicit again; X 198 1 c
//  relaxes the fitted velocity with the factor c.

template<class G>
void sixdof_obj_cfd::momentum_acc(lexer *p, fdm *a, G val, Eigen::Matrix<double,6,6> &N, Eigen::Matrix<double,6,1> &b)
{
    // the normal equations of the least-squares rigid-body motion (U, omega about the CoG) of the face
    // values val(d,x,y,z) in the hull region of the grid (p, a): weights rho (2H - 1) dV at the faces
    // inside the hull (H the forcing weight, 1/2 at the hull), rows g = (e_d, r x e_d); this rank's
    // part, added to N, b
    auto add = [&](int d, double x, double y, double z, double w, double u)
    {
        const Eigen::Vector3d r(x - c_(0), y - c_(1), z - c_(2));
        Eigen::Vector3d e = Eigen::Vector3d::Zero();
        e(d) = 1.0;
        Eigen::Matrix<double,6,1> g;
        g.head<3>() = e;
        g.tail<3>() = r.cross(e);
        N += w*g*g.transpose();
        b += w*u*g;
    };
    
    ULOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;
        
        const double H = 2.0*Hsolidface(p,a,1,0,0) - 1.0;
        
        if(H<=0.0)
        continue;
        
        const double ro = 0.5*(a->ro(i,j,k) + a->ro(i+1,j,k));
        add(0,p->pos1_x(),p->pos1_y(),p->pos1_z(),ro*H*p->DXP[IP]*p->DYN[JP]*p->DZN[KP],val(0,p->pos1_x(),p->pos1_y(),p->pos1_z()));
    }
    
    if(p->j_dir==1)
    VLOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;
        
        const double H = 2.0*Hsolidface(p,a,0,1,0) - 1.0;
        
        if(H<=0.0)
        continue;
        
        const double ro = 0.5*(a->ro(i,j,k) + a->ro(i,j+1,k));
        add(1,p->pos2_x(),p->pos2_y(),p->pos2_z(),ro*H*p->DXN[IP]*p->DYP[JP]*p->DZN[KP],val(1,p->pos2_x(),p->pos2_y(),p->pos2_z()));
    }
    
    WLOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;
        
        const double H = 2.0*Hsolidface(p,a,0,0,1) - 1.0;
        
        if(H<=0.0)
        continue;
        
        const double ro = 0.5*(a->ro(i,j,k) + a->ro(i,j,k+1));
        add(2,p->pos3_x(),p->pos3_y(),p->pos3_z(),ro*H*p->DXN[IP]*p->DYN[JP]*p->DZP[KP],val(2,p->pos3_x(),p->pos3_y(),p->pos3_z()));
    }
}

Eigen::Matrix<double,6,1> sixdof_obj_cfd::momentum_solve(lexer *p, ghostcell *pgc, Eigen::Matrix<double,6,6> &N, Eigen::Matrix<double,6,1> &b, bool &ok)
{
    // the normal equations summed over the ranks, solved for the modes of the grid (2D: surge, heave,
    // pitch)
    for(int q=0; q<6; ++q)
    {
        b(q) = pgc->globalsum(b(q));
        
        for(int r=0; r<6; ++r)
        N(q,r) = pgc->globalsum(N(q,r));
    }
    
    int idx[6], n=0;
    for(int q=0; q<6; ++q)
    if(p->j_dir==1 || q==0 || q==2 || q==4)
    idx[n++] = q;
    
    Eigen::MatrixXd Ns(n,n);
    Eigen::VectorXd bs(n);
    
    for(int q=0; q<n; ++q)
    {
        bs(q) = b(idx[q]);
        
        for(int r=0; r<n; ++r)
        Ns(q,r) = N(idx[q],idx[r]);
    }
    
    Eigen::Matrix<double,6,1> res = Eigen::Matrix<double,6,1>::Zero();
    ok = false;
    
    if(Ns(0,0)>0.0)
    {
        Eigen::LDLT<Eigen::MatrixXd> ldlt(Ns);
        
        if(ldlt.info()==Eigen::Success)
        {
            const Eigen::VectorXd qs = ldlt.solve(bs);
            
            for(int q=0; q<n; ++q)
            res(idx[q]) = qs(q);
            
            ok = true;
        }
    }
    
    return res;
}

void sixdof_obj_cfd::momentum_fit(lexer *p, fdm *a, ghostcell *pgc, field &uvel, field &vvel, field &wvel, int iter)
{
    momentum_fit_begin(p,iter);
    momentum_fit_add(p,a,uvel,vvel,wvel);
    momentum_fit_end(p,pgc);
}

void sixdof_obj_cfd::momentum_fit_begin(lexer *p, int iter)
{
    // the deformation velocity in u* of the stage: the RK combination c1 u_def(t^n) + c2 u_def(t^(s-1))
    // of the forced deformation velocities (TVD RK3, RK2; low-storage RK3 and the first stage:
    // u_def(t^(s-1)))
    const bool ls = (gamma[0]!=0.0);
    mf_c2_ = (ls || iter<=0 || iter>2) ? 1.0 : alpha[iter];
    mf_c1_ = 1.0 - mf_c2_;
    mf_tn_ = p->simtime;
    mf_N_.setZero();
    mf_b_.setZero();
    
    // the first fit of a run from an initial velocity of the body (X 102, X 103): the fluid in the
    // hull has not been forced yet, the body keeps its velocity
    mf_start_ = (mf_first_ && (p->X102>0 || p->X103==1));
    mf_first_ = false;
}

void sixdof_obj_cfd::momentum_fit_add(lexer *p, fdm *a, field &uvel, field &vvel, field &wvel)
{
    // the part of the grid (p, a) with its predicted velocity (level 0, or a patch of the mesh
    // refinement)
    const double c1 = mf_c1_, c2 = mf_c2_;
    std::vector<double> udn[3];
    
    if(fish_!=nullptr && c1>0.0)
    {
        const double tcur = fish_->time();
        const sixdof_deform::correction ccur = fish_->corr();
        deform_time(mf_tn_);
        
        if(dc_on_)
        {
        sixdof_deform::correction cn;
        cn.tx = dcn_[0]; cn.tl = dcn_[1]; cn.psi = dcn_[2];
        cn.vx = dcn_[3]; cn.vl = dcn_[4]; cn.psit = dcn_[5];
        fish_->set_corr(cn);
        }
        
        ULOOP
        if(2.0*Hsolidface(p,a,1,0,0) - 1.0 > 0.0 && (amr_cover==nullptr || !amr_cover[IJK]))
        udn[0].push_back(deform_vel(p->pos1_x(),p->pos1_y(),p->pos1_z())(0));
        
        if(p->j_dir==1)
        VLOOP
        if(2.0*Hsolidface(p,a,0,1,0) - 1.0 > 0.0 && (amr_cover==nullptr || !amr_cover[IJK]))
        udn[1].push_back(deform_vel(p->pos2_x(),p->pos2_y(),p->pos2_z())(1));
        
        WLOOP
        if(2.0*Hsolidface(p,a,0,0,1) - 1.0 > 0.0 && (amr_cover==nullptr || !amr_cover[IJK]))
        udn[2].push_back(deform_vel(p->pos3_x(),p->pos3_y(),p->pos3_z())(2));
        
        deform_time(tcur);
        
        if(dc_on_)
        fish_->set_corr(ccur);
    }
    
    size_t kn[3] = {0,0,0};
    
    auto rel = [&](int d, double x, double y, double z)
    {
        double u = (d==0) ? uvel(i,j,k) : ((d==1) ? vvel(i,j,k) : wvel(i,j,k));
        
        if(fish_!=nullptr)
        {
            const double v = deform_vel(x,y,z)(d);
            
            if(c1>0.0 && kn[d]<udn[d].size())
            u -= c1*udn[d][kn[d]++] + c2*v;
            
            else
            u -= v;
        }
        
        return u;
    };
    
    momentum_acc(p,a,rel,mf_N_,mf_b_);
}

void sixdof_obj_cfd::momentum_fit_end(lexer *p, ghostcell *pgc)
{
    bool ok;
    mom_fit_ = momentum_solve(p,pgc,mf_N_,mf_b_,ok);
    mom_ok_ = ok;
    
    if(mf_start_)
    {
    mom_fit_ = u_fb;
    mom_ok_ = true;
    }
}

void sixdof_obj_cfd::momentum_deform(lexer *p, fdm *a, ghostcell *pgc)
{
    momentum_deform_begin();
    momentum_deform_add(p,a);
    momentum_deform_end(p,pgc);
}

void sixdof_obj_cfd::momentum_deform_begin()
{
    // the rates of the discrete recoil correction at the pose of the stage: minus the rigid part of
    // the deformation velocity (without the rates) in the hull region, in the plane of the
    // undulation (body frame: axial, lateral, rotation about e_x x e_l)
    mf_N_.setZero();
    mf_b_.setZero();
    
    if(!dc_on_)
    return;
    
    sixdof_deform::correction cr = fish_->corr();
    cr.vx = cr.vl = cr.psit = 0.0;
    fish_->set_corr(cr);
}

void sixdof_obj_cfd::momentum_deform_add(lexer *p, fdm *a)
{
    if(!dc_on_)
    return;
    
    momentum_acc(p,a,[&](int d, double x, double y, double z){return deform_vel(x,y,z)(d);},mf_N_,mf_b_);
}

void sixdof_obj_cfd::momentum_deform_end(lexer *p, ghostcell *pgc)
{
    if(!dc_on_)
    return;
    
    bool ok;
    const Eigen::Matrix<double,6,1> q = momentum_solve(p,pgc,mf_N_,mf_b_,ok);
    
    if(!ok)
    return;
    
    sixdof_deform::correction cr = fish_->corr();
    const int l = fish_->lateral();
    const Eigen::Vector3d QB = R_.transpose()*q.head<3>();
    const Eigen::Vector3d WB = R_.transpose()*q.tail<3>();
    Eigen::Vector3d ex = Eigen::Vector3d::Zero(), el = Eigen::Vector3d::Zero();
    ex(0) = 1.0;
    el(l) = 1.0;
    
    cr.vx = -QB(0);
    cr.vl = -QB(l);
    cr.psit = -WB.dot(ex.cross(el));
    fish_->set_corr(cr);
}

void sixdof_obj_cfd::momentum_loads(lexer *p, int iter)
{
    // the loads of the stage from the fitted velocity, with the RK stage of the body:
    // TVD RK3 and RK2: y^s = c1 y^n + c2 (y^(s-1) + dt dy), c2 = alpha (iter 0: c1 = 0);
    // low-storage RK3: y^s = y^(s-1) + gamma dt dy + zeta dt dy_k
    if(!mom_ok_ || iter<0 || iter>2)
    return;

    const double dt = p->dt;
    const bool ls = (gamma[0]!=0.0);
    const double c2 = ls ? gamma[iter] : alpha[iter];
    const double c1 = ls ? 0.0 : (iter==0 ? 0.0 : 1.0 - c2);

    const double rf = p->W1;
    const double rbd = (Vfb>0.0) ? Mass_fb/Vfb : rf;
    const double ratio = p->X198_c*rf/rbd;

    const Eigen::Vector3d U = mom_fit_.head<3>();
    const Eigen::Vector3d W = mom_fit_.tail<3>();
    const Eigen::Vector3d g(p->W20, p->W21, p->W22);

    // momentum without the loads of the stage
    Eigen::Vector3d p0, h0;

    if(ls)
    {
    p0 = rb.p + zeta[iter]*dt*rb.dpk;
    h0 = rb.h + zeta[iter]*dt*rb.dhk;
    }

    else
    {
    p0 = c1*rb.pk + c2*rb.p;
    h0 = c1*rb.hk + c2*rb.h;
    }

    // translation: the body at the fitted velocity (neutrally buoyant), gravity of the body instead
    // of that of the fluid in the hull
    const Eigen::Vector3d Fg = (Mass_fb - rf*Vfb)*g;
    Eigen::Vector3d F = ratio*(Mass_fb*U - p0)/(c2*dt) + Fg;
    
    // low-storage RK3: the fluid of the hull region carries the momentum of the last stage and its own
    // zeta part; the body takes the momentum the fluid gained (ratio of the masses) and its own load Fg
    // with gamma + zeta (before: p0 held zeta dt of the last load of the body, the momentum of the fit
    // of the last stage included, and Fg with gamma alone, 1.68 times Fg over the step)
    if(ls && (ratio!=1.0 || Fg.squaredNorm()>0.0))
    {
    const Eigen::Vector3d ps = rb.p + ratio*(Mass_fb*U - rb.p) + (gamma[iter] + zeta[iter])*dt*Fg;
    F = (ps - p0)/(c2*dt);
    }

    // rotation: the angular momentum of the body frame at the fitted rate; the gyroscopic part of dh
    // from the kinematics with zero moment
    const Eigen::Vector3d M0 = rb.M;
    rb.M.setZero();
    rb.derivatives_rot();
    const Eigen::Vector3d dh0 = rb.dh;
    rb.M = M0;

    const Eigen::Vector3d htarget = rb.I*(R_.transpose()*W);
    Eigen::Vector3d dhM = ratio*(htarget - h0)/(c2*dt) - dh0;
    
    if(ls && ratio!=1.0)
    dhM = (rb.h + ratio*(htarget - rb.h) - h0)/(c2*dt) - dh0;
    const Eigen::Vector3d M = R_*dhM;

    Xe = F(0);
    Ye = F(1);
    Ze = F(2);
    Ke = M(0);
    Me = M(1);
    Ne = M(2);
}
