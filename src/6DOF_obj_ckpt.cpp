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

//  A 6DOF body in checkpoints (P 48, I 45).  The state is the same on all ranks: the rigid body with
//  its stage values and history (sixdof_rigidbody::ckpt), the loads of the last stage (the explicit
//  schemes and the fictitious-mass term use them), the body velocities, the print times, and the
//  state of the motion file, the mooring lines, the load models (ship), the power take-off (X 500) and the
//  nets (X 320).  The arrays of the body on
//  the grid (its level set, Heaviside functions, forcing) are those of the solver's checkpoint.
//  restore puts the trimesh at the pose (geom.transform) and sets the body velocities of the lexer.

#include"6DOF_obj.h"
#include"6DOF_motionext.h"
#include"mooring.h"
#include"net_interface.h"
#include"lexer.h"
#include"ghostcell.h"
#include"ckpt_pack.h"
#include"6DOF_brain.h"

template<class P>
void sixdof_obj::ckpt_fields(lexer *p, P &q)
{
    rb.ckpt(q);
    q.m(u_fb);
    q.b(ext_on_); q.m(ext_v_); q.m(ext_w_);

    for(double *x : {&Sfx_n,&Sfy_n,&Sfz_n,&SKx_n,&SKy_n,&SKz_n,
                     &Xext,&Yext,&Zext,&Kext,&Mext,&Next, &Xe,&Ye,&Ze,&Ke,&Me,&Ne, &Xd,&Yd,&Zd,&Kd,&Md,&Nd,
                     &Uext,&Vext,&Wext,&Pext,&Qext,&Rext, &curr_time,&printtime,&printtimenormal,&wm_v1,&wm_v2})
    q.d(*x);
    q.i(fictmass_ini); q.i(printnormal_count); q.i(triangle_token); q.i(timecount); q.i(timecount_old); q.i(wm_tri0);

    // the print windows (P 35), allocated by print_ini_vtp, print_ini_stl
    if(printtime_wT!=nullptr && (p->X50==1 || p->X50==2))
    q.ad(printtime_wT,p->P35);
    
    // the P 185 windows (SFLOW, NHFLOW)
    if(printtime_wF!=nullptr && p->P185>0 && (p->X50==1 || p->X50==2))
    q.ad(printtime_wF,p->P185);

    q.vd(strip_F_);
    q.m(Aadd_); q.b(am_on_);
    for(vector<double> *a : {&Xme,&Yme,&Zme,&Kme,&Mme,&Nme,&Xne,&Yne,&Zne,&Kne,&Mne,&Nne})
    q.vd(*a);
}

void sixdof_obj::ckpt_state(lexer *p, std::vector<double> &v)
{
    ckpt_put w{v};
    w.d(1.0);   // version
    ckpt_fields(p,w);

    // the parts with their own state, each with its length
    auto part = [&](const std::function<void(std::vector<double>&)> &f)
    {
        std::vector<double> s;
        f(s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    };

    part([&](std::vector<double> &s){ if(pmotion!=nullptr) pmotion->ckpt_state(s); });

    v.push_back(double(pmooring.size()));
    for(mooring *m : pmooring)
    part([&](std::vector<double> &s){ m->ckpt_state(s); });

    v.push_back(double(pload.size()));
    for(sixdof_load *l : pload)
    part([&](std::vector<double> &s){ l->ckpt_state(s); });

    // the power take-off (X 500)
    part([&](std::vector<double> &s){ pto_ckpt_state(s); });

    // the nets (X 320)
    part([&](std::vector<double> &s){ if(pnetinter!=nullptr) pnetinter->ckpt_state(s); });

    // the acceleration of the last stage of the stabilised CFD coupling (am_cfd_)
    part([&](std::vector<double> &s){ if(am_cfd_>0.0) s.assign(am_prev_.data(),am_prev_.data()+6); });
    
    // the deforming body (X 195): the time of the deformation (the checkpoint is read before the
    // time of the run), the power of the deformation, the discrete recoil correction (X 196 2) and
    // its registers of the RK stages; the momentum of the fit (X 198); the state of the coupling
    part([&](std::vector<double> &s)
    {
        s.push_back(fish_!=nullptr ? 1.0 : 0.0);
        if(fish_!=nullptr)
        {
            const sixdof_deform::correction &c = fish_->corr();
            for(double x : {fish_->time(), fish_power_, c.tx, c.tl, c.psi, c.vx, c.vl, c.psit})
            s.push_back(x);
            for(int q=0; q<6; ++q) s.push_back(dcn_[q]);
            for(int q=0; q<6; ++q) s.push_back(dck_[q]);
            for(int q=0; q<3; ++q) s.push_back(dcr_[q]);
        }
        for(int q=0; q<6; ++q) s.push_back(mom_fit_(q));
        s.push_back(mom_ok_ ? 1.0 : 0.0);
        std::vector<double> cp;
        ckpt_coupling(cp);
        s.push_back(double(cp.size()));
        s.insert(s.end(), cp.begin(), cp.end());
    });
    
    // the brain (X 193): its decisions (the gait targets it set) and its state
    part([&](std::vector<double> &s)
    {
        if(brain_==nullptr)
        return;
        s.push_back(double(gait_brain_.size()));
        for(const sixdof_deform::gait_target &g : gait_brain_)
        for(double x : {g.t, g.f, g.a, g.b, g.ts})
        s.push_back(x);
        brain_->ckpt_state(s);
    });
}

bool sixdof_obj::ckpt_restore(lexer *p, ghostcell *pgc, const std::vector<double> &v)
{
    ckpt_get r{v};
    double version=0.0;
    r.d(version);
    if(!r.ok || int(version)!=1)
    return false;

    ckpt_fields(p,r);
    if(!r.ok)
    return false;

    // a part: its length, then its values (the object reads what it knows)
    auto part = [&](const std::function<void(const std::vector<double>&, size_t&)> &f)
    {
        if(!r.has(1))
        return;
        const size_t n = size_t(v[r.k++]);
        if(!r.has(n))
        return;
        if(n>0)
        {
            const std::vector<double> s(v.begin()+r.k,v.begin()+r.k+n);
            size_t j=0;
            f(s,j);
        }
        r.k += n;
    };

    part([&](const std::vector<double> &s, size_t &j){ if(pmotion!=nullptr) pmotion->ckpt_restore(s,j); });

    int nm=0;
    r.i(nm);
    for(int q=0; q<nm && r.ok; ++q)
    part([&](const std::vector<double> &s, size_t &j){ if(q<(int)pmooring.size()) pmooring[q]->ckpt_restore(s,j); });

    int nl=0;
    r.i(nl);
    for(int q=0; q<nl && r.ok; ++q)
    part([&](const std::vector<double> &s, size_t &j){ if(q<(int)pload.size()) pload[q]->ckpt_restore(s,j); });

    // the power take-off (X 500; absent in checkpoints written before it was in them)
    bool pto_ok = true;
    if(r.ok && r.k<v.size())
    part([&](const std::vector<double> &s, size_t&){ pto_ok = pto_ckpt_restore(s); });
    else
    pto_ok = !pto_on_;
    if(!pto_ok)
    return false;

    // the nets (X 320; absent in checkpoints written before they were in them: refused with nets)
    bool net_ok = !(p->X320>0 && pnetinter!=nullptr);
    if(r.ok && r.k<v.size())
    part([&](const std::vector<double> &s, size_t&){ net_ok = (pnetinter==nullptr) ? s.empty() : pnetinter->ckpt_restore(s); });
    if(!net_ok)
    return false;

    // the acceleration of the last stage of the stabilised CFD coupling (absent in older
    // checkpoints: 0)
    if(r.ok && r.k<v.size())
    part([&](const std::vector<double> &s, size_t&){ if(s.size()==6) for(int n=0; n<6; ++n) am_prev_(n) = s[n]; });
    
    // the deforming body, the fit of X 198 and the coupling (absent in older checkpoints: the
    // deformation at the time of the run)
    double tdef = -1.0;
    if(r.ok && r.k<v.size())
    part([&](const std::vector<double> &s, size_t&)
    {
        if(s.empty())
        return;
        size_t j = 0;
        const bool fs = (s[j++]!=0.0);
        if(fs)
        {
            if(s.size()<j+8+15)
            return;
            tdef = s[j++];
            fish_power_ = s[j++];
            sixdof_deform::correction c;
            c.tx = s[j++]; c.tl = s[j++]; c.psi = s[j++]; c.vx = s[j++]; c.vl = s[j++]; c.psit = s[j++];
            for(int q=0; q<6; ++q) dcn_[q] = s[j++];
            for(int q=0; q<6; ++q) dck_[q] = s[j++];
            for(int q=0; q<3; ++q) dcr_[q] = s[j++];
            if(fish_!=nullptr && dc_on_)
            fish_->set_corr(c);
        }
        if(s.size()<j+7)
        return;
        for(int q=0; q<6; ++q) mom_fit_(q) = s[j++];
        mom_ok_ = (s[j++]!=0.0);
        if(j<s.size())
        {
            const size_t nc = size_t(s[j++]);
            if(s.size()>=j+nc)
            ckpt_coupling_restore(std::vector<double>(s.begin()+j, s.begin()+j+nc));
        }
    });
    
    // the brain (X 193): the gait rebuilt with its decisions, its state (refused if missing)
    bool brain_ok = (brain_==nullptr);
    if(r.ok && r.k<v.size())
    part([&](const std::vector<double> &s, size_t&)
    {
        if(brain_==nullptr || s.empty())
        return;
        const size_t nb = size_t(s[0]);
        if(s.size()<1+5*nb)
        return;
        gait_brain_.clear();
        fish_->set_gait(gait_in_);
        for(size_t k=0; k<nb; ++k)
        {
            gait_brain_.push_back({s[1+5*k], s[2+5*k], s[3+5*k], s[4+5*k], s[5+5*k]});
            fish_->push_gait(gait_brain_.back());
        }
        brain_ok = brain_->ckpt_restore(std::vector<double>(s.begin()+1+5*nb, s.end()));
    });
    if(!brain_ok)
    return false;

    // the pose: trimesh, body velocities and the lexer values of the body (CoG, angles, velocities);
    // a deforming body (X 195) at the time of the checkpoint
    if(fish_!=nullptr)
    deform_time(tdef>=0.0 ? tdef : p->simtime);
    
    geom.transform(R_,c_);
    update_fbvel(p,pgc);

    return r.ok;
}
