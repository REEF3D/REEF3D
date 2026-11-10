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

#include"6DOF_obj.h"
#include"6DOF_deform.h"
#include"6DOF_output_dir.h"
#include"lexer.h"
#include"ghostcell.h"
#include"reefckpt.h"
#include"6DOF_brain.h"
#include"6DOF_brain_link.h"
#include<cstdio>

//  Prescribed deformation of a 6DOF body (X 195): undulatory swimming. The deformation acts on the
//  body frame hull (sixdof_geometry::transform), its velocity is added to the rigid-body velocity of
//  the direct forcing and of the wall shear stress, the loads and the rigid-body equations are
//  those of the 6DOF body (self-propulsion: X 11 with free surge, sway and yaw, 2D surge, heave
//  and pitch).

void sixdof_obj::deform_setup(lexer *p, ghostcell *pgc)
{
    if(fish_==nullptr)
    return;

    fish_->setup(geom);
    
    // brain (X 193): before the first time of the gait (the step of the recoil quadrature)
    brain_setup(p);
    
    deform_time(p->simtime);
    geom.deform = fish_;
    
    // X 196 2 with X 198: the discrete recoil correction of the grid (zero at the start)
    if(p->X198>0 && p->X196==2)
    {
    fish_->set_corr(sixdof_deform::correction());
    dc_on_ = true;
    }
    
    // the hull in the deformed shape of the start time (with a ramp: the undeformed one)
    geom.transform(R_,c_);

    if(p->mpirank==0)
    {
        cout<<"6DOF body "<<n6DOF<<": deforming (X 195), L: "<<fish_->length()<<" m, head at x_B = "<<fish_->xhead()
            <<", tail at x_B = "<<fish_->xtail()<<", volume: "<<fish_->volume()<<" m^3 (body: "<<geom.volume()<<")"
            <<", tail amplitude: "<<fish_->tail_amplitude()<<" m, f: "<<fish_->frequency()<<" Hz, lateral axis: "
            <<(fish_->lateral()==1 ? "y" : "z")<<", recoil: "<<(p->X196==0 ? "off" : (dc_on_ ? "on, discrete correction of the grid (X 196 2)" : "on"))<<endl;
        
        if(fish_->gait_on())
        {
            cout<<"6DOF body "<<n6DOF<<": gait control (X 199):";
            
            for(int n=0; n<p->X199; ++n)
            if(p->X199_id[n]==n6DOF)
            cout<<" t "<<p->X199_t[n]<<": f "<<p->X199_f[n]<<" a "<<p->X199_a[n]<<" b "<<p->X199_b[n]<<" t_s "<<p->X199_ts[n]<<";";
            
            cout<<endl;
        }

        char str[1000];
        snprintf(str,sizeof(str),"%s/REEF3D_6DOF_fish_%i.dat",sixdof_output_dir(p),n6DOF);

        if(reefckpt_open(printfish,str,p->I45))
        printfish<<"time \t U_swim [m/s] \t U_lat [m/s] \t heading [deg] \t tail [m] \t ramp \t P_def [W] \t St \t Re"
                 <<(fish_->gait_on() ? " \t f [Hz] \t a \t b" : "")<<endl;
    }
}

void sixdof_obj::deform_time(double t)
{
    if(fish_!=nullptr)
    fish_->set_time(t);
}

void sixdof_obj::deform_stage(lexer *p, int iter)
{
    if(fish_==nullptr)
    return;

    // the body after stage iter is at t^n + c dt: TVD RK3 (1, 1/2, 1), RK2 (1, 1), low-storage RK3 of
    // Wray (gamma != 0: 8/15, 2/3, 1)
    static const double crk3[3]  = {1.0, 0.5, 1.0};
    static const double crkls[3] = {8.0/15.0, 2.0/3.0, 1.0};

    double c = 1.0;

    if(iter>=0 && iter<3)
    {
        if(gamma[0]!=0.0)
        c = crkls[iter];

        else if(alpha[1]!=0.5)
        c = crk3[iter];
    }

    deform_time(p->simtime + c*p->dt);
    
    // X 198: the discrete recoil correction with the RK stage of the body (its rates from the
    // forcing of the last stage): y^s = c1 y^n + c2 (y^(s-1) + dt y'), low-storage RK3
    // y^s = y^(s-1) + gamma dt y' + zeta dt y'_k
    if(dc_on_ && iter>=0 && iter<3)
    {
        sixdof_deform::correction cr = fish_->corr();
        double y[3] = {cr.tx, cr.tl, cr.psi};
        const double r[3] = {cr.vx, cr.vl, cr.psit};
        const double dt = p->dt;
        
        if(gamma[0]!=0.0)
        {
            for(int q=0; q<3; ++q)
            {
            y[q] += gamma[iter]*dt*r[q] + zeta[iter]*dt*dcr_[q];
            dcr_[q] = r[q];
            }
        }
        
        else
        {
            if(iter==0)
            for(int q=0; q<3; ++q)
            dck_[q] = y[q];
            
            const double c2 = alpha[iter];
            const double c1 = (iter==0) ? 0.0 : 1.0 - c2;
            
            for(int q=0; q<3; ++q)
            y[q] = c1*dck_[q] + c2*(y[q] + dt*r[q]);
        }
        
        cr.tx = y[0];
        cr.tl = y[1];
        cr.psi = y[2];
        fish_->set_corr(cr);
    }
}

void sixdof_obj::deform_step_end()
{
    if(!dc_on_)
    return;
    
    const sixdof_deform::correction &cr = fish_->corr();
    dcn_[0] = cr.tx; dcn_[1] = cr.tl; dcn_[2] = cr.psi;
    dcn_[3] = cr.vx; dcn_[4] = cr.vl; dcn_[5] = cr.psit;
}

Eigen::Vector3d sixdof_obj::deform_vel(double x, double y, double z) const
{
    Eigen::Vector3d v = Eigen::Vector3d::Zero();

    if(fish_==nullptr)
    return v;

    // the point in the body frame: axial and lateral coordinate
    const int l = fish_->lateral();
    const double xb = R_(0,0)*(x - c_(0)) + R_(1,0)*(y - c_(1)) + R_(2,0)*(z - c_(2));
    const double lb = R_(0,l)*(x - c_(0)) + R_(1,l)*(y - c_(1)) + R_(2,l)*(z - c_(2));

    double vx, vl;
    fish_->velocity_at(xb,lb,vx,vl);

    // R (e_x vx + e_l vl)
    v(0) = R_(0,0)*vx + R_(0,l)*vl;
    v(1) = R_(1,0)*vx + R_(1,l)*vl;
    v(2) = R_(2,0)*vx + R_(2,l)*vl;

    return v;
}

void sixdof_obj::print_fish(lexer *p, double t)
{
    if(fish_==nullptr || p->mpirank!=0 || p->count%p->X19!=0 || !printfish.is_open())
    return;

    // velocity of the CoG in the body frame: swimming speed along the body axis and lateral
    const Eigen::Vector3d U(p->ufbi, p->vfbi, p->wfbi);
    const Eigen::Vector3d Ub = R_.transpose()*U;
    const int l = fish_->lateral();

    // heading: yaw (3D) or pitch (2D) of the body axis
    const double heading = (l==1) ? atan2(R_(1,0),R_(0,0)) : atan2(-R_(2,0),R_(0,0));

    const double Us = Ub(0);
    const double St = (fabs(Us)>1.0e-12) ? 2.0*fish_->frequency()*fish_->tail_amplitude()/fabs(Us) : 0.0;
    const double Re = fabs(Us)*fish_->length()/p->W2;

    double dxt, dlt;
    fish_->displacement(fish_->xtail(),0.0,dxt,dlt);

    printfish<<t<<" \t "<<Us<<" \t "<<Ub(l)<<" \t "<<heading*(180.0/PI)<<" \t "<<dlt
             <<" \t "<<fish_->ramp()<<" \t "<<fish_power_<<" \t "<<St<<" \t "<<Re;
    
    if(fish_->gait_on())
    printfish<<" \t "<<fish_->frequency()<<" \t "<<fish_->amp_factor()<<" \t "<<fish_->bias();
    
    printfish<<endl;
}

// ---------------------------------------------------------------------------------------------
// brain (X 193 id U_des b_max th_wall th_esc, targets X 192 id t x y z): the senses of each step, the
// decisions once per half tail beat as gait targets of the deforming body
// ---------------------------------------------------------------------------------------------

void sixdof_obj::brain_setup(lexer *p)
{
    int nb = -1;
    
    for(int n=0; n<p->X193; ++n)
    if(p->X193_id[n]==n6DOF)
    nb = n;
    
    if(nb<0 || fish_==nullptr)
    return;
    
    std::vector<sixdof_brain::target_point> tp;
    
    for(int n=0; n<p->X192; ++n)
    if(p->X192_id[n]==n6DOF)
    tp.push_back({p->X192_t[n], p->X192_x[n], p->X192_y[n], p->X192_z[n]});
    
    double f0 = fish_->frequency(), tramp = 0.0;
    
    for(int n=0; n<p->X195; ++n)
    if(p->X195_id[n]==n6DOF)
    {
    f0 = p->X195_f[n];
    tramp = p->X195_tramp[n];
    }
    
    delete brain_;
    brain_ = new sixdof_brain(p->X193_U[nb], p->X193_b[nb], p->X193_D[nb], p->X193_the[nb], f0, fish_->length(), p->W1,
                              tramp + 1.0/f0, fish_->lateral(), tp);
    
    // the largest frequency of the brain fixes the step of the recoil quadrature
    fish_->set_fcap(brain_->fmax());
    fish_->set_gait(gait_in_);
    gait_brain_.clear();
    
    // an external brain (X 191 id port): rank 0 connects to it, sends the hello line
    int port = 0;
    for(int n=0; n<p->X191; ++n)
    if(p->X191_id[n]==n6DOF)
    port = p->X191_port[n];
    
    if(port>0 && p->mpirank==0)
    {
        delete brain_link_;
        brain_link_ = new sixdof_brain_link();
        std::string msg, in;
        
        if(!brain_link_->connect_to(port, 60.0, msg))
        {
            cout<<"6DOF body "<<n6DOF<<": "<<msg<<endl;
            exit(1);
        }
        
        char hello[512];
        snprintf(hello,sizeof(hello),"{\"hello\":1,\"body\":%d,\"L\":%.17g,\"f0\":%.17g,\"rho\":%.17g,\"U_des\":%.17g,\"b_max\":%.17g,\"D\":%.17g,\"th_esc\":%.17g,\"t_start\":%.17g,\"twoD\":%d}",
                 n6DOF, fish_->length(), f0, p->W1, p->X193_U[nb], p->X193_b[nb], p->X193_D[nb], p->X193_the[nb], tramp + 1.0/f0, fish_->lateral()==2 ? 1 : 0);
        
        if(!brain_link_->exchange(hello, in, 60.0))
        {
            cout<<"6DOF body "<<n6DOF<<": the external brain (X 191) did not answer"<<endl;
            exit(1);
        }
        
        cout<<"6DOF body "<<n6DOF<<": external brain (X 191) on port "<<port<<": "<<in<<endl;
    }
    
    if(p->mpirank==0)
    {
        cout<<"6DOF body "<<n6DOF<<": brain (X 193), U_des "<<p->X193_U[nb]<<" m/s, b_max "<<p->X193_b[nb]
            <<", vision "<<p->X193_D[nb]<<" L, escape at "<<p->X193_the[nb]<<" q, "<<tp.size()<<" target points (X 192)"<<endl;
        
        char str[1000];
        snprintf(str,sizeof(str),"%s/REEF3D_6DOF_brain_%i.dat",sixdof_output_dir(p),n6DOF);
        
        if(reefckpt_open(printbrain,str,p->I45))
        printbrain<<"time \t target_x \t target_y \t target_z \t dist [m] \t bearing [deg] \t U [m/s] \t f [Hz] \t a \t b"
                  <<" \t dp_head/q \t p_head/q \t d_vision/L \t mode (0 target, 1 avoid, 2 escape, 3 at target)"
                  <<" \t stride [L] \t current_x \t current_y \t current_z [m/s] \t U_water [m/s]"<<endl;
    }
}

bool sixdof_obj::brain_due(lexer *p) const
{
    return brain_!=nullptr && brain_->due(p->simtime + p->dt);
}

void sixdof_obj::brain_update(lexer *p, ghostcell *pgc, const std::vector<double> *pr, const std::vector<double> *vis)
{
    if(brain_==nullptr)
    return;
    
    // the senses of this step: the body at t^(n+1)
    const int la = fish_->lateral();
    const Eigen::Vector3d ax = R_.col(0), lat = R_.col(la);
    const Eigen::Vector3d U(p->ufbi, p->vfbi, p->wfbi);
    brain_->sense(p->dt, ax, lat, c_, U, pr, ll_s_, ll_side_);
    
    // a decision for the next step
    const double tn = p->simtime + p->dt;
    
    if(!brain_->due(tn))
    return;
    
    fish_->set_time(tn);
    sixdof_brain::decision d = brain_->decide(tn, fish_->frequency(), fish_->amp_factor(), fish_->bias(), vis);
    
    // an external brain (X 191): the senses and the decision of the own brain out, its decision back
    // (rank 0), to all ranks by a sum with zeros elsewhere (exact)
    if(p->X191>0 && d.on)
    {
        double v[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        
        if(p->mpirank==0 && brain_link_!=nullptr)
        {
            const sixdof_brain::perception &q = brain_->perceived();
            std::string s;
            char b[256];
            auto num = [&](const char *k, double x) {snprintf(b,sizeof(b),"\"%s\":%.17g,",k,x); s += b;};
            auto vec = [&](const char *k, const double *x, size_t n)
            {
                s += std::string("\"") + k + "\":[";
                for(size_t m=0; m<n; ++m)
                {
                    snprintf(b,sizeof(b),"%.17g%s",x[m],m+1<n ? "," : "");
                    s += b;
                }
                s += "],";
            };
            s = "{";
            num("t",tn); num("f",fish_->frequency()); num("a",fish_->amp_factor()); num("b",fish_->bias());
            vec("pos",q.c.data(),3); vec("axis",q.ex.data(),3); vec("lat",q.el.data(),3);
            num("U",q.U); num("stride",q.stride);
            if(q.target)
            vec("target",q.T.data(),3);
            vec("rays",q.rays.data(),q.rays.size());
            vec("ray_angles",sixdof_brain::rays().data(),sixdof_brain::rays().size());
            vec("ll_p",q.llp.data(),q.llp.size());
            vec("ll_s",brain_->ll_s().data(),brain_->ll_s().size());
            std::vector<double> sd(brain_->ll_side().begin(), brain_->ll_side().end());
            vec("ll_side",sd.data(),sd.size());
            const sixdof_brain::record &r = brain_->last();
            snprintf(b,sizeof(b),"\"own\":{\"f\":%.17g,\"a\":%.17g,\"b\":%.17g,\"ts\":%.17g,\"mode\":%d,\"second\":%d,\"t2\":%.17g,\"f2\":%.17g,\"a2\":%.17g,\"b2\":%.17g,\"ts2\":%.17g}}",
                     d.f,d.a,d.b,d.ts,r.mode,d.second ? 1 : 0,d.second ? d.t2 : 0.0,d.second ? d.f2 : 0.0,d.second ? d.a2 : 0.0,d.second ? d.b2 : 0.0,d.second ? d.ts2 : 0.0);
            s += b;
            
            std::string in;
            if(!brain_link_->exchange(s, in, 600.0))
            {
                cout<<"6DOF body "<<n6DOF<<": the external brain (X 191) did not answer at t = "<<tn<<endl;
                exit(1);
            }
            
            typedef sixdof_brain_link L;
            v[0] = L::number(in,"f",d.f);
            v[1] = L::number(in,"a",d.a);
            v[2] = L::number(in,"b",d.b);
            v[3] = L::number(in,"ts",0.0);
            v[4] = L::has(in,"t2") ? 1.0 : 0.0;
            v[5] = L::number(in,"t2",0.0);
            v[6] = L::number(in,"f2",0.0);
            v[7] = L::number(in,"a2",1.0);
            v[8] = L::number(in,"b2",0.0);
            v[9] = L::number(in,"ts2",0.0);
            v[10] = L::number(in,"next",0.0);
            v[11] = 1.0;
        }
        
        for(int q=0; q<12; ++q)
        v[q] = pgc->globalsum(v[q]);
        
        if(v[11]!=1.0)
        {
            cout<<"6DOF body "<<n6DOF<<": no external brain (X 191) on rank 0"<<endl;
            exit(1);
        }
        
        // the gait stays within what the recoil quadrature allows (2 f0)
        const double fm = brain_->fmax();
        d.f = std::max(0.0, std::min(fm, v[0])); d.a = v[1]; d.b = v[2]; d.ts = v[3];
        d.second = (v[4]!=0.0);
        if(d.second)
        {
            d.t2 = std::max(tn, v[5]); d.f2 = std::max(0.0, std::min(fm, v[6])); d.a2 = v[7]; d.b2 = v[8]; d.ts2 = v[9];
        }
        
        const double fn = d.second ? d.f2 : d.f;
        brain_->external(d, tn + (v[10]>0.0 ? v[10] : std::min(0.5/std::max(fn,1.0e-3), 0.5)));
    }
    
    if(d.on)
    {
        gait_brain_.push_back({d.t, d.f, d.a, d.b, d.ts});
        fish_->push_gait(gait_brain_.back());
        
        if(d.second)
        {
            gait_brain_.push_back({d.t2, d.f2, d.a2, d.b2, d.ts2});
            fish_->push_gait(gait_brain_.back());
        }
    }
    
    if(p->mpirank==0 && printbrain.is_open())
    {
        const sixdof_brain::record &r = brain_->last();
        printbrain<<r.t<<" \t "<<r.tx<<" \t "<<r.ty<<" \t "<<r.tz<<" \t "<<r.dist<<" \t "<<r.psi*(180.0/PI)<<" \t "<<r.U
                  <<" \t "<<r.Ff<<" \t "<<r.Fa<<" \t "<<r.Fb<<" \t "<<r.dpl<<" \t "<<r.phd<<" \t "<<r.dvis<<" \t "<<r.mode
                  <<" \t "<<r.s<<" \t "<<r.cx<<" \t "<<r.cy<<" \t "<<r.cz<<" \t "<<r.Uw<<endl;
    }
}
