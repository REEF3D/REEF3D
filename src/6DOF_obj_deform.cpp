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
