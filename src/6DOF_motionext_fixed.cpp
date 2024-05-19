/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_motionext_fixed.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

sixdof_motionext_fixed::sixdof_motionext_fixed(lexer *p, ghostcell *pgc)
{
    ini(p,pgc);
}

sixdof_motionext_fixed::~sixdof_motionext_fixed()
{
}

void sixdof_motionext_fixed::ini(lexer *p, ghostcell *pgc)
{
    Uext = p->X210_u;
    Vext = p->X210_v;
    Wext = p->X210_w;

    Pext = p->X211_p;
    Qext = p->X211_q;
    Rext = p->X211_r;
}

void sixdof_motionext_fixed::motionext_trans(lexer *p, ghostcell *pgc, Eigen::Vector3d& dp_, Eigen::Vector3d& dc_)
{
    if (p->X11_u==2)
    {
        dp_(0) = 0.0;
        dc_(0) = Uext*ramp_vel(p);
    }

    if (p->X11_v==2)
    {
        dp_(1) = 0.0;
        dc_(1) = Vext*ramp_vel(p);
    }

    if (p->X11_w==2)
    {
        dp_(2) = 0.0;
        dc_(2) = Wext*ramp_vel(p);
    }
}

void sixdof_motionext_fixed::motionext_rot(lexer *p, Eigen::Vector3d& dh_, Eigen::Vector3d& h_, Eigen::Vector4d& de_, Eigen::Matrix<double, 3, 4>&G_,  Eigen::Matrix3d&I_)
{

    if(p->X11_p==2)
    omega_(0) = Pext*ramp_vel(p);
    
    if(p->X11_q==2)
    omega_(1) = Qext*ramp_vel(p);
    
    if(p->X11_r==2)
    omega_(2) = Rext*ramp_vel(p);
    
    h_ = I_*omega_;
    
    dh_ << 0.0,0.0,0.0;
    
    if(p->mpirank==0)
    cout<<"MOTIONEXT_ROT: "<<h_(0)<<" "<<h_(1)<<" "<<h_(2)<<endl;

    de_ = 0.5*G_.transpose()*I_.inverse()*h_;
}

double sixdof_motionext_fixed::ramp_vel(lexer *p)
{
    double f=1.0;

    if(p->X205==1 && p->X206==1 && p->simtime>=p->X206_ts && p->simtime<p->X206_te)
    f = (p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts);

    if(p->X205==2 && p->X206==1 && p->simtime>=p->X206_ts && p->simtime<p->X206_te)
    f = (p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts)-(1.0/PI)*sin(PI*(p->simtime-p->X206_ts)/(p->X206_te-p->X206_ts));

    if(p->X206==1 && p->simtime<p->X206_ts)
    f=0.0;

    return f;
}

double sixdof_motionext_fixed::ramp_draft(lexer *p)
{
    double f=1.0;

    if(p->X205==1 && p->X207==1 && p->simtime>=p->X207_ts && p->simtime<p->X207_te)
    f = p->simtime/(p->X207_te-p->X207_ts);

    if(p->X205==2 && p->X207==1 && p->simtime>=p->X207_ts && p->simtime<p->X207_te)
    f = p->simtime/(p->X207_te-p->X207_ts) - (1.0/PI)*sin(PI*(p->simtime/(p->X207_te-p->X207_ts)));

    if(p->X207==1 && p->simtime<p->X207_ts)
    f=0.0;

    return f;
}
