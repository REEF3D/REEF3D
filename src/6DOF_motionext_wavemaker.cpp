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

#include"6DOF_motionext_wavemaker.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

sixdof_motionext_wavemaker::sixdof_motionext_wavemaker(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF_motion  Wavemaker "<<endl;
    
    
    ini(p,pgc);
    
    // number of file columns
    colnum = 2;
    
    timecount_old=0;
	timecount=1;
    
    // read file
    read_format_1(p,pgc);
}

sixdof_motionext_wavemaker::~sixdof_motionext_wavemaker()
{
}

void sixdof_motionext_wavemaker::ini(lexer *p, ghostcell *pgc)
{
    Uext = 0.0;
    Vext = 0.0;
    Wext = 0.0;

    Pext = 0.0;
    Qext = 0.0;
    Rext = 0.0;
}

void sixdof_motionext_wavemaker::motionext_trans(lexer *p, ghostcell *pgc, Eigen::Vector3d& dp_, Eigen::Vector3d& dc_)
{
    // find correct time step
    if((p->simtime>data[timecount][0]))
    timecount_old=timecount;
    
	while(p->simtime>data[timecount][0])
	++timecount;
    
    
        Uext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te && timecount<ptnum-1 && timecount_old<ptnum)
        Uext = (data[timecount][1]-data[timecount_old][1])/(data[timecount][0]-data[timecount_old][0]);
        
        if(p->mpirank==0)
        cout<<"6DOF_motion  Uext "<<Uext<<endl;
        
        dp_(0) = 0.0;
        dc_(0) = Uext*ramp_vel(p);

        dp_(1) = 0.0;
        dc_(1) = 0.0;
        
        dp_(2) = 0.0;
        dc_(2) = 0.0;
}

void sixdof_motionext_wavemaker::motionext_rot(lexer *p, Eigen::Vector3d& dh_, Eigen::Vector3d& h_, Eigen::Vector4d& de_, Eigen::Matrix<double, 3, 4>&G_,  Eigen::Matrix3d&I_)
{
        dh_ << 0.0,0.0,0.0;
        
        omega_ << 0.0, 0.0, 0.0;
        
        h_ = I_*omega_;
        
        de_ = 0.5*G_.transpose()*I_.inverse()*h_;
}

double sixdof_motionext_wavemaker::ramp_vel(lexer *p)
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

double sixdof_motionext_wavemaker::ramp_draft(lexer *p)
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
