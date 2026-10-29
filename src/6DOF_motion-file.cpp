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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_motionext_file.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

sixdof_motionext_file::sixdof_motionext_file(lexer *p, ghostcell *pgc, int nb) : id(nb)
{
    ini(p,pgc);
    
    //cout<<"6DOF MOTION START"<<endl;
    
    // number of file columns: X 240 1 is the translation format t, x, y, z (user guide); the
    // columns 1-3 are read by motionext_trans (colnum 3 rejected every 4-column file and read
    // the z column out of range)
    if(p->X240==1)
    colnum = 4;
    
    if(p->X240==2)
    colnum = 4;
    
    if(p->X240==3)
    colnum = 4;
    
    if(p->X240==4)
    colnum = 5;
    
    timecount_old=0;
	timecount=1;
    
    // read file
    //if(p->mpirank==0)
    read_format_1(p,pgc);
}

sixdof_motionext_file::~sixdof_motionext_file()
{
}

void sixdof_motionext_file::ini(lexer *p, ghostcell *pgc)
{
    Uext = 0.0;
    Vext = 0.0;
    Wext = 0.0;

    Pext = 0.0;
    Qext = 0.0;
    Rext = 0.0;
}

void sixdof_motionext_file::motionext_trans(lexer *p, ghostcell *pgc, Eigen::Vector3d& dp_, Eigen::Vector3d& dc_)
{
    // find correct time step
    if (p->simtime <= te)
    {
        if(p->simtime>data[timecount][0])
        timecount_old=timecount;
        
        while(p->simtime>data[timecount][0])
        ++timecount;
    }
    
    
    if(dof[0]==2)
    {
        Uext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Uext = (data[timecount][1]-data[timecount_old][1])/(data[timecount][0]-data[timecount_old][0]);
        
        dp_(0) = 0.0;
        dc_(0) = Uext*ramp_vel(p);
    }

    if (dof[1]==2)
    {
        Vext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Vext = (data[timecount][2]-data[timecount_old][2])/(data[timecount][0]-data[timecount_old][0]);
        
        dp_(1) = 0.0;
        dc_(1) = Vext*ramp_vel(p);
    }
    
    if (dof[2]==2)
    {
        Wext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Wext = (data[timecount][3]-data[timecount_old][3])/(data[timecount][0]-data[timecount_old][0]);
        
        dp_(2) = 0.0;
        dc_(2) = Wext*ramp_vel(p);
    }
}

void sixdof_motionext_file::motionext_rot(lexer *p, Eigen::Vector3d& dh_, Eigen::Vector3d& h_, Eigen::Vector4d& de_, Eigen::Matrix<double, 3, 4>&G_,  Eigen::Matrix3d&I_)
{
    if(dof[3]==2)
    {
        Pext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Pext = (data[timecount][1]-data[timecount_old][1])/(data[timecount][0]-data[timecount_old][0]);
    }

    if(dof[4]==2)
    {
        Qext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Qext = (data[timecount][2]-data[timecount_old][2])/(data[timecount][0]-data[timecount_old][0]);
    }

    if(dof[5]==2)
    {
        Rext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Rext = (data[timecount][3]-data[timecount_old][3])/(data[timecount][0]-data[timecount_old][0]);
    }
    
    if(dof[3] == 2) dh_(0) = 0.0;
    if(dof[4] == 2) dh_(1) = 0.0;
    if(dof[5] == 2) dh_(2) = 0.0;

    omega_ = I_.inverse() * h_;

    if(dof[3] == 2) omega_(0) = Pext*ramp_vel(p);
    if(dof[4] == 2) omega_(1) = Qext*ramp_vel(p);
    if(dof[5] == 2) omega_(2) = Rext*ramp_vel(p);

    h_ = I_*omega_;

    de_ = 0.5*G_.transpose()*omega_;
}

double sixdof_motionext_file::ramp_vel(lexer *p)
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
