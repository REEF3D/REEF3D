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
        
        if(tstage>-1.0e19)
        Uext = stage_vel(1,tstage);
        
        dp_(0) = 0.0;
        dc_(0) = Uext*ramp_vel(p);
    }

    if (dof[1]==2)
    {
        Vext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Vext = (data[timecount][2]-data[timecount_old][2])/(data[timecount][0]-data[timecount_old][0]);
        
        if(tstage>-1.0e19)
        Vext = stage_vel(2,tstage);
        
        dp_(1) = 0.0;
        dc_(1) = Vext*ramp_vel(p);
    }
    
    if (dof[2]==2)
    {
        Wext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Wext = (data[timecount][3]-data[timecount_old][3])/(data[timecount][0]-data[timecount_old][0]);
        
        if(tstage>-1.0e19)
        Wext = stage_vel(3,tstage);
        
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
        
        if(tstage>-1.0e19)
        Pext = stage_vel(1,tstage);
    }

    if(dof[4]==2)
    {
        Qext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Qext = (data[timecount][2]-data[timecount_old][2])/(data[timecount][0]-data[timecount_old][0]);
        
        if(tstage>-1.0e19)
        Qext = stage_vel(2,tstage);
    }

    if(dof[5]==2)
    {
        Rext = 0.0;
        
        if(p->simtime>=ts && p->simtime<=te)
        Rext = (data[timecount][3]-data[timecount_old][3])/(data[timecount][0]-data[timecount_old][0]);
        
        if(tstage>-1.0e19)
        Rext = stage_vel(3,tstage);
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

// The velocity of column col of the file at the time t (the end of a stage of the body integrator):
// the velocity at a sample from the samples before and after it (central difference, one-sided at the
// ends of the file), linear in time between the samples.
//
// The velocity was the difference quotient of the interval of the file at p->simtime, the start of the
// step, for every stage: a step function in time, the same for all stages.  The fluid then took the
// change of the velocity of a whole step in the first stage and none in the others, although the stages
// end at different times (SSP-RK3: t+dt, t+dt/2, t+dt).  With the pressure carried from stage to stage
// the stage pressures alternate (P, 3P, P/2 for a constant acceleration) and the loads of the step see
// about half the pressure of the acceleration (the added mass of a submerged square: 24.6 kg/m with the
// closed projection X 49 against 51-54 kg/m with N 40 5; 66.7 kg/m through the leak of the open one)
double sixdof_motionext_file::stage_vel(int col, double t)
{
    const int n = int(data.size());

    if(n<2 || t<ts || t>te)
    return 0.0;

    kstage = MAX(0,MIN(kstage,n-2));

    while(kstage<n-2 && t>=data[kstage+1][0])
    ++kstage;

    while(kstage>0 && t<data[kstage][0])
    --kstage;

    auto node = [&](int k)
    {
        const int a = MAX(k-1,0), b = MIN(k+1,n-1);
        return (data[b][col]-data[a][col])/(data[b][0]-data[a][0]);
    };

    const double t0 = data[kstage][0], t1 = data[kstage+1][0];
    const double w = (t1>t0) ? MAX(0.0,MIN(1.0,(t-t0)/(t1-t0))) : 0.0;

    return (1.0-w)*node(kstage) + w*node(kstage+1);
}

// The position of column col of the file at the time t, linear between the samples (the interval
// kstage of stage_vel)
double sixdof_motionext_file::stage_pos(int col, double t)
{
    stage_vel(col,t);

    const double t0 = data[kstage][0], t1 = data[kstage+1][0];
    const double w = (t1>t0) ? MAX(0.0,MIN(1.0,(t-t0)/(t1-t0))) : 0.0;

    return (1.0-w)*data[kstage][col] + w*data[kstage+1][col];
}

// The stage of the body integrator begins: the time it ends (tstage), and at the first stage the offset of
// the body to the file, c0 = the position at the start of the stage - the file at that time (before the
// first sample: at the first sample, the body is at rest until then).  The body starts where X 110 / the
// STL put it.  Taken before the stage: after it the body has already moved by dt times the velocity of
// the stage (the offset of the file at the end of the stage was off by that).  The same offset after a
// restart (P 48 / I 45) as long as the body followed the file.
void sixdof_motionext_file::stage_begin(lexer *p, double t, const Eigen::Vector3d &c_)
{
    tstage = t;

    if(c0set || int(data.size())<2 || p->simtime>te)
    return;

    const double ta = MAX(p->simtime,ts);

    for(int q=0; q<3; ++q)
    c0[q] = c_(q) - stage_pos(q+1,ta);

    c0set = true;
}

// The position of the prescribed translations at the time the stage ends, from the file, with the offset
// of stage_begin.  The pose integrated from the velocity at the end of the stage runs ahead of the file by
// up to dt/2 times the change of the velocity (2.5 mm of 30 mm at 80 steps per period); the velocity is the
// target of the forcing, the position is taken from the file.  Not with the ramp of the velocity (X 206)
void sixdof_motionext_file::stage_position(lexer *p, Eigen::Vector3d &c_)
{
    if(!c0set || tstage<-1.0e19 || int(data.size())<2 || tstage<ts || tstage>te || p->X206==1)
    return;

    for(int q=0; q<3; ++q)
    if(dof[q]==2)
    c_(q) = c0[q] + stage_pos(q+1,tstage);
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
