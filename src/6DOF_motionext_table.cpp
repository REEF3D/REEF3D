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

#include"6DOF_motionext.h"
#include"lexer.h"

// Prescribed motion from a table (the motion files X 240 1 and X 240 11) at the stages of the body integrator.
//
// sixdof_obj::motion_stage_begin / motion_stage_end call stage_begin before and stage_position after every
// stage (rk2, rk3, rkls3, rk4, rkmid, the one-step forms) with the time the stage ends.  The velocities of the
// table at that time are the derivatives of the stage (motionext_trans, motionext_rot), and the pose at that
// time replaces the integrated one: the body is where the table puts it, with the velocity of the table, at
// the time of the fluid state that sees it next.
//
// Translations: c = c0 + the column of the table.  Rotations, when no rotational DOF is free: the Euler
// angles (Z-Y-X: psi about z, theta about the new y, phi about the new x, as sixdof_rigidbody::euler_angles)
// a0 + the column for a prescribed DOF, a0 for a fixed one, and the angular velocity of the body from the
// rates of these angles.  The offsets c0, a0 are the pose of the body at the start of the first stage minus
// the table at that time: the body starts where X 110 / X 101 / the STL put it.  A body with a free rotation
// and a prescribed one keeps the integrated orientation (the rate of the angle as the component of the
// angular velocity).

// the velocity of column col at the time t: the velocity at a sample from the samples before and after it
// (central difference, one-sided at the ends of the table), linear in time between the samples
double sixdof_motionext::table_vel(int col, double t)
{
    const int n = int(data.size());

    if(n<2 || col<0 || t<ts || t>te)
    return 0.0;

    ktab = MAX(0,MIN(ktab,n-2));

    while(ktab<n-2 && t>=data[ktab+1][0])
    ++ktab;

    while(ktab>0 && t<data[ktab][0])
    --ktab;

    auto node = [&](int k)
    {
        const int a = MAX(k-1,0), b = MIN(k+1,n-1);
        return (data[b][col]-data[a][col])/(data[b][0]-data[a][0]);
    };

    const double t0 = data[ktab][0], t1 = data[ktab+1][0];
    const double w = (t1>t0) ? MAX(0.0,MIN(1.0,(t-t0)/(t1-t0))) : 0.0;

    return (1.0-w)*node(ktab) + w*node(ktab+1);
}

// the value of column col at the time t (clamped to the table), linear between the samples
double sixdof_motionext::table_pos(int col, double t)
{
    const int n = int(data.size());

    if(n<1 || col<0)
    return 0.0;

    if(n<2 || t<=ts)
    return data[0][col];

    if(t>=te)
    return data[n-1][col];

    table_vel(col,t);    // the interval ktab

    const double t0 = data[ktab][0], t1 = data[ktab+1][0];
    const double w = (t1>t0) ? MAX(0.0,MIN(1.0,(t-t0)/(t1-t0))) : 0.0;

    return (1.0-w)*data[ktab][col] + w*data[ktab+1][col];
}

static void quat_to_euler(const Eigen::Vector4d &e, double *ang)
{
    // as sixdof_rigidbody::euler_angles
    ang[2] = atan2(2.0*(e(1)*e(2) + e(3)*e(0)), 1.0 - 2.0*(e(2)*e(2) + e(3)*e(3)));

    const double arg = 2.0*(e(0)*e(2) - e(1)*e(3));
    ang[1] = (fabs(arg)>=1.0) ? (arg>=0.0 ? 1.0 : -1.0)*PI/2.0 : asin(arg);

    ang[0] = atan2(2.0*(e(2)*e(3) + e(1)*e(0)), 1.0 - 2.0*(e(1)*e(1) + e(2)*e(2)));
}

static void euler_to_quat(const double *ang, Eigen::Vector4d &e)
{
    // as sixdof_rigidbody::quaternion_from_euler (Goldstein p. 604)
    const double cf = cos(0.5*ang[0]), sf = sin(0.5*ang[0]);
    const double ct = cos(0.5*ang[1]), st = sin(0.5*ang[1]);
    const double cp = cos(0.5*ang[2]), sp = sin(0.5*ang[2]);

    e(0) = cf*ct*cp + sf*st*sp;
    e(1) = sf*ct*cp - cf*st*sp;
    e(2) = cf*st*cp + sf*ct*sp;
    e(3) = cf*ct*sp - sf*st*cp;
}

void sixdof_motionext::stage_begin(lexer *p, double t, const Eigen::Vector3d &c, const Eigen::Vector4d &e, bool twoD)
{
    tstage = t;
    twoD_ = twoD;
    ramp_ = (p->X206==1);

    // the offsets of the body to the table, once: before the first stage (after it the body has already moved
    // by the stage), at p->simtime or at the first sample when the table starts later (the body rests until
    // then); the same offsets after a restart as long as the body followed the table
    if(offset_set || int(data.size())<2 || p->simtime>te)
    return;

    const double ta = MAX(p->simtime,ts);

    for(int q=0; q<3; ++q)
    c0[q] = c(q) - (tcol[q]>=0 ? table_pos(tcol[q],ta) : 0.0);

    double ang[3];
    quat_to_euler(e,ang);

    for(int q=0; q<3; ++q)
    a0[q] = ang[q] - (table_rot_col(q) ? table_pos(rcol[q],ta) : 0.0);

    offset_set = true;
}

bool sixdof_motionext::table_stage() const
{
    return offset_set && !ramp_ && int(data.size())>=2 && tstage>=ts && tstage<=te;
}

bool sixdof_motionext::table_rotation() const
{
    bool any=false;

    for(int q=0; q<3; ++q)
    {
        if(twoD_ && q!=1)
        continue;

        if(dof[3+q]==1)
        return false;

        if(dof[3+q]==2)
        any = true;
    }

    return any;
}

// rotation q (phi, theta, psi) follows the table: prescribed, with a column (2D: pitch only)
bool sixdof_motionext::table_rot_col(int q) const
{
    return dof[3+q]==2 && rcol[q]>=0 && !(twoD_ && q!=1);
}

void sixdof_motionext::table_angles(double t, double *ang)
{
    for(int q=0; q<3; ++q)
    ang[q] = a0[q] + (table_rot_col(q) ? table_pos(rcol[q],t) : 0.0);
}

// body rates from the rates of the Z-Y-X Euler angles:
//   p = phi' - psi' sin(theta),  q = theta' cos(phi) + psi' sin(phi) cos(theta),
//   r = -theta' sin(phi) + psi' cos(phi) cos(theta)
void sixdof_motionext::table_omega(Eigen::Vector3d &omega)
{
    double ang[3], rate[3];

    table_angles(tstage,ang);

    for(int q=0; q<3; ++q)
    rate[q] = table_rot_col(q) ? table_vel(rcol[q],tstage) : 0.0;

    const double sf = sin(ang[0]), cf = cos(ang[0]);
    const double st = sin(ang[1]), ct = cos(ang[1]);

    omega(0) = rate[0] - rate[2]*st;
    omega(1) = rate[1]*cf + rate[2]*sf*ct;
    omega(2) = -rate[1]*sf + rate[2]*cf*ct;
}

void sixdof_motionext::stage_position(lexer *p, Eigen::Vector3d &c, Eigen::Vector4d &e)
{
    if(!table_stage())
    return;

    for(int q=0; q<3; ++q)
    if(dof[q]==2)
    c(q) = c0[q] + (tcol[q]>=0 ? table_pos(tcol[q],tstage) : 0.0);

    if(table_rotation())
    {
        double ang[3];
        table_angles(tstage,ang);
        euler_to_quat(ang,e);
    }
}
