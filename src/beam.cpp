/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"beam.h"
#include"ckpt_pack.h"

beam::beam(int number):nBeam(number),iout(0),imas(0),ijac(0)
{
    rtoler = nullptr;
    atoler = nullptr;
    
}

beam::~beam()
{
	delete [] z1;
	delete [] z2;
	delete [] z3;
	delete [] y0;
	delete [] scal;
	delete [] f1;
	delete [] f2;
	delete [] f3;
	delete [] cont;
	delete [] ip1;
	delete [] ip2;
	delete rtoler;
	delete atoler;
	delete [] iphes;
	for (int i = 0; i < lde1; i++) {
		delete [] e1[i];
		delete [] e2r[i];
		delete [] e2i[i];
	}
	for (int i = 0; i < ldjac; i++) delete [] fjac[i];
	for (int i = 0; i < ldmas; i++) delete [] fmas[i];
	delete [] e1;
	delete [] e2r;
	delete [] e2i;
	delete [] fjac;
	delete [] fmas;
	delete [] y;
}

void beam::setConstantLoads(Matrix3Xd& Fext_, Matrix4Xd& Mext_, const Matrix3Xd& c_, const Matrix3Xd& cdot_, const Matrix4Xd& q_, const Matrix4Xd& qdot_)
{
    // Set constant external forces
    Fext = Fext_;

    // Set constant external moments
    Mext = Mext_;
}

void beam::setVariableLoads(Matrix3Xd& Fext_, Matrix4Xd& Mext_, const Matrix3Xd& c_, const Matrix3Xd& cdot_, const Matrix4Xd& q_, const Matrix4Xd& qdot_, const double time)
{
    // Add time dependent loads to external forces
    // Add time dependent loads to external moments
}

void beam::setFieldBC(Matrix3Xd& c_, Matrix3Xd& cdot_, Matrix4Xd& q_, Matrix4Xd& q0_, Matrix4Xd& qdot_, Matrix3Xd& f_, Matrix4Xd& m0_, Matrix3Xd& rhs_cdot_, double time , int ind)
{
}

void beam::setStateDot(double *state)
{
    int index = 0;
    for (int i = 0; i < Ne+1; i++) 
    {
        state[index++] = cdot(0,i);
        state[index++] = cdot(1,i);
        state[index++] = cdot(2,i);
    }
    for (int i = 0; i < Ne+2; i++) 
    {
        state[index++] = qdot(0,i);
        state[index++] = qdot(1,i);
        state[index++] = qdot(2,i);
        state[index++] = qdot(3,i);
    }
    for (int i = 0; i < Ne+1; i++) 
    {
        state[index++] = rhs_cdot(0,i);
        state[index++] = rhs_cdot(1,i);
        state[index++] = rhs_cdot(2,i);
    }
    for (int i = 0; i < Ne+2; i++) 
    {
        state[index++] = rhs_qdot(0,i);
        state[index++] = rhs_qdot(1,i);
        state[index++] = rhs_qdot(2,i);
        state[index++] = rhs_qdot(3,i);
    }
}

void beam::setState(double *state)
{
    int index = 0;
    for (int i = 0; i < Ne+1; i++) 
    {
        state[index++] = c(0,i);
        state[index++] = c(1,i);
        state[index++] = c(2,i);
    }
    for (int i = 0; i < Ne+2; i++) 
    {
        state[index++] = q(0,i);
        state[index++] = q(1,i);
        state[index++] = q(2,i);
        state[index++] = q(3,i);
    }
    for (int i = 0; i < Ne+1; i++) 
    {
        state[index++] = cdot(0,i);
        state[index++] = cdot(1,i);
        state[index++] = cdot(2,i);
    }
    for (int i = 0; i < Ne+2; i++) 
    {
        state[index++] = qdot(0,i);
        state[index++] = qdot(1,i);
        state[index++] = qdot(2,i);
        state[index++] = qdot(3,i);
    }
}

void beam::getState(double *state)
{
    int index = 0;
    for (int i = 0; i < Ne+1; i++) 
    {
        c(0,i) = state[index++];
        c(1,i) = state[index++];
        c(2,i) = state[index++];
    }
    for (int i = 0; i < Ne+2; i++) 
    {
        q(0,i) = state[index++];
        q(1,i) = state[index++];
        q(2,i) = state[index++];
        q(3,i) = state[index++];
    }
    for (int i = 0; i < Ne+1; i++) 
    {
        cdot(0,i) = state[index++];
        cdot(1,i) = state[index++];
        cdot(2,i) = state[index++];
    }
    for (int i = 0; i < Ne+2; i++) 
    {
        qdot(0,i) = state[index++];
        qdot(1,i) = state[index++];
        qdot(2,i) = state[index++];
        qdot(3,i) = state[index++];
    }
}

// checkpoints (P 48, I 45)
void beam::ckpt_beam_state(std::vector<double> &v) const
{
    ckpt_put w{v};
    w.i(n_dim);
    w.ad(y,n_dim);
    for(const Matrix3Xd *a : {&Fext,&c,&cdot,&cdotdot,&rhs_cdot,&f})
    w.m(*a);
    for(const Matrix4Xd *a : {&Mext,&q,&qdot,&rhs_qdot,&f0,&m0})
    w.m(*a);
    for(double x : {delta_t,hold,xold,x,x_end,xd,printtime})
    w.d(x);
}

void beam::ckpt_beam_restore(const std::vector<double> &v, size_t &k)
{
    ckpt_get r{v};
    r.k = k;
    int n=0;
    r.i(n);
    if(!r.ok || n!=n_dim)
    return;
    r.ad(y,n_dim);
    for(Matrix3Xd *a : {&Fext,&c,&cdot,&cdotdot,&rhs_cdot,&f})
    r.m(*a);
    for(Matrix4Xd *a : {&Mext,&q,&qdot,&rhs_qdot,&f0,&m0})
    r.m(*a);
    for(double *a : {&delta_t,&hold,&xold,&x,&x_end,&xd,&printtime})
    r.d(*a);
    k = r.k;
}
