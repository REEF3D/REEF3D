/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"fnpf_weno7.h"
#include"lexer.h"
#include"vec.h"
#include"field.h"
#include"slice.h"
#include"fnpf_discrete_weights.h"

fnpf_weno7::fnpf_weno7(lexer* p) :  ddweno_f_nug(p), epsilon(1.0e-10)
{
    p->DYD=p->DXD;

}

fnpf_weno7::~fnpf_weno7()
{
}

double fnpf_weno7::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    grad=0.0;
    
    if(0.5*(ivel1+ivel2)>0.0)
    grad=ddwenox(f,1.0);
    
    if(0.5*(ivel1+ivel2)<0.0)
    grad=ddwenox(f,-1.0);
    
    return grad;
}

double fnpf_weno7::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    grad=0.0;
    
    if(0.5*(jvel1+jvel2)>0.0)
    grad=ddwenoy(f,1.0);
    
    if(0.5*(jvel1+jvel2)<0.0)
    grad=ddwenoy(f,-1.0);
    
    return grad;
}

double fnpf_weno7::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    grad=0.0;
    
/*
    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3] && p->flag4[IJKm4]>0 && p->flag4[IJKm5])
    {
        if(i+p->origin_i>0)
        grad = (-(49.0/20.0)*f(i,j,k+1) + 6.0*f(i,j,k) - 7.5*f(i,j,k-1) + (20.0/3.0)*f(i,j,k-2) - (15.0/4.0)*f(i,j,k-3) + (6.0/5.0)*f(i,j,k-4) - (1.0/6.0)*f(i,j,k-5))
          /(-(49.0/20.0)*p->ZP[KP1] + 6.0*p->ZP[KP] - 7.5*p->ZP[KM1] + (20.0/3.0)*p->ZP[KM2] - (15.0/4.0)*p->ZP[KM3] + (6.0/5.0)*p->ZP[KM4] - (1.0/6.0)*p->ZP[KM5]);
              
        if(i+p->origin_i==0)
        grad = (-(49.0/20.0)*f(i,j,k) + 6.0*f(i,j,k-1) - 7.5*f(i,j,k-2) + (20.0/3.0)*f(i,j,k-3) - (15.0/4.0)*f(i,j,k-4) + (6.0/5.0)*f(i,j,k-5) - (1.0/6.0)*f(i,j,k-6))
          /(-(49.0/20.0)*p->ZP[KP] + 6.0*p->ZP[KM1] - 7.5*p->ZP[KM2] + (20.0/3.0)*p->ZP[KM3] - (15.0/4.0)*p->ZP[KM4] + (6.0/5.0)*p->ZP[KM5] - (1.0/6.0)*p->ZP[KM6]);
              
        //cout<<" return 6"<<endl;
            
        return grad;
    }*/

    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3]>0)
    {
        if(i+p->origin_i>0)
        grad = (-(25.0/12.0)*f(i,j,k+1) + 4.0*f(i,j,k) - 3.0*f(i,j,k-1) + (4.0/3.0)*f(i,j,k-2) - 0.25*f(i,j,k-3))
              /(-(25.0/12.0)*p->ZP[KP1] + 4.0*p->ZP[KP] - 3.0*p->ZP[KM1] + (4.0/3.0)*p->ZP[KM2] - 0.25*p->ZP[KM3]);
              
        if(i+p->origin_i==0)
        grad = (-(25.0/12.0)*f(i,j,k) + 4.0*f(i,j,k-1) - 3.0*f(i,j,k-2) + (4.0/3.0)*f(i,j,k-3) - 0.25*f(i,j,k-4))
              /(-(25.0/12.0)*p->ZP[KP] + 4.0*p->ZP[KM1] - 3.0*p->ZP[KM2] + (4.0/3.0)*p->ZP[KM3] - 0.25*p->ZP[KM4]);
              
        //cout<<" return 4"<<endl;
            
        return grad;
    }
    
    else
    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0)
    {
        if(i+p->origin_i>0)
        grad = (-1.5*f(i,j,k+1) + 2.0*f(i,j,k) - 0.5*f(i,j,k-1))/(-1.5*p->ZP[KP1] + 2.0*p->ZP[KP] - 0.5*p->ZP[KM1]);
              
        if(i+p->origin_i==0)
        grad = (-1.5*f(i,j,k) + 2.0*f(i,j,k-1) - 0.5*f(i,j,k-2))/(-1.5*p->ZP[KP] + 2.0*p->ZP[KM1] - 0.5*p->ZP[KM2]);
            
        //cout<<" return 2"<<endl;    
        
        return grad;
    }
    
    
    else
    {
        if(i+p->origin_i>0)
        grad = (f(i,j,k) - f(i,j,k-1))/(p->ZP[KM1]);
              
        if(i+p->origin_i==0)
        grad = (f(i,j,k) - f(i,j,k-1))/(p->ZP[KM1]);
            
        //cout<<" return 1"<<endl;    
            
        return grad;
    }
}

double fnpf_weno7::sx(lexer *p, slice &f, double ivel)
{
    grad=0.0;
    
    if(ivel>0.0)
    {
	iqmin(p,f);
	is();
	alpha();
	weight();
	
	grad =  w1*(-(1.0/4.0)*q1 + (13.0/12.0)*q2 - (23.0/12.0)*q3 + (25.0/12.0)*q4)
          + w2*((1.0/12.0)*q2 - (5.0/12.0)*q3 + (13.0/12.0)*q4 + (1.0/4.0)*q5)
          + w3*(-(1.0/12.0)*q3 + (7.0/12.0)*q4 + (7.0/12.0)*q5 - (1.0/12.0)*q6)
          + w4*((1.0/4.0)*q4 + (13.0/12.0)*q5 - (5.0/12.0)*q6 + (1.0/12.0)*q7);
	}
    
    if(ivel<0.0)
    {
	iqmax(p,f);
	is();
	alpha();
	weight();
	
	grad =  w1*(-(1.0/4.0)*q1 + (13.0/12.0)*q2 - (23.0/12.0)*q3 + (25.0/12.0)*q4)
          + w2*((1.0/12.0)*q2 - (5.0/12.0)*q3 + (13.0/12.0)*q4 + (1.0/4.0)*q5)
          + w3*(-(1.0/12.0)*q3 + (7.0/12.0)*q4 + (7.0/12.0)*q5 - (1.0/12.0)*q6)
          + w4*((1.0/4.0)*q4 + (13.0/12.0)*q5 - (5.0/12.0)*q6 + (1.0/12.0)*q7);
	}
    
    
    return grad;
}

double fnpf_weno7::sy(lexer *p, slice &f, double jvel)
{
    grad=0.0;
    
    if(jvel>0.0)
    {
	jqmin(p,f);
	is();
	alpha();
	weight();
	
	grad =  w1*(-(1.0/4.0)*q1 + (13.0/12.0)*q2 - (23.0/12.0)*q3 + (25.0/12.0)*q4)
          + w2*((1.0/12.0)*q2 - (5.0/12.0)*q3 + (13.0/12.0)*q4 + (1.0/4.0)*q5)
          + w3*(-(1.0/12.0)*q3 + (7.0/12.0)*q4 + (7.0/12.0)*q5 - (1.0/12.0)*q6)
          + w4*((1.0/4.0)*q4 + (13.0/12.0)*q5 - (5.0/12.0)*q6 + (1.0/12.0)*q7);
	}
    
    if(jvel<0.0)
    {
	jqmax(p,f);
	is();
	alpha();
	weight();
	
	grad =  w1*(-(1.0/4.0)*q1 + (13.0/12.0)*q2 - (23.0/12.0)*q3 + (25.0/12.0)*q4)
          + w2*((1.0/12.0)*q2 - (5.0/12.0)*q3 + (13.0/12.0)*q4 + (1.0/4.0)*q5)
          + w3*(-(1.0/12.0)*q3 + (7.0/12.0)*q4 + (7.0/12.0)*q5 - (1.0/12.0)*q6)
          + w4*((1.0/4.0)*q4 + (13.0/12.0)*q5 - (5.0/12.0)*q6 + (1.0/12.0)*q7);
	}
    
    return grad;   
}

double fnpf_weno7::sz(lexer *p, double *f)
{
    grad = (-(49.0/20.0)*f[FIJK] + 6.0*f[FIJKm1] - 7.5*f[FIJKm2] + (20.0/3.0)*f[FIJKm3] - 3.75*f[FIJKm4] + (6.0/5.0)*f[FIJKm5] - (1.0/6.0)*f[FIJKm6])
          /(-(49.0/20.0)*p->ZN[KP] + 6.0*p->ZN[KM1] - 7.5*p->ZN[KM2] + (20.0/3.0)*p->ZN[KM3] - 3.75*p->ZN[KM4] + (6.0/5.0)*p->ZN[KM5] - (1.0/6.0)*p->ZN[KM6]);
    
    return grad;   
}

// --------------

void fnpf_weno7::iqmin(lexer *p, slice& f)
{	
    q1 = (f(i-3,j) - f(i-4,j))/p->DXD;
	q2 = (f(i-2,j) - f(i-3,j))/p->DXD;
	q3 = (f(i-1,j) - f(i-2,j))/p->DXD;
	q4 = (f(i,j)   - f(i-1,j))/p->DXD;
	q5 = (f(i+1,j) - f(i,j  ))/p->DXD;
    q6 = (f(i+2,j) - f(i+1,j))/p->DXD;
    q7 = (f(i+3,j) - f(i+2,j))/p->DXD;
}

void fnpf_weno7::jqmin(lexer *p, slice& f)
{	
    q1 = (f(i,j-3) - f(i,j-4))/p->DYD;
	q2 = (f(i,j-2) - f(i,j-3))/p->DYD;
	q3 = (f(i,j-1) - f(i,j-2))/p->DYD;
	q4 = (f(i,j)   - f(i,j-1))/p->DYD;
	q5 = (f(i,j+1) - f(i,j  ))/p->DYD;
    q6 = (f(i,j+2) - f(i,j+1))/p->DYD;
    q7 = (f(i,j+3) - f(i,j+2))/p->DYD;
}

void fnpf_weno7::iqmax(lexer *p, slice& f)
{	
    q1 = (f(i+4,j) - f(i+3,j))/p->DXD;
	q2 = (f(i+3,j) - f(i+2,j))/p->DXD;
	q3 = (f(i+2,j) - f(i+1,j))/p->DXD;
	q4 = (f(i+1,j) - f(i,j  ))/p->DXD;
	q5 = (f(i,j)   - f(i-1,j))/p->DXD;
    q6 = (f(i-1,j) - f(i-2,j))/p->DXD;
    q7 = (f(i-2,j) - f(i-3,j))/p->DXD;
}

void fnpf_weno7::jqmax(lexer *p, slice& f)
{	
	q1 = (f(i,j+4) - f(i,j+3))/p->DYD;
    q2 = (f(i,j+3) - f(i,j+2))/p->DYD;
	q3 = (f(i,j+2) - f(i,j+1))/p->DYD;
	q4 = (f(i,j+1) - f(i,j  ))/p->DYD;
	q5 = (f(i,j)   - f(i,j-1))/p->DYD;
    q6 = (f(i,j-1) - f(i,j-2))/p->DYD;
    q7 = (f(i,j-2) - f(i,j-3))/p->DYD;
}

void fnpf_weno7::is()
{
	is1 = q1*(547.0*q1 - 3882.0*q2 + 4642.0*q3 - 1854.0*q4) + q2*(7043.0*q2 - 17246.0*q3 + 7042.0*q4) 
        + q3*(11003.0*q3 - 9402.0*q4) + 2107.0*q4*q4;
        
    is2 = q2*(267.0*q2 - 1642.0*q3 + 1602.0*q4 - 494.0*q5) + q3*(2843.0*q3 - 5966.0*q4 + 1922.0*q5) 
        + q4*(3443.0*q4 - 2522.0*q5) + 547.0*q5*q5;
        
    is3 = q3*(547.0*q3 - 2522.0*q4 + 1922.0*q5 - 494.0*q6) + q4*(3443.0*q4 - 5966.0*q5 + 1602.0*q6) 
        + q5*(2843.0*q2 - 1642.0*q3) + 267.0*q6*q6;
        
    is4 = q4*(2107.0*q4 - 9402.0*q5 + 7042.0*q6 - 1854.0*q7) + q5*(11003.0*q5 - 17246.0*q6 + 4642.0*q7) 
        + q6*(7043.0*q6 - 3882.0*q7) + 547.0*q7*q7;
}

void fnpf_weno7::alpha()
{
	alpha1=(1.0/35.0)/pow(epsilon+is1,2.0);
	alpha2=(12.0/35.0)/pow(epsilon+is2,2.0);
	alpha3=(18.0/35.0)/pow(epsilon+is3,2.0);
    alpha4=(4.0/35.0)/pow(epsilon+is4,2.0);
}

void fnpf_weno7::weight()
{
	w1=alpha1/(alpha1+alpha2+alpha3+alpha4);
	w2=alpha2/(alpha1+alpha2+alpha3+alpha4);
	w3=alpha3/(alpha1+alpha2+alpha3+alpha4);
    w4=alpha4/(alpha1+alpha2+alpha3+alpha4);
}

