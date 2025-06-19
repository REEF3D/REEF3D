/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_strain.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_strain::nhflow_strain(lexer *p, fdm_nhf *d)	: nhflow_gradient(p),epsi(p->F45*p->DXM)
{
    p->Darray(PK,p->imax*p->jmax*(p->kmax+2));
    p->Darray(PK0,p->imax*p->jmax*(p->kmax+2));
    p->Darray(PK_b,p->imax*p->jmax*(p->kmax+2));
    
}

nhflow_strain::~nhflow_strain()
{
}

void nhflow_strain::wallf_update(lexer *p, fdm_nhf *d, ghostcell *pgc, int *WALLF)
{

	LOOP
	WALLF[IJK]=0;
    
    LOOP
    if(p->DF[IJK]>0)
    {  
        if(p->flag4[Im1JK]<0 && p->IO[Im1JK]!=1)
        WALLF[IJK]=1;

        if(p->flag4[Ip1JK]<0  && p->IO[Ip1JK]!=2)
        WALLF[IJK]=1;
        
        if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0)
        WALLF[IJK]=1;
        
        if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0)
        WALLF[IJK]=1;
        
        if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
        WALLF[IJK]=1;

        if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
        WALLF[IJK]=1;
    }
}

void nhflow_strain::Pk_update(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    LOOP
    {
        if(p->j_dir==1)
        {
        s11 = dudx(d->U);
        s22 = dvdy(d->V);
        s33 = dwdz(d->W);
        s12 = (dudy(d->U) + dvdx(d->V));
        s13 = (dudz(d->U) + dwdx(d->W));
        s23 = (dvdz(d->V) + dwdy(d->W));
        }
        
        if(p->j_dir==0)
        {
        s11 = dudx(d->U);
        s22 = 0.0;
        s33 = dwdz(d->W);
        s12 = 0.0;
        s13 = (dudz(d->U) + dwdx(d->W));
        s23 = 0.0;
        }

        PK0[IJK] = d->EV0[IJK]*(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);
        PK[IJK]  =  d->EV[IJK]*(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);
    }
}

void nhflow_strain::Pk_b_update(lexer *p, fdm_nhf *d, ghostcell *pgc)
{ 
    if(p->A566==1)
    LOOP
    {
        
    val = 0.0;
    
    if(k==p->knoz-1)
    val = (1.0/0.85)*(1.0/p->W1)*d->EV0[IJK]*(p->W22*(p->W3 - p->W1)/(p->DZP[KP1]*d->WL(i,j)));
    PK_b[IJK] = val;
    }
}

double nhflow_strain::sij(lexer *p, fdm_nhf *d, int ii, int jj)
{
	double s=0.0;
/*
	if(ii==1 && jj==1)
	s = 2.0*pudx(p,a);

	if((ii==1 && jj==2) || (ii==2 && jj==1))
	s = pudy(p,a) + pvdx(p,a);

	if((ii==1 && jj==3) ||( ii==3 && jj==1))
	s = pudz(p,a) + pwdx(p,a);

	if(ii==2 && jj==2)
	s = pvdy(p,a);

	if((ii==2 && jj==3) || (ii==3 && jj==2))
	s = pvdz(p,a) + pwdy(p,a);

	if(ii==3 && jj==3)
	s = 2.0*pwdz(p,a);
*/
	return 0.5*s;
}


double nhflow_strain::qij(lexer *p, fdm_nhf *d, int ii, int jj)
{
	double q=0.0;
/*
	if((ii==1 && jj==1) || (ii==2 && jj==2) || (ii==3 && jj==3))
	q = 0.0;

	if(ii==1 && jj==2)
	q = pudy(p,a) - pvdx(p,a);

	if(ii==2 && jj==1)
	q = -pudy(p,a) + pvdx(p,a);

	if(ii==1 && jj==3)
	q = pudz(p,a) - pwdx(p,a);

	if(ii==3 && jj==1)
	q = -pudz(p,a) + pwdx(p,a);

	if(ii==2 && jj==3)
	q = pvdz(p,a) - pwdy(p,a);

	if(ii==3 && jj==2)
	q = -pvdz(p,a) + pwdy(p,a);
*/
	return q;
}

double nhflow_strain::pk(lexer *p, fdm_nhf *d)
{ 
    return PK[IJK];
}
double nhflow_strain::strainterm(lexer *p, fdm_nhf *d)
{
	double s=0.0;
    
    if(p->j_dir==1)
    {
	s11 = dudx(d->U);
	s22 = dvdy(d->V);
	s33 = dwdz(d->W);
	s12 = (dudy(d->U) + dvdx(d->V));
	s13 = (dudz(d->U) + dwdx(d->W));
	s23 = (dvdz(d->V) + dwdy(d->W));
    }
    
    if(p->j_dir==0)
    {
	s11 = dudx(d->U);
	s22 = 0.0;
	s33 = dwdz(d->W);
	s12 = 0.0;
	s13 = (dudz(d->U) + dwdx(d->W));
	s23 = 0.0;
    }

    s = sqrt(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);
    
	return s;
}

double nhflow_strain::strainterm(lexer *p, double *U, double *V, double *W)
{
	double s=0.0;
    
    if(p->j_dir==1)
    {
	s11 = dudx(U);
	s22 = dvdy(V);
	s33 = dwdz(W);
	s12 = (dudy(U) + dvdx(V));
	s13 = (dudz(U) + dwdx(W));
	s23 = (dvdz(V) + dwdy(W));
    }
    
    if(p->j_dir==0)
    {
	s11 = dudx(U);
	s22 = 0.0;
	s33 = dwdz(W);
	s12 = 0.0;
	s13 = (dudz(U) + dwdx(W));
	s23 = 0.0;
    }

    s = sqrt(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);
    
	return s;
}

double nhflow_strain::rotationterm(lexer *p, fdm_nhf *d)
{
	double r=0.0;
    
    if(p->j_dir==1)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = (dudy(d->U) - dvdx(d->V));
	r13 = (dudz(d->U) - dwdx(d->W));
	r23 = (dvdz(d->V) - dwdy(d->W));
    }
    
    if(p->j_dir==0)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = 0.0;
	r13 = (dudz(d->U) - dwdx(d->W));
	r23 = 0.0;
    }

    r = sqrt(r12*r12 + r13*r13 + r23*r23);

	return r;
}

double nhflow_strain::rotationterm(lexer *p, double *U, double *V, double *W)
{
	double r=0.0;
    
    if(p->j_dir==1)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = (dudy(U) - dvdx(V));
	r13 = (dudz(U) - dwdx(W));
	r23 = (dvdz(V) - dwdy(W));
    }
    
    if(p->j_dir==0)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = 0.0;
	r13 = (dudz(U) - dwdx(W));
	r23 = 0.0;
    }

    r = sqrt(2.0*r11*r11 + 2.0*r22*r22 + 2.0*r33*r33 + r12*r12 + r13*r13 + r23*r23);

	return r;
}

double nhflow_strain::magSqrSd(lexer *p, fdm_nhf *d)
{
	double Sd=0.0;
	/*double IV_SR=0.0;
	double Strain=0.0;	
	double Omega=0.0;
	

	if(p->j_dir==1)
    {
	s11 = pudx(p,a);
	s22 = pvdy(p,a);
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a));
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a));
    }
    
    if(p->j_dir==0)
    {
	s11 = pudx(p,a);
	s22 = 0.0;
	s33 = pwdz(p,a);
	s12 = 0.0;
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = 0.0;
    }
	
	if(p->j_dir==1)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = (pudy(p,a) - pvdx(p,a));
	r13 = (pudz(p,a) - pwdx(p,a));
	r23 = (pvdz(p,a) - pwdy(p,a));
    }
    
    if(p->j_dir==0)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = 0.0;
	r13 = (pudz(p,a) - pwdx(p,a));
	r23 = 0.0;
    }
	
	ss11 = (s11*s11 + 0.25*s12*s12 + 0.25*s13*s13);
	ss22 = (0.25*s12*s12 + s22*s22 + 0.25*s23*s23);
	ss33 = (0.25*s13*s13 + 0.25*s23*s23 + s33*s33);
	ss12 = (0.5*s11*s12 + 0.5*s12*s22 + 0.25*s13*s23);
	ss13 = (0.5*s11*s13 + 0.25*s12*s23 + 0.5*s13*s33);
	ss23 = (0.25*s12*s13 + 0.5*s22*s23 + 0.5*s23*s33);
	
	rr11 = -0.25*(r12*r12 + r13*r13);
	rr22 = -0.25*(r12*r12 + r23*r23);
	rr33 = -0.25*(r13*r13 + r23*r23);
	rr12 = -0.25*r13*r23;
	rr13 = 0.25*r12*r23;
	rr23 = -0.25*r12*r13;
	
	IV_SR = ss11*rr11 + 2.0*ss12*rr12 + 2.0*ss13*rr13 + ss22*rr22 + 2.0*ss23*rr23 + ss33*rr33;	

	Strain = nhflow_strainterm(p,a);	
	Omega = rotationterm(p,a);
	
    Sd = ((1.0/6.0)*((pow(Strain, 2.0)*pow(Strain, 2.0)) + (pow(Omega, 2.0)*pow(Omega, 2.0)))) + ((2.0/3.0)*(pow(Strain, 2.0)*pow(Omega, 2.0))) + (2.0*IV_SR);				
//    Sd = ((1.0/6.0)*((pow(nhflow_strainterm(p,a), 2.0)*pow(nhflow_strainterm(p,a), 2.0)) + (pow(rotationterm(p,a), 2.0)*pow(rotationterm(p,a), 2.0)))) + ((2.0/3.0)*(pow(nhflow_strainterm(p,a), 2.0)*pow(rotationterm(p,a), 2.0))) + (2.0*IV_SR);
*/
	return Sd;
}

double nhflow_strain::magSqrSd(lexer *p, double *U, double *V, double *W)
{
    double Sd=0.0;
	/*
	double IV_SR=0.0;
	double Strain=0.0;	
	double Omega=0.0;
	

	if(p->j_dir==1)
    {
	s11 = pudx(p,u);
	s22 = pvdy(p,v);
	s33 = pwdz(p,w);
	s12 = (pudy(p,u) + pvdx(p,v));
	s13 = (pudz(p,u) + pwdx(p,w));
	s23 = (pvdz(p,v) + pwdy(p,w));
    }
    
    if(p->j_dir==0)
    {
	s11 = pudx(p,u);
	s22 = 0.0;
	s33 = pwdz(p,w);
	s12 = 0.0;
	s13 = (pudz(p,u) + pwdx(p,w));
	s23 = 0.0;
    }
	
	if(p->j_dir==1)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = (pudy(p,u) - pvdx(p,v));
	r13 = (pudz(p,u) - pwdx(p,w));
	r23 = (pvdz(p,v) - pwdy(p,w));
    }
    
    if(p->j_dir==0)
    {
	r11 = 0.0;
	r22 = 0.0;
	r33 = 0.0;
	r12 = 0.0;
	r13 = (pudz(p,u) - pwdx(p,w));
	r23 = 0.0;
    }
	
	ss11 = (s11*s11 + 0.25*s12*s12 + 0.25*s13*s13);
	ss22 = (0.25*s12*s12 + s22*s22 + 0.25*s23*s23);
	ss33 = (0.25*s13*s13 + 0.25*s23*s23 + s33*s33);
	ss12 = (0.5*s11*s12 + 0.5*s12*s22 + 0.25*s13*s23);
	ss13 = (0.5*s11*s13 + 0.25*s12*s23 + 0.5*s13*s33);
	ss23 = (0.25*s12*s13 + 0.5*s22*s23 + 0.5*s23*s33);
	
	rr11 = -0.25*(r12*r12 + r13*r13);
	rr22 = -0.25*(r12*r12 + r23*r23);
	rr33 = -0.25*(r13*r13 + r23*r23);
	rr12 = -0.25*r13*r23;
	rr13 = 0.25*r12*r23;
	rr23 = -0.25*r12*r13;
	
	IV_SR = ss11*rr11 + 2.0*ss12*rr12 + 2.0*ss13*rr13 + ss22*rr22 + 2.0*ss23*rr23 + ss33*rr33;

	Strain = nhflow_strainterm(p,u,v,w);	
	Omega = rotationterm(p,u,v,w);
	
	Sd = ((1.0/6.0)*((pow(Strain, 2.0)*pow(Strain, 2.0)) + (pow(Omega, 2.0)*pow(Omega, 2.0)))) + ((2.0/3.0)*(pow(Strain, 2.0)*pow(Omega, 2.0))) + (2.0*IV_SR);		
//    Sd = ((1.0/6.0)*((pow(nhflow_strainterm(p,u,v,w), 2.0)*pow(nhflow_strainterm(p,u,v,w), 2.0)) + (pow(rotationterm(p,u,v,w), 2.0)*pow(rotationterm(p,u,v,w), 2.0)))) + ((2.0/3.0)*(pow(nhflow_strainterm(p,u,v,w), 2.0)*pow(rotationterm(p,u,v,w), 2.0))) + (2.0*IV_SR);
*/
	return Sd;
}


double nhflow_strain::strainplain(lexer *p, fdm_nhf *d)
{
	double s=0.0;
/*
	s11 = pudx(p,a);
	s22 = pvdy(p,a)*p->y_dir;
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a))*p->y_dir;
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a))*p->y_dir;

    s=fabs(s11)+fabs(s22)+fabs(s33)+0.5*fabs(s12)+0.5*fabs(s13)+0.5*fabs(s13);
*/
	return s;
}
