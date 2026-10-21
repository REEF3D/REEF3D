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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_strain.h"
#include"nhflow_wall.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_strain::nhflow_strain(lexer *p, fdm_nhf *d)	: nhflow_gradient(p),epsi(p->F45*p->DXM)
{
    p->Darray(PK,p->imax*p->jmax*(p->kmax+2));
    p->Darray(PK0,p->imax*p->jmax*(p->kmax+2));
    p->Darray(PK_b,p->imax*p->jmax*(p->kmax+2));
    
    TQ11=TQ22=TQ33=TQ12=TQ13=TQ23=nullptr;
    
    // extra stresses of the QCR (A 563) or the EARSM (A 559)
    if(p->A563>0.0 || p->A559>0)
    {
    p->Darray(TQ11,p->imax*p->jmax*(p->kmax+2));
    p->Darray(TQ22,p->imax*p->jmax*(p->kmax+2));
    p->Darray(TQ33,p->imax*p->jmax*(p->kmax+2));
    p->Darray(TQ12,p->imax*p->jmax*(p->kmax+2));
    p->Darray(TQ13,p->imax*p->jmax*(p->kmax+2));
    p->Darray(TQ23,p->imax*p->jmax*(p->kmax+2));
    }
}

nhflow_strain::~nhflow_strain()
{
}

void nhflow_strain::wallf_update(lexer *p, fdm_nhf *d, ghostcell *pgc, int *WALLF)
{
    double dist,ks,ut;
    
	LOOP
	WALLF[IJK]=0;
    
    // wall cells: production and dissipation come from the wall function there (nhflow_wall.h:
    // the bed always, side and top walls with A519 2, none with B11 0)
    LOOP
    if(nhflow_turb_wall(p,d,i,j,k,dist,ks,ut)==1)
    WALLF[IJK]=1;
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
        
        // HLES split (A 561 1): the RANS model carries the 3D (bed-generated) turbulence, its production comes from
        // the vertical shear only; the horizontal shear is mixed by the resolved eddies and nu_sgs (hles_update)
        if(p->A561==1)
        {
        s11 = s22 = s33 = s12 = 0.0;
        s13 = dudz(d->U);
        s23 = (p->j_dir==1) ? dvdz(d->V) : 0.0;
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
    val = (1.0/0.85)*(1.0/p->W1)*d->EV0[IJK]*(p->W22*(p->W3 - p->W1)/(p->DZN[KP]*d->WL(i,j)));
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
	s = 2.0*pvdy(p,a);

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

double nhflow_strain::strainterm_v(lexer *p, fdm_nhf *d)
{
    const double uz = dudz(d->U);
    const double vz = (p->j_dir==1) ? dvdz(d->V) : 0.0;
    
	return sqrt(uz*uz + vz*vz);
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

double nhflow_strain::Sij2(lexer *p, fdm_nhf *d)
{
    double s = strainterm(p,d);
    
    s = s*s;
    
    return s;
}

double nhflow_strain::Qij2(lexer *p, fdm_nhf *d)
{
    double q = rotationterm(p,d);
    
    q = q*q;
    
    return q;
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

// ---------------------------------------------------------------------------------------------------------------
// HLES split (A 561 > 0, A 562 C_s nu_b): horizontal and vertical eddy viscosity of the momentum equations
//   vertical   (stresses with a z index: xz, yz, zz and the whole w equation): nu_v = d->EV (k-omega / k-eps)
//   horizontal (stresses xx, yy, xy):                                          nu_h = d->EVH
//     A 561 1: nu_h = nu_t + (C_s Delta_h)^2 |S_h| + nu_b, the RANS production from the vertical shear only (Pk_update)
//     A 561 2: nu_h =        (C_s Delta_h)^2 |S_h| + nu_b, full RANS production
//   |S_h| = sqrt(2 S11^2 + 2 S22^2 + 4 S12^2) (horizontal strain along the z = const plane), Delta_h = sqrt(dx dy)
// With a nu_h well below the isotropic RANS value the horizontal shear layers become unsteady and their large,
// quasi 2D eddies are resolved on the horizontal grid (horizontal large eddy simulation, as Delft3D HLES /
// Hinterberger, Froehlich & Rodi 2007), the vertical resolution is not changed.
// ---------------------------------------------------------------------------------------------------------------
void nhflow_strain::hles_update(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->A561==0)
    return;
    
    LOOP
    {
        double nusgs = 0.0;
        
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
            double sh11 = dudx(d->U);
            double sh22 = 0.0;
            double sh12 = 0.0;
            double dh2  = p->DXN[IP]*p->DXN[IP];
            
            if(p->j_dir==1)
            {
            sh22 = dvdy(d->V);
            sh12 = dudy(d->U) + dvdx(d->V);
            dh2  = (p->geo_curv==1) ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP];
            }
            
            const double Sh = sqrt(2.0*sh11*sh11 + 2.0*sh22*sh22 + sh12*sh12);
            
            nusgs = p->A562_cs*p->A562_cs*dh2*Sh;
        }
        
        if(p->A561==1)
        d->EVH[IJK] = d->EV[IJK] + nusgs + p->A562_nub;
        
        if(p->A561==2)
        d->EVH[IJK] = nusgs + p->A562_nub;
        
        if(p->DF[IJK]<0)
        d->EVH[IJK] = 0.0;
    }
    
    pgc->start24V(p,d->EVH,24);
}

// ---------------------------------------------------------------------------------------------------------------
// QCR (A 563 C_cr1 > 0, quadratic constitutive relation, Spalart 2000): anisotropy of the normal stresses, which
// drives the secondary currents of the second kind (corners, narrow channels); explicit, lagged by one step
// ---------------------------------------------------------------------------------------------------------------
void nhflow_strain::qcr_update(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->A563<=0.0 || TQ11==nullptr)
    return;
    
    double g[3][3],tau[3][3],O[3][3],dt[3][3];
    double *vel[3] = {d->U,d->V,d->W};
    
    LOOP
    {
        TQ11[IJK]=TQ22[IJK]=TQ33[IJK]=TQ12[IJK]=TQ13[IJK]=TQ23[IJK]=0.0;
        
        if(p->wet[IJ]==0 || p->DF[IJK]<0)
        continue;
        
        // g[a][b] = du_a/dx_b (Cartesian, sigma chain rule in nhflow_gradient)
        for(int a=0; a<3; ++a)
        {
            g[a][0] = dudx(vel[a]);
            g[a][1] = (p->j_dir==1) ? dudy(vel[a]) : 0.0;
            g[a][2] = dudz(vel[a]);
        }
        
        if(p->j_dir==0)
        g[1][0] = g[1][2] = 0.0;
        
        double gn = 0.0;
        for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
        gn += g[a][b]*g[a][b];
        
        gn = sqrt(gn);
        
        if(gn<1.0e-10)
        continue;
        
        const double nut = d->EV[IJK];
        
        for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
        {
            tau[a][b] = nut*(g[a][b] + g[b][a]);
            O[a][b]   = (g[a][b] - g[b][a])/gn;
        }
        
        for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
        {
            double s=0.0;
            for(int c=0; c<3; ++c)
            s += O[a][c]*tau[b][c] + O[b][c]*tau[a][c];
            
            dt[a][b] = -p->A563*s;
        }
        
        TQ11[IJK] = dt[0][0];
        TQ22[IJK] = dt[1][1];
        TQ33[IJK] = dt[2][2];
        TQ12[IJK] = dt[0][1];
        TQ13[IJK] = dt[0][2];
        TQ23[IJK] = dt[1][2];
    }
    
    pgc->start24V(p,TQ11,24);
    pgc->start24V(p,TQ22,24);
    pgc->start24V(p,TQ33,24);
    pgc->start24V(p,TQ12,24);
    pgc->start24V(p,TQ13,24);
    pgc->start24V(p,TQ23,24);
}

// momentum is solved for UH: source D (d dtau_xx/dx + d dtau_xy/dy + d dtau_xz/dz) at z = const
double nhflow_strain::qcr_fx(lexer *p, fdm_nhf *d)
{
    if(p->wet[IJ]==0 || p->DF[IJK]<0)
    return 0.0;
    
    return d->WL(i,j)*(dudx(TQ11) + (p->j_dir==1 ? dudy(TQ12) : 0.0) + dudz(TQ13));
}

double nhflow_strain::qcr_fy(lexer *p, fdm_nhf *d)
{
    if(p->wet[IJ]==0 || p->DF[IJK]<0 || p->j_dir==0)
    return 0.0;
    
    return d->WL(i,j)*(dudx(TQ12) + dudy(TQ22) + dudz(TQ23));
}

double nhflow_strain::qcr_fz(lexer *p, fdm_nhf *d)
{
    if(p->wet[IJ]==0 || p->DF[IJK]<0)
    return 0.0;
    
    return d->WL(i,j)*(dudx(TQ13) + (p->j_dir==1 ? dudy(TQ23) : 0.0) + dudz(TQ33));
}

// ---------------------------------------------------------------------------------------------------------------
// EARSM (A 559 1), explicit algebraic Reynolds stress model of Wallin & Johansson (2000), J. Fluid Mech. 403
//   tau = max(1/(Cmu omega), C_tau sqrt(nu/(Cmu k omega))), C_tau = 6
//   S = tau/2 (grad u + grad u^T),  O = tau/2 (grad u - grad u^T)  (O_ij = tau/2 (du_i/dx_j - du_j/dx_i))
//   a = b1 S + b3 (O^2 - IIO/3 I) + b4 (S O - O S) + b6 (S O^2 + O^2 S - IIO S - 2/3 IV I) + b9 (O S O^2 - O^2 S O)
//   N: exact root of the 2D cubic (c1' = 1.8) with the 3D correction, Q = 5/6 (N^2 - 2 IIO)(2 N^2 - IIO)
//   b1 = -N (2 N^2 - IIO)/Q, b3 = -12 IV/(N Q), b4 = -2 (N^2 - 2 IIO)/Q, b6 = -6 N/Q, b9 = 6/Q
//   (b1 with this T6 form is the exact solution of N a = -6/5 S + (a O - O a), checked numerically; it gives
//    Cmu_eff = -b1/2 = 0.089 in equilibrium shear tau*gamma = 3.33. The form b1 = -N(2N^2 - 7 IIO)/Q belongs to
//    T6 without the -IIO S term; mixing the two doubles nu_t in shear)
// Reynolds stress u_i u_j = k (a_ij + 2/3 delta_ij):
//   linear part  -> nu_t = -b1 k tau/2 in the implicit diffusion (d->EV, d->EV0, used for the production too:
//                   the extra anisotropy is production-free in 2D mean flows)
//   extra part   -> TQ_ij = -k a_ex_ij, momentum source D d(TQ_ij)/dx_j (qcr_fx/fy/fz), explicit, lagged
// ---------------------------------------------------------------------------------------------------------------
void nhflow_strain::earsm_update(lexer *p, fdm_nhf *d, ghostcell *pgc, double *KIN, double *EPS)
{
    const double c1 = 1.8;
    const double ctau = 6.0;
    double g[3][3],S[3][3],O[3][3],S2[3][3],O2[3][3],SO[3][3],OS[3][3],SO2[3][3],O2S[3][3],OSO2[3][3],O2SO[3][3],T[3][3];
    double *vel[3] = {d->U,d->V,d->W};
    
    auto mul = [](double A[3][3], double B[3][3], double C[3][3])
    {
        for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
        {
            C[a][b]=0.0;
            for(int c=0; c<3; ++c)
            C[a][b] += A[a][c]*B[c][b];
        }
    };
    
    LOOP
    {
        TQ11[IJK]=TQ22[IJK]=TQ33[IJK]=TQ12[IJK]=TQ13[IJK]=TQ23[IJK]=0.0;
        
        const double kval = MAX(KIN[IJK],0.0);
        const double wval = EPS[IJK];
        
        if(p->wet[IJ]==0 || p->DF[IJK]<0 || wval<1.0e-20 || kval<1.0e-20)
        {
        d->EV0[IJK] = d->EV[IJK] = 0.0001*d->VISC[IJK];
        continue;
        }
        
        // turbulent time scale with the Kolmogorov limit
        const double tau = MAX(1.0/(p->cmu*wval), ctau*sqrt(d->VISC[IJK]/(p->cmu*kval*wval)));
        
        for(int a=0; a<3; ++a)
        {
            g[a][0] = dudx(vel[a]);
            g[a][1] = (p->j_dir==1) ? dudy(vel[a]) : 0.0;
            g[a][2] = dudz(vel[a]);
        }
        
        if(p->j_dir==0)
        g[1][0] = g[1][2] = 0.0;
        
        for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
        {
            S[a][b] = 0.5*tau*(g[a][b] + g[b][a]);
            O[a][b] = 0.5*tau*(g[a][b] - g[b][a]);
        }
        
        // the strain without the trace (incompressible, discrete divergence is not exactly zero)
        const double trS = (S[0][0]+S[1][1]+S[2][2])/3.0;
        S[0][0]-=trS; S[1][1]-=trS; S[2][2]-=trS;
        
        mul(S,S,S2);
        mul(O,O,O2);
        mul(S,O,SO);
        mul(O,S,OS);
        mul(S,O2,SO2);
        mul(O2,S,O2S);
        mul(O,SO2,OSO2);   // O S O^2
        mul(O2S,O,O2SO);   // O^2 S O
        
        double IIS=0.0, IIO=0.0, IV=0.0, V=0.0;
        for(int a=0; a<3; ++a)
        {
            IIS += S2[a][a];
            IIO += O2[a][a];
            IV  += SO2[a][a];
            for(int b=0; b<3; ++b)
            V += S2[a][b]*O2[b][a];
        }
        
        // N: 2D cubic (exact), then the 3D correction
        const double P1 = (c1*c1/27.0 + 0.45*IIS - (2.0/3.0)*IIO)*c1;
        const double P2 = P1*P1 - pow(c1*c1/9.0 + 0.9*IIS + (2.0/3.0)*IIO, 3.0);
        double N;
        
        if(P2>=0.0)
        N = c1/3.0 + cbrt(P1 + sqrt(P2)) + cbrt(P1 - sqrt(P2));
        
        else
        {
        const double arg = MAX(-1.0, MIN(1.0, P1/sqrt(P1*P1 - P2)));
        N = c1/3.0 + 2.0*pow(P1*P1 - P2, 1.0/6.0)*cos(acos(arg)/3.0);
        }
        
        const double D = 20.0*pow(N,4.0)*(N - 0.5*c1) - IIO*(10.0*N*N*N + 15.0*c1*N*N) + 10.0*c1*IIO*IIO;
        
        if(fabs(D)>1.0e-12)
        N += 162.0*(IV*IV + (V - 0.5*IIS*IIO)*N*N)/D;
        
        N = MAX(N, c1/3.0);
        
        const double Q  = (5.0/6.0)*(N*N - 2.0*IIO)*(2.0*N*N - IIO);
        const double b1 = -N*(2.0*N*N - IIO)/Q;
        const double b3 = -12.0*IV/(N*Q);
        const double b4 = -2.0*(N*N - 2.0*IIO)/Q;
        const double b6 = -6.0*N/Q;
        const double b9 = 6.0/Q;
        
        // extra anisotropy
        for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
        {
            const double del = (a==b) ? 1.0 : 0.0;
            
            T[a][b] = b3*(O2[a][b] - IIO*del/3.0)
                    + b4*(SO[a][b] - OS[a][b])
                    + b6*(SO2[a][b] + O2S[a][b] - IIO*S[a][b] - (2.0/3.0)*IV*del)
                    + b9*(OSO2[a][b] - O2SO[a][b]);
        }
        
        double ev = -0.5*b1*kval*tau;
        double kf = kval;
        
        // URANS filter (A 560 22), as the standard k-omega: on nu_t and the extra stress
        if(p->A560==22)
        {
            double dxm;
            
            if(p->j_dir==0)
            dxm = pow(p->DXN[IP]*p->DZN[KP]*d->WL(i,j), (1.0/2.0));

            else
            dxm = pow((p->geo_curv==1 ? p->AREA2D[IJ] : p->DXN[IP]*p->DYN[JP])*p->DZN[KP]*d->WL(i,j), (1.0/3.0));

            const double f = MIN(1.0, dxm*p->A568*p->cmu*wval/sqrt(kval));
            
            ev *= f;
            kf *= f;
        }
        
        d->EV0[IJK] = d->EV[IJK] = MAX(ev, 0.0001*d->VISC[IJK]);
        
        TQ11[IJK] = -kf*T[0][0];
        TQ22[IJK] = -kf*T[1][1];
        TQ33[IJK] = -kf*T[2][2];
        TQ12[IJK] = -kf*T[0][1];
        TQ13[IJK] = -kf*T[0][2];
        TQ23[IJK] = -kf*T[1][2];
    }
    
    pgc->start24V(p,TQ11,24);
    pgc->start24V(p,TQ22,24);
    pgc->start24V(p,TQ33,24);
    pgc->start24V(p,TQ12,24);
    pgc->start24V(p,TQ13,24);
    pgc->start24V(p,TQ23,24);
}
