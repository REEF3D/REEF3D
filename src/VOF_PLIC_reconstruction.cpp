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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::reconstructPlane(fdm* a, lexer* p)
{
	//- Calculating interface normal vector n
	
	//calcNormalFO(a, p);
	//calcNormalWENO(a, p);
	//calcNormalPhi(a, p);
	calcNormalLS(a, p);


	//- Scale n_i according to cell size
	/*
	nx(i,j,k) *= p->DXN[IP];
	ny(i,j,k) *= p->DYN[JP];
	nz(i,j,k) *= p->DZN[KP];

	
    //- Ensure that n_i > 0
	
    int invx = -1;
    int invy = -1;
    int invz = -1;
    
    if (nx(i,j,k) < 0.0)
    {
        nx(i,j,k) *= -1.0;
        invx = 1;
    }
    if (ny(i,j,k) < 0.0)
    {
        ny(i,j,k) *= -1.0;
        invy = 1;
    }
    if (nz(i,j,k) < 0.0)
    {
        nz(i,j,k) *= -1.0;
        invz = 1;
    }
    */

	//- Normalise plane

    double sum = sqrt(nx(i,j,k)*nx(i,j,k) + ny(i,j,k)*ny(i,j,k) + nz(i,j,k)*nz(i,j,k));// + pow(p->DXN[IP]*p->DYN[JP]*p->DZN[KP],1.0/3.0)*1e-5;
    
    nx(i,j,k) /= sum;
    ny(i,j,k) /= sum;
    nz(i,j,k) /= sum;		


	//-  Calculating alpha from n and vof according to Scardovelli p.234
		
	alpha(i, j, k) = calcAlpha2(a->vof(i,j,k),p->DXN[IP],p->DYN[JP],p->DZN[KP],nx(i, j, k), ny(i, j, k), nz(i, j, k));
   // alpha(i, j, k) =calcAlpha(a,p,nx(i,j,k),ny(i,j,k),nz(i,j,k));
	

	//- Return to original plane
	/*
	nx(i, j, k) *= -invx;
	ny(i, j, k) *= -invy;
	nz(i, j, k) *= -invz;
     
	
	alpha(i,j,k) += 
		min(0.0, nx(i, j, k)) + min(0.0, ny(i, j, k)) + min(0.0, nz(i, j, k));	
         */
}


double VOF_PLIC::calcAlpha
(
	fdm* a, lexer* p,
	double& nnx,
	double& nny,
	double& nnz
)
{
    double nx,ny,nz;
    if(nnx<0.0)
        nx=-nnx;
    else
        nx=nnx;
        
    if(nny<0.0)
        ny=-nny;
    else
        ny=nny;
        
    if(nnz<0.0)
        nz=-nnz;
    else
        nz=nnz;
    
	//- Assign n_i to m_i such that m1*dx1 < m2*dx2 < m3*dx3
	double m1,m2,m3,dx1,dx2,dx3;
    
    if(nx*p->DXN[IP]<=ny*p->DYN[JP])
    {
        m1=nx;
        dx1=p->DXN[IP];
        m2=ny;
        dx2=p->DYN[JP];
	}
    else
    {
        m1=ny;
        dx1=p->DYN[JP];
        m2=nx;
        dx2=p->DXN[IP];
    }
    
    if(nz*p->DZN[KP]>=m2*dx2)
    {
        m3=nz;
        dx3=p->DZN[KP];
    }
    else if(nz*p->DZN[KP]<m1*dx1)
    {
        m3=m2;
        dx3=dx2;
        m2=m1;
        dx2=dx1;
        m1=nz;
        dx1=p->DZN[KP];
    }
    else
    {
        m3=m2;
        dx3=dx2;
        m2=nz;
        dx2=p->DZN[KP];
    }

	//- Calculate ranges of functions V1, V2, V3
	
	double m12 = m1 + m2;
	double pr = max(6.0*m1*m2*m3, 1.0e-20);
	
	double V1  = pow(m1, 3.0)/pr;
	double V2  = V1 + 0.5*(m2 - m1)/(m3 + 1e-20);
	
	double V3, mm;
	if (m3 < m12)	// if m = m3
	{
		mm = m3;
		V3 = 
			(m3*m3*(3.0*m12 - m3) + m1*m1*(m1 - 3.0*m3) + m2*m2*(m2 - 3.0*m3))
			/pr;
	}
	else	// if m = m12
	{
		mm = m12;
		V3 = 0.5*mm/(m3 + 1e-20);
	}
    
	
	//- Limit vof such that 0 < vof < 0.5
	
	double vofLim = min(a->vof(i, j, k), 1.0 - a->vof(i, j, k));
     	 
	 
	//- Calculate alpha from vofLim and V_i

	double alpha;
	if (vofLim <= V1)	// 0 < V < V1
	{
		alpha = pow(pr*vofLim, 1.0/3.0);
	}
	else if (vofLim < V2)	// V1 < V < V2
	{
		alpha = 0.5*(m1 + sqrt(m1*m1 + 8.0*m2*m3*(vofLim-V1)));
	}
	else if (vofLim < V3)	// V2 < V < V3
	{
		double p = 2.0*m1*m2;
		double cs = cos(acos((1.5*m1*m2*(m12 - 2.0*m3*vofLim))/(p*sqrt(p)))/3.0);
		
		alpha = sqrt(p)*(sqrt(3.0*(1.0 - cs*cs)) - cs) + m12;
	}
	else if (m12 < m3)	// V32 < V < 0.5
	{
         alpha = m3*vofLim + 0.5*mm;
	}
	else	// V31 < V < 0.5
	{
         double p = m1*(m2 + m3) + m2*m3 - 0.25;
         double cs = cos(acos((1.5*m1*m2*m3*(0.5 - vofLim))/(p*sqrt(p)))/3.0);
		 
         alpha = sqrt(p)*(sqrt(3.0*(1.0 - cs*cs)) - cs) + 0.5;
	}


	//- Inverse result if 0.5 < vof < 1.0
	
	if (a->vof(i,j,k) > 0.5)  
	{
		alpha = 1.0 - alpha;
	}

	return alpha;

	
}

double VOF_PLIC::calcAlpha2( double Frac, double dx, double dy, double dz, double& nx, double& ny, double& nz)
{
    //making sure all vector components are positive
    double n_x,n_y,n_z,ret;
    n_x=abs(nx);
    n_y=abs(ny);
    n_z=abs(nz);
    
    //sorting by n_i * d_i
    double n_a,n_b,n_g,d_a,d_b,d_g;
    
    if(n_x*dx<=n_y*dy)
    {
        n_a=n_x;
        d_a=dx;
        n_b=n_y;
        d_b=dy;
	}
    else
    {
        n_a=n_y;
        d_a=dy;
        n_b=n_x;
        d_b=dx;
    }
    
    if(n_z*dz>=n_b*d_b)
    {
        n_g=n_z;
        d_g=dz;
    }
    else if(n_z*dz<n_a*d_a)
    {
        n_g=n_b;
        d_g=d_b;
        n_b=n_a;
        d_b=d_a;
        n_a=n_z;
        d_a=dz;
    }
    else
    {
        n_g=n_b;
        d_g=d_b;
        n_b=n_z;
        d_b=dz;
    }
    
    double m_a,m_b;
    m_a=-n_a/n_g;
    m_b=-n_b/n_g;
    
    //calculate n_crits and Vol_crits
    
    double Vol_,Vol_ab,Vol_bc,Vol_cd,Vol_de,Vol_ef,Vol_fg;
    double n_, n_ab,n_bc,n_cd,n_de,n_ef,n_fg;
    
    n_ab=-m_a*d_a;
    n_bc=-m_b*d_b;
    n_cd=-m_a*d_a-m_b*d_b;
    n_de=d_g;
    n_ef=d_g-m_a*d_a;
    n_fg=d_g-m_b*d_b;
    
    Vol_=Frac*dx*dy*dz;
    
    Vol_ab=(n_ab*n_ab*n_ab)/(6*m_a*m_b);
    
    Vol_bc=-1.0/(2*m_b)*(d_a*n_bc*n_bc+m_a*d_a*d_a/2*n_bc+m_a*m_a*d_a*d_a*d_a/3);
    
    Vol_cd=d_a*d_b*n_cd+0.5*m_a*d_a*d_a*d_b+0.5*m_b*d_a*d_b*d_b;
    
    Vol_de=d_a*d_b*n_de+0.5*m_a*d_a*d_a*d_b+0.5*m_b*d_a*d_b*d_b;
    
    Vol_ef=d_a*d_b*n_ef+0.5*m_a*d_a*d_a*d_b+0.5*m_b*d_a*d_b*d_b
            -((n_ef-d_g)*(n_ef-d_g)*(n_ef-d_g))/(6*m_a*m_b);
            
    Vol_fg=d_a*d_b*n_fg+0.5*m_a*d_a*d_a*d_b+0.5*m_b*d_a*d_b*d_b
            +1.0/(2*m_b)*(d_a*(n_fg-d_g)*(n_fg-d_g)+m_a*d_a*d_a/2*(n_fg-d_g)+m_a*m_a*d_a*d_a*d_a/3);
    
    //check which case the Volume is giving
    int c_ase;
    if(Vol_<=Vol_ab)
        c_ase=0;
    else if(Vol_<=Vol_bc)
        c_ase=1;
    else if(Vol_<=Vol_cd)
        c_ase=2;
    else if(Vol_<=Vol_de)
        c_ase=3;
    else if(Vol_<=Vol_ef)
        c_ase=4;
    else if(Vol_<=Vol_fg)
        c_ase=5;
    else if (Vol_<=d_a*d_b*d_g)
        c_ase=6;
    else
    {
        cout<<"volume too large"<<endl;
        c_ase=6;
    }
    //calculate n_
    double n_tilde,p_cub,q_cub,D_cub,u1_cub,u2_cub;
    double n_tilde_1,n_tilde_2,n_1,n_2,p_quad,q_quad;
    switch(c_ase)
    {
        case 0:
            n_=cbrt(6*m_a*m_b*Vol_);
            break;
            
        case 1:
            n_1=-m_a*d_a/4+sqrt(m_a*m_a*d_a*d_a/16-m_a*m_a*d_a*d_a/3-2*m_b*Vol_/d_a);
            n_2=-m_a*d_a/4-sqrt(m_a*m_a*d_a*d_a/16-m_a*m_a*d_a*d_a/3-2*m_b*Vol_/d_a);
            
            if((n_1>=n_ab && n_1<=n_bc) && (n_2<n_ab || n_2 > n_bc || (n_2<n_1+n_1/1000 && n_2>n_1-n_1/1000)))
                n_=n_1;
            else if((n_2>=n_ab && n_2<=n_bc) && (n_1<n_ab || n_1 > n_bc || (n_1<n_2+n_2/1000 && n_1>n_2-n_2/1000)))
                n_=n_2;
            else
            {
                cout<<"mucho problemo case b: "<<"n_ab="<<n_ab<< "n_bc="<<n_bc<<" n_1="<<n_1<<" n_2="<<n_2<<endl;
                n_=n_ab;
            }
            break;
        
        case 2:
            p_cub=-6*m_a*m_b*d_a*d_b;
            q_cub=3*m_a*m_a*m_b*d_a*d_a*d_b+3*m_a*m_b*m_b*d_a*d_b*d_b+6*m_a*m_b*Vol_;
            D_cub=q_cub*q_cub/4+p_cub*p_cub*p_cub/27;
            if(D_cub>-0.0000001&&D_cub<=0.0)
                D_cub=0.000001;
            u1_cub=-q_cub/2+sqrt(D_cub);
            u2_cub=-q_cub/2-sqrt(D_cub);
            n_tilde=cbrt(u1_cub)+cbrt(u2_cub);
            n_=n_tilde-m_a*d_a-m_b*d_b;
            if(D_cub<=0.0)
            {
                cout<<"mucho problemo case c: "<<"p="<<p_cub<<" q="<<q_cub<<" D="<<D_cub<<endl;
                n_=n_bc;
            }
            break;
            
        case 3:
            n_=Vol_/(d_a*d_b)-0.5*m_a*d_a-0.5*m_b*d_b;
            break;
            
        case 4:
            p_cub=-6*m_a*m_b*d_a*d_b;
            q_cub=-3*m_a*m_a*m_b*d_a*d_a*d_b-3*m_a*m_b*m_b*d_a*d_b*d_b-6*m_a*m_b*d_a*d_b*d_g+6*m_a*m_b*Vol_;
            D_cub=q_cub*q_cub/4+p_cub*p_cub*p_cub/27;
            if(D_cub>-0.000001 && D_cub<=0.0)
                D_cub=0.0000001;
            u1_cub=-q_cub/2+sqrt(D_cub);
            u2_cub=-q_cub/2-sqrt(D_cub);
            n_tilde=cbrt(u1_cub)+cbrt(u2_cub);
            n_=d_g+n_tilde;
            if(D_cub<=0.0)
            {
                cout<<"mucho problemo case e: "<<"p="<<p_cub<<" q="<<q_cub<<" D="<<D_cub<<endl;
                n_=n_de;
            }
            break;
            
        case 5:
            p_quad=0.5*m_a*d_a+2*m_b*d_b;
            q_quad=1/3*d_a*d_a*m_a*m_a+m_b*m_a*d_a+m_b*m_b*d_b*d_b-2*m_b/d_a*Vol_;
            n_tilde_1=-p_quad/2+sqrt(p_quad*p_quad/4-q_quad);
            n_tilde_2=-p_quad/2-sqrt(p_quad*p_quad/4-q_quad);
            n_1=d_g+n_tilde_1;
            n_2=d_g+n_tilde_2;
            
            if((n_1>=n_ef && n_1<=n_fg) && (n_2<n_ef || n_2 > n_fg || (n_2<n_1+n_1/1000 && n_2>n_1-n_1/1000)))
                n_=n_1;
            else if((n_2>=n_ef && n_2<=n_fg) && (n_1<n_ef || n_1 > n_fg || (n_1<n_2+n_2/1000 && n_1>n_2-n_2/1000)))
                n_=n_2;
            else
            {
                cout<<"mucho problemo case f: "<<"n_ef="<<n_ef<< "n_fg="<<n_fg<<" n_1="<<n_1<<" n_2="<<n_2<<endl;
                n_=n_ef;
            }
            break;
            
        case 6:
            n_tilde=cbrt(6*m_a*m_b*(d_a*d_b*d_g-Vol_));
            n_=d_g-m_a*d_a-m_b*d_b-n_tilde;
            break;
    }
    ret=n_*n_g;
    return ret;
}
