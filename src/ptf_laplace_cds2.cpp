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

#include"ptf_laplace_cds2.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"solver.h"
#include"convection.h"
#include"fnpf_weno5.h"

ptf_laplace_cds2::ptf_laplace_cds2(lexer *p, fdm_ptf* a, ghostcell *pgc) : bc(p)
{
    // bc ini
    SLICELOOP4
	bc(i,j) = 0;

    pgc->gcsl_start4int(p,bc,50);

    if(p->B98>=3)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];

    bc(i-1,j) = 1;
    }

    if(p->B99>=3)
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];

    bc(i+1,j) = 2;
    }
}

ptf_laplace_cds2::~ptf_laplace_cds2()
{
}

void ptf_laplace_cds2::start(lexer* p, fdm_ptf* a, ghostcell *pgc, solver *psolv, field &f, slice &Fifsf, slice &eta)
{
    


	n=0;
    FLUIDLOOP
	{

    if(p->flag4[IJK]>0)
    {
	a->M.p[n]  =  1.0/(p->DXP[IP]*p->DXN[IP])
                + 1.0/(p->DXP[IM1]*p->DXN[IP])

                + 1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                + 1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir

                + 1.0/(p->DZP[KP]*p->DZN[KP])
                + 1.0/(p->DZP[KM1]*p->DZN[KP]);

   	a->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP]);
	a->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP]);

	a->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
	a->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

	a->M.t[n] = -1.0/(p->DZP[KP]*p->DZN[KP]);
	a->M.b[n] = -1.0/(p->DZP[KM1]*p->DZN[KP]);

	a->rhsvec.V[n] = 0.0;
    }
    
    

    if(p->flag4[IJK]<0)
    {
	a->M.p[n] =  1.0;

   	a->M.n[n] = 0.0;
	a->M.s[n] = 0.0;

	a->M.w[n] = 0.0;
	a->M.e[n] = 0.0;

	a->M.t[n] = 0.0;
	a->M.b[n] = 0.0;

	a->rhsvec.V[n] = 0.0;
    }
    
	++n;
	}

    // Boundary Conditions
    n=0;
	FLUIDLOOP
	{
        if(p->flag4[IJK]>0)
        {

            // south
            if(p->flag4[Im1JK]<AIR && bc(i-1,j)==0)
            {        
            a->M.p[n] += a->M.s[n];
            a->M.s[n] = 0.0;
            }
            

            if(p->flag4[Im1JK]<AIR && bc(i-1,j)==1)
            {
            a->rhsvec.V[n] += a->M.s[n]*a->u(i-1,j,k)*p->DXP[IP];
            a->M.p[n] += a->M.s[n];
            a->M.s[n] = 0.0;
            }

            // north
            if(p->flag4[Ip1JK]<AIR && bc(i+1,j)==0)
            {
           // if(p->A323==78 && p->flag4[Ip1Jp1K]>AIR)
           //     cout<<"SW corner going north, eta S: "<<eta(i,j)<<" eta P: "<<eta(i+1,j)<<"eta W: "<<eta(i+1,j+1)<<" || Fifsf S: "<<Fifsf(i,j)<<" Fifsf P: "<<Fifsf(i+1,j)<<" Fifsf W: "<<Fifsf(i+1,j+1)<<endl;
              //  cout<<"eta S: "<<eta(i,j)<<" eta P: "<<eta(i+1,j)<<" || Fifsf S: "<<Fifsf(i,j)<<" Fifsf P: "<<Fifsf(i+1,j)<<endl;
                
            a->M.p[n] += a->M.n[n];
            a->M.n[n] = 0.0;
            }

            if(p->flag4[Ip1JK]<AIR && bc(i+1,j)==2)
            {
            a->rhsvec.V[n] -= a->M.n[n]*a->u(i+1,j,k)*p->DXP[IP1];
            a->M.p[n] += a->M.n[n];
            a->M.n[n] = 0.0;
            }

            // east
            if(p->flag4[IJm1K]<AIR)
            {
            //if(p->A323==78 && p->flag4[Im1Jm1K]>AIR)
             //   cout<<"SW corner going east, eta W: "<<eta(i,j)<<" eta P: "<<eta(i,j-1)<<"eta S: "<<eta(i-1,j-1)<<" || Fifsf W: "<<Fifsf(i,j)<<" Fifsf P: "<<Fifsf(i,j-1)<<" Fifsf S: "<<Fifsf(i-1,j-1)<<endl;
            a->M.p[n] += a->M.e[n];
            a->M.e[n] = 0.0;
            }

            // west
            if(p->flag4[IJp1K]<AIR)
            {
            a->M.p[n] += a->M.w[n];
            a->M.w[n] = 0.0;
            }
            
            // top
            if(p->flag4[IJKp1]<AIR)
            {
            a->M.p[n] += a->M.t[n];
            a->M.t[n] = 0.0;
            }

        // FSFBC
            // south
            if(p->flag4[Im1JK]==AIR)
            {
                // -----------
                if(p->A323==1)
                {
                a->rhsvec.V[n] -= a->M.s[n]*f(i-1,j,k);
                a->M.s[n] = 0.0;
                }
                
                // -----------
                if(p->A323==2)
                {
                double lsv0,lsv1,Fival;

                lsv0 = fabs(a->phi(i,j,k));
                lsv1 = fabs(a->phi(i-1,j,k));

                lsv0 = fabs(lsv0)>1.0e-6?lsv0:1.0e20;
                
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k)));
                Fival = teta*Fifsf(i-1,j) + (1.0-teta)*Fifsf(i,j);

                a->rhsvec.V[n] -= a->M.s[n]*Fival*(1.0 + lsv1/lsv0);
                a->M.p[n] -= a->M.s[n]*lsv1/lsv0;
                a->M.s[n] = 0.0;
                }
                
                // -----------
                if(p->A323==3)
                {
                double x0,x1,x2,y2;
                double x,y,Fival;
                double Lx0,Lx1,Lx2;
                
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k)));
                Fival = teta*Fifsf(i-1,j) + (1.0-teta)*Fifsf(i,j);

                x0 = -fabs(a->phi(i+1,j,k));
                x1 = -fabs(a->phi(i,j,k));
                x2 = 0.0;
                
                y2 = Fival;

                x = fabs(a->phi(i-1,j,k));

                Lx0 = ((x-x1)/(x0-x1)) * ((x-x2)/(x0-x2));
                Lx1 = ((x-x0)/(x1-x0)) * ((x-x2)/(x1-x2));
                Lx2 = ((x-x0)/(x2-x0)) * ((x-x1)/(x2-x1));

                a->rhsvec.V[n]  -= a->M.s[n]*Lx2*y2;
                a->M.p[n]       += a->M.s[n]*Lx1;
                a->M.n[n]       += a->M.s[n]*Lx0;
                a->M.s[n]       = 0.0;
                }
                
                // -----------
                if(p->A323==4)
                {
                double Fival,teta;

                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k)));
                Fival = teta*Fifsf(i-1,j) + (1.0-teta)*Fifsf(i,j);
                
                teta = teta>1.0e-6?teta:1.0e20;
                
                a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                a->rhsvec.V[n] -= a->M.s[n]*Fival;
                a->M.s[n] = 0.0;
                }
                
              
                if(p->A323==5)
                {
                    double Fival,teta,zpos_p,zpos_s,x_pos;
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    x_pos=((0.0-zpos_p)*p->DXP[IM1])/(zpos_s-zpos_p);
                    
                    teta = x_pos/p->DXP[IM1];
                    Fival = teta*Fifsf(i-1,j) + (1.0-teta)*Fifsf(i,j);
                
                    teta = teta>1.0e-6?teta:1.0e20;
                
                    a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                    a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.s[n]*Fival;
                    a->M.s[n] = 0.0;
                  
                }
                
                if(p->A323==6 || p->A323==7 || p->A323==77)
                {
                    double zpos_p,zpos_s,zpos_n;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    if(x1_quad<=0.0 && x1_quad>=0.0-(p->DXP[IM1]) && !(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1]))
                        x_quad=x1_quad;
                    else if(x2_quad<=0.0 && x2_quad>=0.0-(p->DXP[IM1]) && !(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=-p->DXP[IM1];
                    }
                    teta=(0.0-x_quad)/p->DXP[IM1];
                    teta = teta>1.0e-6?teta:1.0e20;
                    Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                
                    a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                    a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.s[n]*Fival;
                    a->M.s[n] = 0.0;
                }
                
                if(p->A323==78)
                {
                    double zpos_p,zpos_s,zpos_ss,zpos_n;
                    double s_Fival,p_Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double s_a_quad,s_b_quad,s_b_quad_num,s_b_quad_denom,s_c_quad,s_x1_quad,s_x2_quad,s_x_quad,s_p_quad,s_q_quad;
                    double s_a_fi,s_b_fi_denom,s_b_fi_num,s_c_fi,s_b_fi;
                    double p_a_quad,p_b_quad,p_b_quad_num,p_b_quad_denom,p_c_quad,p_x1_quad,p_x2_quad,p_x_quad,p_p_quad,p_q_quad;
                    double p_a_fi,p_b_fi_denom,p_b_fi_num,p_c_fi,p_b_fi;
                    double x_quad, Fival;
                    
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_ss=eta(i-2,j)-(p->ZP[KP]-p->F60);
                        
                    s_b_quad_denom=(-1.0*p->DXP[IM1])+(-1.0*p->DXP[IM2])-((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(-1.0*p->DXP[IM1]);
                    s_b_quad_num=zpos_ss-zpos_p+zpos_p*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1])-zpos_s*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1]);
                    s_b_quad=s_b_quad_num/s_b_quad_denom;
                    s_c_quad=zpos_p;
                    s_a_quad=zpos_s/(p->DXP[IM1]*p->DXP[IM1])-zpos_p/(p->DXP[IM1]*p->DXP[IM1])-s_b_quad/(-1.0*p->DXP[IM1]);
                    s_p_quad=s_b_quad/s_a_quad;
                    s_q_quad=s_c_quad/s_a_quad;
                    s_x1_quad=(0.0-s_p_quad)/2.0+sqrt(((s_p_quad/2.0)*(s_p_quad/2.0)-s_q_quad));
                    s_x2_quad=(0.0-s_p_quad)/2.0-sqrt(((s_p_quad/2.0)*(s_p_quad/2.0)-s_q_quad));
                        
                    s_b_fi_denom=(-1.0*p->DXP[IM1])+(-1.0*p->DXP[IM2])-((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(-1.0*p->DXP[IM1]);
                    s_b_fi_num=Fifsf(i-2,j)-Fifsf(i,j)+Fifsf(i,j)*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i-1,j)*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1]);
                    s_b_fi=s_b_fi_num/s_b_fi_denom;
                    s_c_fi=Fifsf(i,j);
                    s_a_fi=Fifsf(i-1,j)/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i,j)/(p->DXP[IM1]*p->DXP[IM1])-s_b_fi/(-1.0*p->DXP[IM1]);

                    p_b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_quad=p_b_quad_num/p_b_quad_denom;
                    p_a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-p_b_quad/p->DXP[IP];
                    
                    p_c_quad=zpos_p;
                    p_p_quad=p_b_quad/p_a_quad;
                    p_q_quad=p_c_quad/p_a_quad;
                    p_x1_quad=(0.0-p_p_quad)/2.0+sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    p_x2_quad=(0.0-p_p_quad)/2.0-sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    
                    p_b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_fi=p_b_fi_num/p_b_fi_denom;
                    p_a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-p_b_fi/p->DXP[IP];
                    p_c_fi=Fifsf(i,j);
                    
                    
                    if(s_x1_quad<=0.0 && s_x1_quad>=0.0-(p->DXP[IM1]) && !(s_x2_quad<=0.0 && s_x2_quad>=0.0-p->DXP[IM1]))
                        s_x_quad=s_x1_quad;
                    else if(s_x2_quad<=0.0 && s_x2_quad>=0.0-(p->DXP[IM1]) && !(s_x1_quad<=0.0 && s_x1_quad>=0.0-p->DXP[IM1]))
                        s_x_quad=s_x2_quad;
                    else if(s_x1_quad <= s_x2_quad+1.0e-06 && s_x1_quad >= s_x2_quad-1.0e-06)
                    {
                        cout<<"equal s s"<< "s x1:"<<s_x1_quad<<"s x2:"<<s_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_ss:"<<zpos_ss<<" DX_s:"<<p->DXP[IM1]<<" DX_ss:"<<p->DXP[IM2]<<endl;
                        s_x_quad=s_x1_quad;
                    }
                    else
                    {   
                        cout<<"mucho problemo s s: "<< "s x1:"<<s_x1_quad<<"s x2:"<<s_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_ss:"<<zpos_ss<<" DX_s:"<<p->DXP[IM1]<<" DX_ss:"<<p->DXP[IM2]<<endl;
                        s_x_quad=-p->DXP[IM1];
                    }


                    if(p_x1_quad<=0.0 && p_x1_quad>=0.0-(p->DXP[IM1]) && !(p_x2_quad<=0.0 && p_x2_quad>=0.0-p->DXP[IM1]))
                        p_x_quad=p_x1_quad;
                    else if(p_x2_quad<=0.0 && p_x2_quad>=0.0-(p->DXP[IM1]) && !(p_x1_quad<=0.0 && p_x1_quad>=0.0-p->DXP[IM1]))
                        p_x_quad=p_x2_quad;
                    else if(p_x1_quad <= p_x2_quad+1.0e-06 && p_x1_quad >= p_x2_quad-1.0e-06)
                    {   
                        cout<<"equal s p"<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=p_x1_quad;
                    }
                    else
                    {   cout<<"mucho problemo s p: "<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=-p->DXP[IM1];
                    }
                    
                    x_quad=(s_x_quad+p_x_quad)/2;
                    
                    teta=(0.0-x_quad)/p->DXP[IM1];
                    
                    s_Fival=s_a_fi*x_quad*x_quad+s_b_fi*x_quad+s_c_fi;
                    p_Fival=p_a_fi*x_quad*x_quad+p_b_fi*x_quad+p_c_fi;
                    Fival=(s_Fival+p_Fival)/2;
                
                    if(p->flag4[Ip1JK]<=AIR)
                    {   
                        teta = teta>1.0e-6?teta:1.0e20;
                
                        a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                        a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.s[n]*Fival;
                        a->M.s[n] = 0.0;
                        
                    }
                    
                     else
                    {
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                        
                    Z_t=p->DXP[IM1];
                    Z_b=p->DXP[IP];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.n[n]+=(M_b_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n]=0.0;
                    }
                    
                    
                }
                
                if(p->A323==79)
                {
                    double zpos_p,zpos_s,zpos_ss,zpos_n;
                    double s_Fival,p_Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double s_a_quad,s_b_quad,s_b_quad_num,s_b_quad_denom,s_c_quad,s_x1_quad,s_x2_quad,s_x_quad,s_p_quad,s_q_quad;
                    double s_a_fi,s_b_fi_denom,s_b_fi_num,s_c_fi,s_b_fi;
                    double p_a_quad,p_b_quad,p_b_quad_num,p_b_quad_denom,p_c_quad,p_x1_quad,p_x2_quad,p_x_quad,p_p_quad,p_q_quad;
                    double p_a_fi,p_b_fi_denom,p_b_fi_num,p_c_fi,p_b_fi;
                    double x_quad, Fival;
                    
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_ss=eta(i-2,j)-(p->ZP[KP]-p->F60);
                        
                    s_b_quad_denom=(-1.0*p->DXP[IM1])+(-1.0*p->DXP[IM2])-((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(-1.0*p->DXP[IM1]);
                    s_b_quad_num=zpos_ss-zpos_p+zpos_p*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1])-zpos_s*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1]);
                    s_b_quad=s_b_quad_num/s_b_quad_denom;
                    s_c_quad=zpos_p;
                    s_a_quad=zpos_s/(p->DXP[IM1]*p->DXP[IM1])-zpos_p/(p->DXP[IM1]*p->DXP[IM1])-s_b_quad/(-1.0*p->DXP[IM1]);
                    s_p_quad=s_b_quad/s_a_quad;
                    s_q_quad=s_c_quad/s_a_quad;
                    s_x1_quad=(0.0-s_p_quad)/2.0+sqrt(((s_p_quad/2.0)*(s_p_quad/2.0)-s_q_quad));
                    s_x2_quad=(0.0-s_p_quad)/2.0-sqrt(((s_p_quad/2.0)*(s_p_quad/2.0)-s_q_quad));
                        
                    s_b_fi_denom=(-1.0*p->DXP[IM1])+(-1.0*p->DXP[IM2])-((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(-1.0*p->DXP[IM1]);
                    s_b_fi_num=Fifsf(i-2,j)-Fifsf(i,j)+Fifsf(i,j)*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i-1,j)*((p->DXP[IM1]+p->DXP[IM2])*(p->DXP[IM1]+p->DXP[IM2]))/(p->DXP[IM1]*p->DXP[IM1]);
                    s_b_fi=s_b_fi_num/s_b_fi_denom;
                    s_c_fi=Fifsf(i,j);
                    s_a_fi=Fifsf(i-1,j)/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i,j)/(p->DXP[IM1]*p->DXP[IM1])-s_b_fi/(-1.0*p->DXP[IM1]);

                    p_b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_quad=p_b_quad_num/p_b_quad_denom;
                    p_a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-p_b_quad/p->DXP[IP];
                    
                    p_c_quad=zpos_p;
                    p_p_quad=p_b_quad/p_a_quad;
                    p_q_quad=p_c_quad/p_a_quad;
                    p_x1_quad=(0.0-p_p_quad)/2.0+sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    p_x2_quad=(0.0-p_p_quad)/2.0-sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    
                    p_b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_fi=p_b_fi_num/p_b_fi_denom;
                    p_a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-p_b_fi/p->DXP[IP];
                    p_c_fi=Fifsf(i,j);
                    
                    
                    if(s_x1_quad<=0.0 && s_x1_quad>=0.0-(p->DXP[IM1]) && !(s_x2_quad<=0.0 && s_x2_quad>=0.0-p->DXP[IM1]))
                        s_x_quad=s_x1_quad;
                    else if(s_x2_quad<=0.0 && s_x2_quad>=0.0-(p->DXP[IM1]) && !(s_x1_quad<=0.0 && s_x1_quad>=0.0-p->DXP[IM1]))
                        s_x_quad=s_x2_quad;
                    else if(s_x1_quad <= s_x2_quad+1.0e-06 && s_x1_quad >= s_x2_quad-1.0e-06)
                    {
                        cout<<"equal s s"<< "s x1:"<<s_x1_quad<<"s x2:"<<s_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_ss:"<<zpos_ss<<" DX_s:"<<p->DXP[IM1]<<" DX_ss:"<<p->DXP[IM2]<<endl;
                        s_x_quad=s_x1_quad;
                    }
                    else
                    {   
                        cout<<"mucho problemo s s: "<< "s x1:"<<s_x1_quad<<"s x2:"<<s_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_ss:"<<zpos_ss<<" DX_s:"<<p->DXP[IM1]<<" DX_ss:"<<p->DXP[IM2]<<endl;
                        s_x_quad=-p->DXP[IM1];
                    }


                    if(p_x1_quad<=0.0 && p_x1_quad>=0.0-(p->DXP[IM1]) && !(p_x2_quad<=0.0 && p_x2_quad>=0.0-p->DXP[IM1]))
                        p_x_quad=p_x1_quad;
                    else if(p_x2_quad<=0.0 && p_x2_quad>=0.0-(p->DXP[IM1]) && !(p_x1_quad<=0.0 && p_x1_quad>=0.0-p->DXP[IM1]))
                        p_x_quad=p_x2_quad;
                    else if(p_x1_quad <= p_x2_quad+1.0e-06 && p_x1_quad >= p_x2_quad-1.0e-06)
                    {   
                        cout<<"equal s p"<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=p_x1_quad;
                    }
                    else
                    {   cout<<"mucho problemo s p: "<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=-p->DXP[IM1];
                    }
                    
                    x_quad=(s_x_quad+p_x_quad)/2;
                    
                    teta=(0.0-x_quad)/p->DXP[IM1];
                    teta = teta>1.0e-6?teta:1.0e20;
                    
                    s_Fival=s_a_fi*x_quad*x_quad+s_b_fi*x_quad+s_c_fi;
                    p_Fival=p_a_fi*x_quad*x_quad+p_b_fi*x_quad+p_c_fi;
                    Fival=(s_Fival+p_Fival)/2;
                
                    a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                    a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.s[n]*Fival;
                    a->M.s[n] = 0.0;
                }
                
                
                if(p->A323==8)
                {
                    double zpos_p,zpos_s, delta_x;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    
                    delta_x=p->DXP[IM1];
                    xpos_zero=((-zpos_p)*(-delta_x))/(zpos_s-zpos_p);
                    teta=xpos_zero/(-delta_x);
                    Fival = teta*Fifsf(i-1,j) + (1.0-teta)*Fifsf(i,j);
                    
                    if(p->flag4[Ip1JK]==AIR)
                    {
                        teta = teta>1.0e-6?teta:1.0e20;
                
                        a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                        a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.s[n]*Fival;
                        a->M.s[n] = 0.0;
                  
                    }
                    
                    else
                    {
                        if(teta<1.0e-7)
                            teta=1.0e-7;
                        Z_t=p->DXP[IM1];
                        Z_b=p->DXP[IP];
                    
                        denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                        M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                        M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                        a->M.p[n]+=(M_p_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.n[n]+=(M_b_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                        a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.s[n]=0.0;
                        
                    }
                }
                
                if(p->A323==9)
                {
                    double zpos_p,zpos_s,zpos_n;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    double x_quad_b,Fival_b;
                    
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    if(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1] && !(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1]))
                        x_quad=x1_quad;
                    else if(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1] && !(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=-0.5*p->DXP[IM1];
                    }
                    
                    teta=(0.0-x_quad)/p->DXP[IM1];
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                        
                    if(p->flag4[Ip1JK]==AIR)
                    {   
                        //if(p->flag4[IJKp1]==AIR)
                         //   cout<<"Air s LRT x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                       // else
                        //    cout<<"AIR s LR!!! x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        
                        if(x_quad==x1_quad && (x2_quad>=0.0 && x2_quad<=p->DXP[IP]))
                            x_quad_b=x2_quad;
                        else if(x_quad==x2_quad && (x1_quad>=0.0 && x1_quad<=p->DXP[IP]))
                            x_quad_b=x1_quad;
                        else if (x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                        {
                            cout<<"equal_b s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                            x_quad_b=x1_quad;
                        }
                        else
                        {
                            cout<<"tja s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                            x_quad_b=0.5*p->DXP[IP];
                        }
                        
                        Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                        Fival_b=a_fi*x_quad_b*x_quad_b+b_fi*x_quad_b+c_fi;
                        
                        Z_t=p->DXP[IM1];
                        Z_b=x_quad_b;
                    
                        denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                        M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                        M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                        a->M.p[n]+=(M_p_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                        a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IM1]*p->DXN[IP]);
                        a->rhsvec.V[n]-=(Fival_b*M_b_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.n[n]=0.0;
                        a->M.s[n]=0.0;
                        
                    }
                    
                    else
                    {
                    
                    Fival = a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                    
                    Z_t=p->DXP[IM1];
                    Z_b=p->DXP[IP];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.n[n]+=(M_b_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n]=0.0;
                    }
                }
                
                if(p->A323==10)
                {
                    double zpos_p,zpos_s,zpos_n;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    double x_quad_b,Fival_b;
                    
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    if(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1] && !(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1]))
                        x_quad=x1_quad;
                    else if(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1] && !(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=-0.5*p->DXP[IM1];
                    }
                    
                    teta=(0.0-x_quad)/p->DXP[IM1];
                    
                    if(p->flag4[Ip1JK]==AIR)
                    {   
                        //if(p->flag4[IJKp1]==AIR)
                        //    cout<<"Air s LRT x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                       // else
                       //     cout<<"AIR s LR!!! x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        
                        Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                        
                        teta = teta>1.0e-6?teta:1.0e20;
                
                        a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                        a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.s[n]*Fival;
                        a->M.s[n] = 0.0;
                        
                    }
                    
                     else
                    {
                    
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                        
                    Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                    
                    Z_t=p->DXP[IM1];
                    Z_b=p->DXP[IP];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.n[n]+=(M_b_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n]=0.0;
                    }
                }
                
                if(p->A323==11)
                {
                    double zpos_p,zpos_s,zpos_n;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad,x_defin;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    
                    double coeff_l,coeff_r,zpos_ss;
                    double b_plus_denom,b_plus_num,b_plus,a_plus,c_plus,q_plus,p_plus,x1_plus,x2_plus,x_plus;
                    double b_fi_plus_denom,b_fi_plus_num,b_fi_plus,a_fi_plus,c_fi_plus,Fi_plus,Fi_quad;
                    
                    coeff_l=0.5;
                    coeff_r=0.5;
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_ss=eta(i-2,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    b_plus_denom=p->DXP[IM2]*p->DXP[IM2]/p->DXP[IM1]+p->DXP[IM2];
                    b_plus_num=zpos_p*p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-zpos_ss-zpos_s*(p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-1.0);
                    b_plus=b_plus_num/b_plus_denom;
                    a_plus=zpos_p/(p->DXP[IM1]*p->DXP[IM1])-zpos_s/(p->DXP[IM1]*p->DXP[IM1])-b_quad/p->DXP[IM1];
                    c_plus=zpos_s;
                    p_plus=b_plus/a_plus;
                    q_plus=c_plus/a_plus;
                    x1_plus=(0.0-p_plus)/2.0+sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    x2_plus=(0.0-p_plus)/2.0-sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    
                    b_fi_plus_denom=p->DXP[IM2]*p->DXP[IM2]/p->DXP[IM1]+p->DXP[IM2];
                    b_fi_plus_num=Fifsf(i,j)*p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i-2,j)-Fifsf(i-1,j)*(p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-1.0);
                    b_fi_plus=b_fi_plus_num/b_fi_plus_denom;
                    a_fi_plus=Fifsf(i,j)/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i-1,j)/(p->DXP[IM1]*p->DXP[IM1])-b_fi_plus/p->DXP[IM1];
                    c_fi=Fifsf(i-1,j);
                    
                    if(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1] && !(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1]))
                        x_quad=x1_quad;
                    else if(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1] && !(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=-0.5*p->DXP[IM1];
                    }
                    
                    if(x1_plus<=0.0 && x1_plus>=0.0-p->DXP[IM1] && !(x2_plus<=0.0 && x2_plus>=0.0-p->DXP[IM1]))
                        x_plus=x1_plus;
                    else if(x2_plus<=0.0 && x2_plus>=0.0-p->DXP[IM1] && !(x1_plus<=0.0 && x1_plus>=0.0-p->DXP[IM1]))
                        x_plus=x2_plus;
                    else if(x1_plus <= x2_plus+1.0e-6 && x1_plus >= x2_plus-1.0e-6)
                    {
                        cout<<"equal plus s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_plus=x1_plus;
                    }
                    else
                    {
                        cout<<"mucho problemo plus s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_plus=-0.5*p->DXP[IM1];
                    }
                    
                    x_defin=coeff_l*x_plus+coeff_r*x_quad;
                    
                    teta=(0.0-x_defin)/p->DXP[IM1];
                    teta = teta>1.0e-6?teta:1.0e20;
                    Fi_quad=a_fi*x_defin*x_defin+b_fi*x_defin+c_fi;
                    Fi_plus=a_fi_plus*x_defin*x_defin+b_fi_plus*x_defin+c_fi_plus;
                    Fival=coeff_l*Fi_plus+coeff_r*Fi_quad;
                
                    a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                    a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.s[n]*Fival;
                    a->M.s[n] = 0.0;
                }
                
                if(p->A323==12)
                {
                    double zpos_p,zpos_s,zpos_n;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad,x_defin;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    
                    double coeff_l,coeff_r,zpos_ss;
                    double b_plus_denom,b_plus_num,b_plus,a_plus,c_plus,q_plus,p_plus,x1_plus,x2_plus,x_plus;
                    double b_fi_plus_denom,b_fi_plus_num,b_fi_plus,a_fi_plus,c_fi_plus,Fi_plus,Fi_quad;
                    
                    coeff_l=0.5;
                    coeff_r=0.5;
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_ss=eta(i-2,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    b_plus_denom=p->DXP[IM2]*p->DXP[IM2]/p->DXP[IM1]+p->DXP[IM2];
                    b_plus_num=zpos_p*p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-zpos_ss-zpos_s*(p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-1.0);
                    b_plus=b_plus_num/b_plus_denom;
                    a_plus=zpos_p/(p->DXP[IM1]*p->DXP[IM1])-zpos_s/(p->DXP[IM1]*p->DXP[IM1])-b_quad/p->DXP[IM1];
                    c_plus=zpos_s;
                    p_plus=b_plus/a_plus;
                    q_plus=c_plus/a_plus;
                    x1_plus=(0.0-p_plus)/2.0+sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    x2_plus=(0.0-p_plus)/2.0-sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    
                    b_fi_plus_denom=p->DXP[IM2]*p->DXP[IM2]/p->DXP[IM1]+p->DXP[IM2];
                    b_fi_plus_num=Fifsf(i,j)*p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i-2,j)-Fifsf(i-1,j)*(p->DXP[IM2]*p->DXP[IM2]/(p->DXP[IM1]*p->DXP[IM1])-1.0);
                    b_fi_plus=b_fi_plus_num/b_fi_plus_denom;
                    a_fi_plus=Fifsf(i,j)/(p->DXP[IM1]*p->DXP[IM1])-Fifsf(i-1,j)/(p->DXP[IM1]*p->DXP[IM1])-b_fi_plus/p->DXP[IM1];
                    c_fi=Fifsf(i-1,j);
                    
                    if(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1] && !(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1]))
                        x_quad=x1_quad;
                    else if(x2_quad<=0.0 && x2_quad>=0.0-p->DXP[IM1] && !(x1_quad<=0.0 && x1_quad>=0.0-p->DXP[IM1]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=-0.5*p->DXP[IM1];
                    }
                    
                    if(x1_plus<=0.0 && x1_plus>=0.0-p->DXP[IM1] && !(x2_plus<=0.0 && x2_plus>=0.0-p->DXP[IM1]))
                        x_plus=x1_plus;
                    else if(x2_plus<=0.0 && x2_plus>=0.0-p->DXP[IM1] && !(x1_plus<=0.0 && x1_plus>=0.0-p->DXP[IM1]))
                        x_plus=x2_plus;
                    else if(x1_plus <= x2_plus+1.0e-6 && x1_plus >= x2_plus-1.0e-6)
                    {
                        cout<<"equal plus s"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_plus=x1_plus;
                    }
                    else
                    {
                        cout<<"mucho problemo plus s: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_plus=-0.5*p->DXP[IM1];
                    }
                    
                    x_defin=coeff_l*x_plus+coeff_r*x_quad;
                    
                    teta=(0.0-x_defin)/p->DXP[IM1];
                    Fi_quad=a_fi*x_defin*x_defin+b_fi*x_defin+c_fi;
                    Fi_plus=a_fi_plus*x_defin*x_defin+b_fi_plus*x_defin+c_fi_plus;
                    Fival=coeff_l*Fi_plus+coeff_r*Fi_quad;
                    
                    if(p->flag4[Ip1JK]==AIR)
                    {   
                        Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                        
                        teta = teta>1.0e-6?teta:1.0e20;
                
                        a->M.p[n] -= 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                           
                        a->M.s[n] += 1.0/(p->DXP[IM1]*p->DXN[IP]);
                        a->M.s[n] -= 1.0/(teta*p->DXP[IM1]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.s[n]*Fival;
                        a->M.s[n] = 0.0;
                        
                    }
                    
                     else
                    {
                    
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                        
                    Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                    
                    Z_t=p->DXP[IM1];
                    Z_b=p->DXP[IP];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.n[n]+=(M_b_num/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IM1]*p->DXN[IP]);
                    a->M.s[n]=0.0;
                    }
                }
            
            }

            // north
            if(p->flag4[Ip1JK]==AIR)
            {
                // -----------
                if(p->A323==1)
                {
                a->rhsvec.V[n] -= a->M.n[n]*f(i+1,j,k);
                a->M.n[n] = 0.0;
                }
                
                // -----------
                if(p->A323==2)
                {
                double lsv0,lsv1,Fival;

                lsv0 = fabs(a->phi(i,j,k));
                lsv1 = fabs(a->phi(i+1,j,k));

                lsv0 = fabs(lsv0)>1.0e-6?lsv0:1.0e20;
                
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k)));
                Fival = teta*Fifsf(i+1,j) + (1.0-teta)*Fifsf(i,j);


                a->rhsvec.V[n] -= a->M.n[n]*Fival*(1.0 + lsv1/lsv0);
                a->M.p[n] -= a->M.n[n]*lsv1/lsv0;
                a->M.n[n] = 0.0;
                }
                
                
                // -----------
                if(p->A323==3)
                {
                double x0,x1,x2,y2;
                double x,y,Fival;
                double Lx0,Lx1,Lx2;
                
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k)));
                Fival = teta*Fifsf(i+1,j) + (1.0-teta)*Fifsf(i,j);

                x0 = -fabs(a->phi(i-1,j,k));
                x1 = -fabs(a->phi(i,j,k));
                x2 = 0.0;

                y2 = Fival;

                x = fabs(a->phi(i+1,j,k));

                Lx0 = ((x-x1)/(x0-x1)) * ((x-x2)/(x0-x2));
                Lx1 = ((x-x0)/(x1-x0)) * ((x-x2)/(x1-x2));
                Lx2 = ((x-x0)/(x2-x0)) * ((x-x1)/(x2-x1));

                a->rhsvec.V[n]  -= a->M.n[n]*Lx2*y2;
                a->M.p[n]       += a->M.n[n]*Lx1;
                a->M.s[n]       += a->M.n[n]*Lx0;
                a->M.n[n]       = 0.0;
                }
                
                // -----------
                if(p->A323==4)
                {
                double Fival,teta;

                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k)))+ 0.0001*p->DXN[IP]/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k)));
              
                Fival = teta*Fifsf(i+1,j) + (1.0-teta)*Fifsf(i,j);
                
                teta = teta>1.0e-6?teta:1.0e20;
            
                a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                a->rhsvec.V[n] -= a->M.n[n]*Fival;
                a->M.n[n] = 0.0;
                }
                
                if(p->A323==5)
                {
                    double Fival,teta;
                    double zpos_p,zpos_n,x_pos;
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    x_pos=((0.0-zpos_p)*p->DXP[IP])/(zpos_n-zpos_p);
                    
                    teta = x_pos/p->DXP[IP];
                    Fival = teta*Fifsf(i+1,j) + (1.0-teta)*Fifsf(i,j);
                
                    teta = teta>1.0e-6?teta:1.0e20;
            
                    a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.n[n]*Fival;
                    a->M.n[n] = 0.0;
                  
                }
                
                if(p->A323==6 || p->A323==7 || p->A323==77)
                {
                    double zpos_p,zpos_n,zpos_s;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_denom,b_fi_num,c_fi,b_fi;
                    double x_quad_b,Fival_b;
                    
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    if(x1_quad>=0.0 && x1_quad<=p->DXP[IP] && !(x2_quad>=0.0 && x2_quad<=p->DXP[IP]))
                        x_quad=x1_quad;
                    else if(x2_quad>=0.0 && x2_quad<=p->DXP[IP] && !(x1_quad>=0.0 && x1_quad<=p->DXP[IP]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=p->DXP[IP];
                    }
                    
                    teta=(x_quad)/p->DXP[IP];
                    Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                    
                    teta = teta>1.0e-6?teta:1.0e20;
            
                    a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.n[n]*Fival;
                    a->M.n[n] = 0.0;
                    
                }
                
                if(p->A323==78)
                {
                    
                    double zpos_p,zpos_n,zpos_nn,zpos_s;
                    double n_Fival,p_Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double n_a_quad,n_b_quad,n_b_quad_num,n_b_quad_denom,n_c_quad,n_x1_quad,n_x2_quad,n_x_quad,n_p_quad,n_q_quad;
                    double n_a_fi,n_b_fi_denom,n_b_fi_num,n_c_fi,n_b_fi;
                    double p_a_quad,p_b_quad,p_b_quad_num,p_b_quad_denom,p_c_quad,p_x1_quad,p_x2_quad,p_x_quad,p_p_quad,p_q_quad;
                    double p_a_fi,p_b_fi_denom,p_b_fi_num,p_c_fi,p_b_fi;
                    double x_quad, Fival;
                    
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    zpos_nn=eta(i+2,j)-(p->ZP[KP]-p->F60);
                        
                    n_b_quad_denom=p->DXP[IP]+p->DXP[IP1]-((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/p->DXP[IP];
                    n_b_quad_num=zpos_nn-zpos_p+zpos_p*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP])-zpos_n*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP]);
                    n_b_quad=n_b_quad_num/n_b_quad_denom;
                    n_c_quad=zpos_p;
                    n_a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-n_b_quad/p->DXP[IP];
                    n_p_quad=n_b_quad/n_a_quad;
                    n_q_quad=n_c_quad/n_a_quad;
                    n_x1_quad=(0.0-n_p_quad)/2.0+sqrt(((n_p_quad/2.0)*(n_p_quad/2.0)-n_q_quad));
                    n_x2_quad=(0.0-n_p_quad)/2.0-sqrt(((n_p_quad/2.0)*(n_p_quad/2.0)-n_q_quad));
                        
                    n_b_fi_denom=p->DXP[IP]+p->DXP[IP1]-((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/p->DXP[IP];
                    n_b_fi_num=Fifsf(i+2,j)-Fifsf(i,j)+Fifsf(i,j)*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP])-Fifsf(i+1,j)*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP]);
                    n_b_fi=n_b_fi_num/n_b_fi_denom;
                    n_c_fi=Fifsf(i,j);
                    n_a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-n_b_fi/p->DXP[IP];
                    
                    p_b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_quad=p_b_quad_num/p_b_quad_denom;
                    p_a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-p_b_quad/p->DXP[IP];
                    p_c_quad=zpos_p;
                    p_p_quad=p_b_quad/p_a_quad;
                    p_q_quad=p_c_quad/p_a_quad;
                    p_x1_quad=(0.0-p_p_quad)/2.0+sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    p_x2_quad=(0.0-p_p_quad)/2.0-sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    
                    p_b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_fi=p_b_fi_num/p_b_fi_denom;
                    p_a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-p_b_fi/p->DXP[IP];
                    p_c_fi=Fifsf(i,j);
                    
                    
                    if(n_x1_quad>=0.0 && n_x1_quad<=p->DXP[IP] && !(n_x2_quad>=0.0 && n_x2_quad<=p->DXP[IP]))
                        n_x_quad=n_x1_quad;
                    else if(n_x2_quad>=0.0 && n_x2_quad<=p->DXP[IP] && !(n_x1_quad>=0.0 && n_x1_quad<=p->DXP[IP]))
                        n_x_quad=n_x2_quad;
                    else if(n_x1_quad <= n_x2_quad+1.0e-06 && n_x1_quad >= n_x2_quad-1.0e-06)
                    {   
                        cout<<"equal n n"<< "n x1:"<<n_x1_quad<<"n x2:"<<n_x2_quad<<" zpos_nn:" << zpos_nn<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_nn:"<<p->DXP[IP1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        n_x_quad=n_x1_quad;
                    }
                    else
                    {   
                        cout<<"mucho problemo n n: "<< "n x1:"<<n_x1_quad<<"n x2:"<<n_x2_quad<<" zpos_nn:" << zpos_nn<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_nn:"<<p->DXP[IP1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        n_x_quad=p->DXP[IP];
                    }
                    
                    if(p_x1_quad>=0.0 && p_x1_quad<=p->DXP[IP] && !(p_x2_quad>=0.0 && p_x2_quad<=p->DXP[IP]))
                        p_x_quad=p_x1_quad;
                    else if(p_x2_quad>=0.0 && p_x2_quad<=p->DXP[IP] && !(p_x1_quad>=0.0 && p_x1_quad<=p->DXP[IP]))
                        p_x_quad=p_x2_quad;
                    else if(p_x1_quad <= p_x2_quad+1.0e-06 && p_x1_quad >= p_x2_quad-1.0e-06)
                    {   
                        cout<<"equal n p"<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=p_x1_quad;
                    }
                    else
                    {   
                        cout<<"mucho problemo n p: "<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=p->DXP[IP];
                    }
                    
                    x_quad=(p_x_quad+n_x_quad)/2;
                    teta=(x_quad)/p->DXP[IP];
                    
                    n_Fival=n_a_fi*x_quad*x_quad+n_b_fi*x_quad+n_c_fi;
                    p_Fival=p_a_fi*x_quad*x_quad+p_b_fi*x_quad+p_c_fi;
                    Fival=(n_Fival+p_Fival)/2;
                    
                    if(p->flag4[Im1JK]<=AIR)
                    {
                       
                        teta = teta>1.0e-6?teta:1.0e20;
            
                        a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.n[n]*Fival;
                        a->M.n[n] = 0.0;
                        
                        
                    }
                    
                    else
                    {
                        
                    if(teta<1.0e-7)
                        teta=1.0e-7;
        
                    Z_t=p->DXP[IP];
                    Z_b=p->DXP[IM1];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.s[n]+=(M_b_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n]=0.0;
                    }
                }
                
                if(p->A323==79)
                {
                    
                    double zpos_p,zpos_n,zpos_nn,zpos_s;
                    double n_Fival,p_Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double n_a_quad,n_b_quad,n_b_quad_num,n_b_quad_denom,n_c_quad,n_x1_quad,n_x2_quad,n_x_quad,n_p_quad,n_q_quad;
                    double n_a_fi,n_b_fi_denom,n_b_fi_num,n_c_fi,n_b_fi;
                    double p_a_quad,p_b_quad,p_b_quad_num,p_b_quad_denom,p_c_quad,p_x1_quad,p_x2_quad,p_x_quad,p_p_quad,p_q_quad;
                    double p_a_fi,p_b_fi_denom,p_b_fi_num,p_c_fi,p_b_fi;
                    double x_quad, Fival;
                    
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    zpos_nn=eta(i+2,j)-(p->ZP[KP]-p->F60);
                        
                    n_b_quad_denom=p->DXP[IP]+p->DXP[IP1]-((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/p->DXP[IP];
                    n_b_quad_num=zpos_nn-zpos_p+zpos_p*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP])-zpos_n*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP]);
                    n_b_quad=n_b_quad_num/n_b_quad_denom;
                    n_c_quad=zpos_p;
                    n_a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-n_b_quad/p->DXP[IP];
                    n_p_quad=n_b_quad/n_a_quad;
                    n_q_quad=n_c_quad/n_a_quad;
                    n_x1_quad=(0.0-n_p_quad)/2.0+sqrt(((n_p_quad/2.0)*(n_p_quad/2.0)-n_q_quad));
                    n_x2_quad=(0.0-n_p_quad)/2.0-sqrt(((n_p_quad/2.0)*(n_p_quad/2.0)-n_q_quad));
                        
                    n_b_fi_denom=p->DXP[IP]+p->DXP[IP1]-((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/p->DXP[IP];
                    n_b_fi_num=Fifsf(i+2,j)-Fifsf(i,j)+Fifsf(i,j)*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP])-Fifsf(i+1,j)*((p->DXP[IP]+p->DXP[IP1])*(p->DXP[IP]+p->DXP[IP1]))/(p->DXP[IP]*p->DXP[IP]);
                    n_b_fi=n_b_fi_num/n_b_fi_denom;
                    n_c_fi=Fifsf(i,j);
                    n_a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-n_b_fi/p->DXP[IP];
                    
                    p_b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_quad=p_b_quad_num/p_b_quad_denom;
                    p_a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-p_b_quad/p->DXP[IP];
                    p_c_quad=zpos_p;
                    p_p_quad=p_b_quad/p_a_quad;
                    p_q_quad=p_c_quad/p_a_quad;
                    p_x1_quad=(0.0-p_p_quad)/2.0+sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    p_x2_quad=(0.0-p_p_quad)/2.0-sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    
                    p_b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    p_b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    p_b_fi=p_b_fi_num/p_b_fi_denom;
                    p_a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-p_b_fi/p->DXP[IP];
                    p_c_fi=Fifsf(i,j);
                    
                    
                    if(n_x1_quad>=0.0 && n_x1_quad<=p->DXP[IP] && !(n_x2_quad>=0.0 && n_x2_quad<=p->DXP[IP]))
                        n_x_quad=n_x1_quad;
                    else if(n_x2_quad>=0.0 && n_x2_quad<=p->DXP[IP] && !(n_x1_quad>=0.0 && n_x1_quad<=p->DXP[IP]))
                        n_x_quad=n_x2_quad;
                    else if(n_x1_quad <= n_x2_quad+1.0e-06 && n_x1_quad >= n_x2_quad-1.0e-06)
                    {   
                        cout<<"equal n n"<< "n x1:"<<n_x1_quad<<"n x2:"<<n_x2_quad<<" zpos_nn:" << zpos_nn<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_nn:"<<p->DXP[IP1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        n_x_quad=n_x1_quad;
                    }
                    else
                    {   
                        cout<<"mucho problemo n n: "<< "n x1:"<<n_x1_quad<<"n x2:"<<n_x2_quad<<" zpos_nn:" << zpos_nn<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_nn:"<<p->DXP[IP1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        n_x_quad=p->DXP[IP];
                    }
                    
                    if(p_x1_quad>=0.0 && p_x1_quad<=p->DXP[IP] && !(p_x2_quad>=0.0 && p_x2_quad<=p->DXP[IP]))
                        p_x_quad=p_x1_quad;
                    else if(p_x2_quad>=0.0 && p_x2_quad<=p->DXP[IP] && !(p_x1_quad>=0.0 && p_x1_quad<=p->DXP[IP]))
                        p_x_quad=p_x2_quad;
                    else if(p_x1_quad <= p_x2_quad+1.0e-06 && p_x1_quad >= p_x2_quad-1.0e-06)
                    {   
                        cout<<"equal n p"<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=p_x1_quad;
                    }
                    else
                    {   
                        cout<<"mucho problemo n p: "<< "p x1:"<<p_x1_quad<<"p x2:"<<p_x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        p_x_quad=p->DXP[IP];
                    }
                    
                    x_quad=(p_x_quad+n_x_quad)/2;
                    teta=(x_quad)/p->DXP[IP];
                    
                    n_Fival=n_a_fi*x_quad*x_quad+n_b_fi*x_quad+n_c_fi;
                    p_Fival=p_a_fi*x_quad*x_quad+p_b_fi*x_quad+p_c_fi;
                    Fival=(n_Fival+p_Fival)/2;
                    
                    teta = teta>1.0e-6?teta:1.0e20;
            
                    a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.n[n]*Fival;
                    a->M.n[n] = 0.0;
                }
            
                if(p->A323==8)
                {
                    double zpos_p,zpos_n, delta_x;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    delta_x=p->DXP[IP];
                    xpos_zero=((-zpos_p)*delta_x)/(zpos_n-zpos_p);
                    teta=xpos_zero/delta_x;
                    Fival=teta*Fifsf(i+1,j)+(1.0-teta)*Fifsf(i,j);
                    
                    if(p->flag4[Im1JK]==AIR)
                    {
                        teta = teta>1.0e-6?teta:1.0e20;
            
                        a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.n[n]*Fival;
                        a->M.n[n] = 0.0;
                    }
                    
                    else
                    {
                        if(teta<1.0e-7)
                            teta=1.0e-7;
                        Z_t=p->DXP[IP];
                        Z_b=p->DXP[IM1];
                    
                        denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                        M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                        M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                        a->M.p[n]+=(M_p_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                        a->M.s[n]+=(M_b_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                        a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IP]*p->DXN[IP]);
                        a->M.n[n]=0.0;
                    }
                }
                
                if(p->A323==9)
                {
                    double zpos_p,zpos_n,zpos_s;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    double x_quad_b,Fival_b;
                    
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    if(x1_quad>=0.0 && x1_quad<=p->DXP[IP] && !(x2_quad>=0.0 && x2_quad<=p->DXP[IP]))
                        x_quad=x1_quad;
                    else if(x2_quad>=0.0 && x2_quad<=p->DXP[IP] && !(x1_quad>=0.0 && x1_quad<=p->DXP[IP]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=0.5*p->DXP[IP];
                    }
                    
                    teta=(x_quad)/p->DXP[IP];
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                        
                    if(p->flag4[Im1JK]==AIR)
                    {
                        //if(p->flag4[IJKp1]==AIR)
                        //    cout<<"Air n LRT x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                       // else
                        //    cout<<"AIR n LR!!! x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        
                        if(x_quad==x1_quad && (x2_quad<=0.0 && x2_quad>=(0.0-p->DXP[IM1])))
                            x_quad_b=x2_quad;
                        else if(x_quad==x2_quad && (x1_quad<=0.0 && x1_quad>=(0.0-p->DXP[IM1])))
                            x_quad_b=x1_quad;
                        else if (x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                        {
                            cout<<"equal_b n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                            x_quad_b=x1_quad;
                        }
                        else
                        {
                            cout<<"tja n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                            x_quad_b=-0.5*p->DXP[IM1];
                        }
                        
                        Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                        Fival_b=a_fi*x_quad_b*x_quad_b+b_fi*x_quad_b+c_fi;
                        
                        Z_t=p->DXP[IP];
                        Z_b=(0.0-x_quad_b);
                    
                        denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                        M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                        M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                        a->M.p[n]+=(M_p_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                        a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IP]*p->DXN[IP]);
                        a->rhsvec.V[n]-=(Fival_b*M_b_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                        a->M.n[n]=0.0;
                        a->M.s[n]=0.0;
                        
                    }
                    
                    else
                    {
                    
                    Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                    
                    Z_t=p->DXP[IP];
                    Z_b=p->DXP[IM1];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.s[n]+=(M_b_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n]=0.0;
                    }
                }
                
                if(p->A323==10)
                {
                    double zpos_p,zpos_n,zpos_s;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    double x_quad_b,Fival_b;
                    
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    if(x1_quad>=0.0 && x1_quad<=p->DXP[IP] && !(x2_quad>=0.0 && x2_quad<=p->DXP[IP]))
                        x_quad=x1_quad;
                    else if(x2_quad>=0.0 && x2_quad<=p->DXP[IP] && !(x1_quad>=0.0 && x1_quad<=p->DXP[IP]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=0.5*p->DXP[IP];
                    }
                    
                    teta=(x_quad)/p->DXP[IP];
                        
                    if(p->flag4[Im1JK]==AIR)
                    {
                       // if(p->flag4[IJKp1]==AIR)
                         //   cout<<"Air n LRT x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                       // else
                       //     cout<<"AIR n LR!!! x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
            
                       Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                        
                        teta = teta>1.0e-6?teta:1.0e20;
            
                        a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.n[n]*Fival;
                        a->M.n[n] = 0.0;
                        
                        
                    }
                    
                    else
                    {
                        
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                    Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                    Z_t=p->DXP[IP];
                    Z_b=p->DXP[IM1];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.s[n]+=(M_b_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n]=0.0;
                    }
                }
                
                if(p->A323==11)
                {
                    double zpos_p,zpos_n,zpos_s;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    double x_quad_b,Fival_b;
                    
                    double coeff_l,coeff_r,zpos_nn;
                    double b_plus_denom,b_plus_num,b_plus,a_plus,c_plus,q_plus,p_plus,x1_plus,x2_plus,x_plus,x_defin;
                    double b_fi_plus_denom,b_fi_plus_num,b_fi_plus,a_fi_plus,c_fi_plus,Fi_plus,Fi_quad;
                    
                    coeff_l=0.5;
                    coeff_r=0.5;
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    zpos_nn=eta(i+2,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    b_plus_denom=p->DXP[IP]*p->DXP[IP]/p->DXP[IP1]+p->DXP[IP];
                    b_plus_num=zpos_nn*p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-zpos_p-zpos_n*(p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-1.0);
                    b_plus=b_plus_num/b_plus_denom;
                    a_plus=zpos_nn/(p->DXP[IP1]*p->DXP[IP1])-zpos_n/(p->DXP[IP1]*p->DXP[IP1])-b_plus/p->DXP[IP1];
                    c_plus=zpos_n;
                    p_plus=b_plus/a_plus;
                    q_plus=c_plus/a_plus;
                    x1_plus=(0.0-p_plus)/2.0+sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    x2_plus=(0.0-p_plus)/2.0-sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    
                    b_fi_plus_denom=p->DXP[IP]*p->DXP[IP]/p->DXP[IP1]+p->DXP[IP];
                    b_fi_plus_num=Fifsf(i+2,j)*p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-Fifsf(i,j)-Fifsf(i+1,j)*(p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-1.0);
                    b_fi_plus=b_fi_plus_num/b_fi_plus_denom;
                    a_fi_plus=Fifsf(i+2,j)/(p->DXP[IP1]*p->DXP[IP1])-Fifsf(i+1,j)/(p->DXP[IP1]*p->DXP[IP1])-b_fi_plus/p->DXP[IP1];
                    c_fi=Fifsf(i+1,j);
                    
                    if(x1_quad>=0.0 && x1_quad<=p->DXP[IP] && !(x2_quad>=0.0 && x2_quad<=p->DXP[IP]))
                        x_quad=x1_quad;
                    else if(x2_quad>=0.0 && x2_quad<=p->DXP[IP] && !(x1_quad>=0.0 && x1_quad<=p->DXP[IP]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=0.5*p->DXP[IP];
                    }
                    
                    if(x1_plus>=0.0 && x1_plus<=p->DXP[IP] && !(x2_plus>=0.0 && x2_plus<=p->DXP[IP]))
                        x_plus=x1_plus;
                    else if(x2_plus>=0.0 && x2_plus<=p->DXP[IP] && !(x1_plus>=0.0 && x1_plus<=p->DXP[IP]))
                        x_plus=x2_plus;
                    else if(x1_plus <= x2_plus+1.0e-06 && x1_plus >= x2_plus-1.0e-06)
                    {
                        cout<<"equal plus n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;;
                        x_plus=x1_plus;
                    }
                    else
                    {
                        cout<<"mucho problemo n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_plus=0.5*p->DXP[IP];
                    }
                    
                    x_defin=coeff_l*x_quad+coeff_r*x_plus;
                    
                    teta=(0.0-x_defin)/p->DXP[IM1];
                    teta = teta>1.0e-6?teta:1.0e20;
                    Fi_quad=a_fi*x_defin*x_defin+b_fi*x_defin+c_fi;
                    Fi_plus=a_fi_plus*x_defin*x_defin+b_fi_plus*x_defin+c_fi_plus;
                    Fival=coeff_l*Fi_quad+coeff_r*Fi_plus;
                    
                    teta = teta>1.0e-6?teta:1.0e20;
            
                    a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                    a->rhsvec.V[n] -= a->M.n[n]*Fival;
                    a->M.n[n] = 0.0;
                }
                
                if(p->A323==12)
                {
                    double zpos_p,zpos_n,zpos_s;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,x1_quad,x2_quad,x_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    double x_quad_b,Fival_b;
                    
                    double coeff_l,coeff_r,zpos_nn;
                    double b_plus_denom,b_plus_num,b_plus,a_plus,c_plus,q_plus,p_plus,x1_plus,x2_plus,x_plus,x_defin;
                    double b_fi_plus_denom,b_fi_plus_num,b_fi_plus,a_fi_plus,c_fi_plus,Fi_plus,Fi_quad;
                    
                    coeff_l=0.5;
                    coeff_r=0.5;
                    
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_s=eta(i-1,j)-(p->ZP[KP]-p->F60);
                    zpos_n=eta(i+1,j)-(p->ZP[KP]-p->F60);
                    zpos_nn=eta(i+2,j)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_quad_num=zpos_n*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-zpos_s-zpos_p*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_n/(p->DXP[IP]*p->DXP[IP])-zpos_p/(p->DXP[IP]*p->DXP[IP])-b_quad/p->DXP[IP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    x1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    x2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DXP[IM1]*p->DXP[IM1]/p->DXP[IP]+p->DXP[IM1];
                    b_fi_num=Fifsf(i+1,j)*p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-Fifsf(i-1,j)-Fifsf(i,j)*(p->DXP[IM1]*p->DXP[IM1]/(p->DXP[IP]*p->DXP[IP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i+1,j)/(p->DXP[IP]*p->DXP[IP])-Fifsf(i,j)/(p->DXP[IP]*p->DXP[IP])-b_fi/p->DXP[IP];
                    c_fi=Fifsf(i,j);
                    
                    b_plus_denom=p->DXP[IP]*p->DXP[IP]/p->DXP[IP1]+p->DXP[IP];
                    b_plus_num=zpos_nn*p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-zpos_p-zpos_n*(p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-1.0);
                    b_plus=b_plus_num/b_plus_denom;
                    a_plus=zpos_nn/(p->DXP[IP1]*p->DXP[IP1])-zpos_n/(p->DXP[IP1]*p->DXP[IP1])-b_plus/p->DXP[IP1];
                    c_plus=zpos_n;
                    p_plus=b_plus/a_plus;
                    q_plus=c_plus/a_plus;
                    x1_plus=(0.0-p_plus)/2.0+sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    x2_plus=(0.0-p_plus)/2.0-sqrt(((p_plus/2.0)*(p_plus/2.0)-q_plus));
                    
                    b_fi_plus_denom=p->DXP[IP]*p->DXP[IP]/p->DXP[IP1]+p->DXP[IP];
                    b_fi_plus_num=Fifsf(i+2,j)*p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-Fifsf(i,j)-Fifsf(i+1,j)*(p->DXP[IP]*p->DXP[IP]/(p->DXP[IP1]*p->DXP[IP1])-1.0);
                    b_fi_plus=b_fi_plus_num/b_fi_plus_denom;
                    a_fi_plus=Fifsf(i+2,j)/(p->DXP[IP1]*p->DXP[IP1])-Fifsf(i+1,j)/(p->DXP[IP1]*p->DXP[IP1])-b_fi_plus/p->DXP[IP1];
                    c_fi=Fifsf(i+1,j);
                    
                    if(x1_quad>=0.0 && x1_quad<=p->DXP[IP] && !(x2_quad>=0.0 && x2_quad<=p->DXP[IP]))
                        x_quad=x1_quad;
                    else if(x2_quad>=0.0 && x2_quad<=p->DXP[IP] && !(x1_quad>=0.0 && x1_quad<=p->DXP[IP]))
                        x_quad=x2_quad;
                    else if(x1_quad <= x2_quad+1.0e-06 && x1_quad >= x2_quad-1.0e-06)
                    {
                        cout<<"equal n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;;
                        x_quad=x1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_quad=0.5*p->DXP[IP];
                    }
                    
                    if(x1_plus>=0.0 && x1_plus<=p->DXP[IP] && !(x2_plus>=0.0 && x2_plus<=p->DXP[IP]))
                        x_plus=x1_plus;
                    else if(x2_plus>=0.0 && x2_plus<=p->DXP[IP] && !(x1_plus>=0.0 && x1_plus<=p->DXP[IP]))
                        x_plus=x2_plus;
                    else if(x1_plus <= x2_plus+1.0e-06 && x1_plus >= x2_plus-1.0e-06)
                    {
                        cout<<"equal plus n"<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;;
                        x_plus=x1_plus;
                    }
                    else
                    {
                        cout<<"mucho problemo plus n: "<< " x1:"<<x1_quad<<" x2:"<<x2_quad<<" zpos_s:" << zpos_s<<" zpos_p:"<<zpos_p<<" zpos_n:"<<zpos_n<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                        x_plus=0.5*p->DXP[IP];
                    }
                    
                    x_defin=coeff_l*x_quad+coeff_r*x_plus;
                    
                    teta=(0.0-x_defin)/p->DXP[IM1];
                    Fi_quad=a_fi*x_defin*x_defin+b_fi*x_defin+c_fi;
                    Fi_plus=a_fi_plus*x_defin*x_defin+b_fi_plus*x_defin+c_fi_plus;
                    Fival=coeff_l*Fi_quad+coeff_r*Fi_plus;
                    
                    if(p->flag4[Im1JK]==AIR)
                    {
                       // if(p->flag4[IJKp1]==AIR)
                         //   cout<<"Air n LRT x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
                       // else
                       //     cout<<"AIR n LR!!! x:"<<p->XP[IP]<<" z-FS:"<<p->ZP[KP]-p->F60<<" eta_s:" << eta(i-1,j)<<" eta_p:"<<eta(i,j)<<" eta_n:"<<eta(i+1,j)<<" DX_s:"<<p->DXP[IM1]<<" DX_n:"<<p->DXP[IP]<<endl;
            
                       Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                        
                        teta = teta>1.0e-6?teta:1.0e20;
            
                        a->M.p[n] -= 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.p[n] += 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->M.n[n] += 1.0/(p->DXP[IP]*p->DXN[IP]);
                        a->M.n[n] -= 1.0/(teta*p->DXP[IP]*p->DXN[IP]);
                
                        a->rhsvec.V[n] -= a->M.n[n]*Fival;
                        a->M.n[n] = 0.0;
                        
                        
                    }
                    
                    else
                    {
                        
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                    Fival=a_fi*x_quad*x_quad+b_fi*x_quad+c_fi;
                    Z_t=p->DXP[IP];
                    Z_b=p->DXP[IM1];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.s[n]+=(M_b_num/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DXP[IP]*p->DXN[IP]);
                    a->M.n[n]=0.0;
                    }
                }
            }

            // east
            if(p->flag4[IJm1K]==AIR)
            {
                if(p->A323<=21)
                {
                    a->rhsvec.V[n] -= a->M.e[n]*f(i,j-1,k);
                    a->M.e[n] = 0.0;
                }
                
                if(p->A323==77)
                {
                    double zpos_p,zpos_w,zpos_e;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,y1_quad,y2_quad,y_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
            
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_e=eta(i,j-1)-(p->ZP[KP]-p->F60);
                    zpos_w=eta(i,j+1)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    b_quad_num=zpos_w*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-zpos_e-zpos_p*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-b_quad/p->DYP[JP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    y1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    y2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    b_fi_num=Fifsf(i,j+1)*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j-1)-Fifsf(i,j)*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-b_fi/p->DYP[JP];
                    c_fi=Fifsf(i,j);
                    
                    if(y1_quad<=0.0 && y1_quad>=0.0-(p->DYP[JM1]) && !(y2_quad<=0.0 && y2_quad>=0.0-p->DYP[JM1]))
                        y_quad=y1_quad;
                    else if(y2_quad<=0.0 && y2_quad>=0.0-(p->DYP[JM1]) && !(y1_quad<=0.0 && y1_quad>=0.0-p->DYP[JM1]))
                        y_quad=y2_quad;
                    else if(y1_quad <= y2_quad+1.0e-06 && y1_quad >= y2_quad-1.0e-06)
                    {
                        cout<<"equal e"<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        y_quad=y1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo e: "<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DX_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        y_quad=-p->DYP[JM1];
                    }
                    teta=(0.0-y_quad)/p->DYP[JM1];
                    teta = teta>1.0e-6?teta:1.0e20;
                    Fival=a_fi*y_quad*y_quad+b_fi*y_quad+c_fi;
                
                    a->M.p[n] -= 1.0/(p->DYP[JM1]*p->DYN[JP]);
                    a->M.p[n] += 1.0/(teta*p->DYP[JM1]*p->DYN[JP]);
                           
                    a->M.e[n] += 1.0/(p->DYP[JM1]*p->DYN[JP]);
                    a->M.e[n] -= 1.0/(teta*p->DYP[JM1]*p->DYN[JP]);
                
                    a->rhsvec.V[n] -= a->M.e[n]*Fival;
                    a->M.e[n] = 0.0;
                }
                
                if(p->A323==78)
                {
                   double zpos_p,zpos_w,zpos_e,zpos_ee;
                    double Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,y1_quad,y2_quad,y_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    
                    if(p->flag4[IJp1K]<OUTFLOW)
                    {
                        zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                        zpos_e=eta(i,j-1)-(p->ZP[KP]-p->F60);
                        zpos_ee=eta(i,j-2)-(p->ZP[KP]-p->F60);
                        
                        b_quad_denom=(-1.0*p->DYP[JM1])+(-1.0*p->DYP[JM2])-((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(-1.0*p->DYP[JM1]);
                        b_quad_num=zpos_ee-zpos_p+zpos_p*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1])-zpos_e*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1]);
                        b_quad=b_quad_num/b_quad_denom;
                        c_quad=zpos_p;
                        a_quad=zpos_e/(p->DYP[JM1]*p->DYP[JM1])-zpos_p/(p->DYP[JM1]*p->DYP[JM1])-b_quad/(-1.0*p->DYP[JM1]);
                        p_quad=b_quad/a_quad;
                        q_quad=c_quad/a_quad;
                        y1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                        y2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                        
                        b_fi_denom=(-1.0*p->DYP[JM1])+(-1.0*p->DYP[JM2])-((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(-1.0*p->DYP[JM1]);
                        b_fi_num=Fifsf(i,j-2)-Fifsf(i,j)+Fifsf(i,j)*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1])-Fifsf(i,j-1)*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1]);
                        b_fi=b_fi_num/b_fi_denom;
                        c_fi=Fifsf(i,j);
                        a_fi=Fifsf(i,j-1)/(p->DYP[JM1]*p->DYP[JM1])-Fifsf(i,j)/(p->DYP[JM1]*p->DYP[JM1])-b_fi/(-1.0*p->DYP[JM1]);
                    }
                    
                    else
                    {
                        zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                        zpos_e=eta(i,j-1)-(p->ZP[KP]-p->F60);
                        zpos_w=eta(i,j+1)-(p->ZP[KP]-p->F60);
                    
                        b_quad_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                        b_quad_num=zpos_w*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-zpos_e-zpos_p*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                        b_quad=b_quad_num/b_quad_denom;
                        a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-b_quad/p->DYP[JP];
                        c_quad=zpos_p;
                        p_quad=b_quad/a_quad;
                        q_quad=c_quad/a_quad;
                        y1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                        y2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                        
                        b_fi_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                        b_fi_num=Fifsf(i,j+1)*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j-1)-Fifsf(i,j)*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                        b_fi=b_fi_num/b_fi_denom;
                        a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-b_fi/p->DYP[JP];
                        c_fi=Fifsf(i,j);
                    }
                    
                    if(y1_quad<=0.0 && y1_quad>=0.0-(p->DYP[JM1]) && !(y2_quad<=0.0 && y2_quad>=0.0-p->DYP[JM1]))
                        y_quad=y1_quad;
                    else if(y2_quad<=0.0 && y2_quad>=0.0-(p->DYP[JM1]) && !(y1_quad<=0.0 && y1_quad>=0.0-p->DYP[JM1]))
                        y_quad=y2_quad;
                    else if(y1_quad <= y2_quad+1.0e-06 && y1_quad >= y2_quad-1.0e-06)
                    {
                    //    cout<<"equal e"<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        y_quad=y1_quad;
                    }
                    else
                    {
                    //    cout<<"mucho problemo e: "<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DX_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        y_quad=-p->DYP[JM1];
                    }
                    teta=(0.0-y_quad)/p->DYP[JM1];
                    teta = teta>1.0e-6?teta:1.0e20;
                    Fival=a_fi*y_quad*y_quad+b_fi*y_quad+c_fi;
                
                    a->M.p[n] -= 1.0/(p->DYP[JM1]*p->DYN[JP]);
                    a->M.p[n] += 1.0/(teta*p->DYP[JM1]*p->DYN[JP]);
                           
                    a->M.e[n] += 1.0/(p->DYP[JM1]*p->DYN[JP]);
                    a->M.e[n] -= 1.0/(teta*p->DYP[JM1]*p->DYN[JP]);
                
                    a->rhsvec.V[n] -= a->M.e[n]*Fival;
                    a->M.e[n] = 0.0;
                }
                
                if(p->A323==79)
                {
                    double zpos_p,zpos_w,zpos_e,zpos_ee;
                    double e_Fival,p_Fival,xpos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double e_a_quad,e_b_quad,e_b_quad_num,e_b_quad_denom,e_c_quad,e_y1_quad,e_y2_quad,e_y_quad,e_p_quad,e_q_quad;
                    double e_a_fi,e_b_fi_b,e_b_fi_denom,e_b_fi_num,e_c_fi,e_b_fi;
                    double p_a_quad,p_b_quad,p_b_quad_num,p_b_quad_denom,p_c_quad,p_y1_quad,p_y2_quad,p_y_quad,p_p_quad,p_q_quad;
                    double p_a_fi,p_b_fi_denom,p_b_fi_num,p_c_fi,p_b_fi;
                    double Fival,y_quad;
                    
                    zpos_w=eta(i,j+1)-(p->ZP[KP]-p->F60);
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_e=eta(i,j-1)-(p->ZP[KP]-p->F60);
                    zpos_ee=eta(i,j-2)-(p->ZP[KP]-p->F60);
                        
                    e_b_quad_denom=(-1.0*p->DYP[JM1])+(-1.0*p->DYP[JM2])-((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(-1.0*p->DYP[JM1]);
                    e_b_quad_num=zpos_ee-zpos_p+zpos_p*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1])-zpos_e*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1]);
                    e_b_quad=e_b_quad_num/e_b_quad_denom;
                    e_c_quad=zpos_p;
                    e_a_quad=zpos_e/(p->DYP[JM1]*p->DYP[JM1])-zpos_p/(p->DYP[JM1]*p->DYP[JM1])-e_b_quad/(-1.0*p->DYP[JM1]);
                    e_p_quad=e_b_quad/e_a_quad;
                    e_q_quad=e_c_quad/e_a_quad;
                    e_y1_quad=(0.0-e_p_quad)/2.0+sqrt(((e_p_quad/2.0)*(e_p_quad/2.0)-e_q_quad));
                    e_y2_quad=(0.0-e_p_quad)/2.0-sqrt(((e_p_quad/2.0)*(e_p_quad/2.0)-e_q_quad));
                        
                    e_b_fi_denom=(-1.0*p->DYP[JM1])+(-1.0*p->DYP[JM2])-((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(-1.0*p->DYP[JM1]);
                    e_b_fi_num=Fifsf(i,j-2)-Fifsf(i,j)+Fifsf(i,j)*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1])-Fifsf(i,j-1)*((p->DYP[JM1]+p->DYP[JM2])*(p->DYP[JM1]+p->DYP[JM2]))/(p->DYP[JM1]*p->DYP[JM1]);
                    e_b_fi=e_b_fi_num/e_b_fi_denom;
                    e_c_fi=Fifsf(i,j);
                    e_a_fi=Fifsf(i,j-1)/(p->DYP[JM1]*p->DYP[JM1])-Fifsf(i,j)/(p->DYP[JM1]*p->DYP[JM1])-e_b_fi/(-1.0*p->DYP[JM1]);
                    
                    p_b_quad_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    p_b_quad_num=zpos_w*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-zpos_e-zpos_p*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    p_b_quad=p_b_quad_num/p_b_quad_denom;
                    p_a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-p_b_quad/p->DYP[JP];
                    p_c_quad=zpos_p;
                    p_p_quad=p_b_quad/p_a_quad;
                    p_q_quad=p_c_quad/p_a_quad;
                    p_y1_quad=(0.0-p_p_quad)/2.0+sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    p_y2_quad=(0.0-p_p_quad)/2.0-sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                        
                    p_b_fi_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    p_b_fi_num=Fifsf(i,j+1)*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j-1)-Fifsf(i,j)*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    p_b_fi=p_b_fi_num/p_b_fi_denom;
                    p_a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-p_b_fi/p->DYP[JP];
                    p_c_fi=Fifsf(i,j);
                    
                    if(e_y1_quad<=0.0 && e_y1_quad>=0.0-(p->DYP[JM1]) && !(e_y2_quad<=0.0 && e_y2_quad>=0.0-p->DYP[JM1]))
                        e_y_quad=e_y1_quad;
                    else if(e_y2_quad<=0.0 && e_y2_quad>=0.0-(p->DYP[JM1]) && !(e_y1_quad<=0.0 && e_y1_quad>=0.0-p->DYP[JM1]))
                        e_y_quad=e_y2_quad;
                    else if(e_y1_quad <= e_y2_quad+1.0e-06 && e_y1_quad >= e_y2_quad-1.0e-06)
                    {
                        cout<<"equal e e"<< "e y1:"<<e_y1_quad<<"e y2:"<<e_y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_ee:"<<zpos_ee<<" DY_e:"<<p->DYP[JM1]<<" DY_ee:"<<p->DYP[JM2]<<endl;
                        e_y_quad=e_y1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo e e: "<< "e y1:"<<e_y1_quad<<"e y2:"<<e_y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_ee:"<<zpos_ee<<" DX_e:"<<p->DYP[JM1]<<" DY_ee:"<<p->DYP[JM2]<<endl;
                        e_y_quad=-p->DYP[JM1];
                    }
                    
                    if(p_y1_quad<=0.0 && p_y1_quad>=0.0-(p->DYP[JM1]) && !(p_y2_quad<=0.0 && p_y2_quad>=0.0-p->DYP[JM1]))
                       p_y_quad=p_y1_quad;
                    else if(p_y2_quad<=0.0 && p_y2_quad>=0.0-(p->DYP[JM1]) && !(p_y1_quad<=0.0 && p_y1_quad>=0.0-p->DYP[JM1]))
                        p_y_quad=p_y2_quad;
                    else if(p_y1_quad <= p_y2_quad+1.0e-06 && p_y1_quad >= p_y2_quad-1.0e-06)
                    {
                        cout<<"equal e p"<< "p y1:"<<p_y1_quad<<"p y2:"<<p_y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        p_y_quad=p_y1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo e p: "<< "p y1:"<<p_y1_quad<<"p y2:"<<p_y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DX_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        p_y_quad=-p->DYP[JM1];
                    }
                    
                    y_quad=(e_y_quad+p_y_quad)/2;
                    teta=(0.0-y_quad)/p->DYP[JM1];
                    teta = teta>1.0e-6?teta:1.0e20;
                    
                    e_Fival=e_a_fi*y_quad*y_quad+e_b_fi*y_quad+e_c_fi;
                    p_Fival=p_a_fi*y_quad*y_quad+p_b_fi*y_quad+p_c_fi;
                    Fival=(e_Fival+p_Fival)/2;
                
                    a->M.p[n] -= 1.0/(p->DYP[JM1]*p->DYN[JP]);
                    a->M.p[n] += 1.0/(teta*p->DYP[JM1]*p->DYN[JP]);
                           
                    a->M.e[n] += 1.0/(p->DYP[JM1]*p->DYN[JP]);
                    a->M.e[n] -= 1.0/(teta*p->DYP[JM1]*p->DYN[JP]);
                
                    a->rhsvec.V[n] -= a->M.e[n]*Fival;
                    a->M.e[n] = 0.0;
                }
            }

            // west
            if(p->flag4[IJp1K]==AIR)
            {
                if(p->A323<=21)
                {
                    a->rhsvec.V[n] -= a->M.w[n]*f(i,j+1,k);
                    a->M.w[n] = 0.0;
                }
                
                if(p->A323==77)
                {
                    double zpos_p,zpos_w,zpos_e;
                    double Fival,ypos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,y1_quad,y2_quad,y_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
            
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_e=eta(i,j-1)-(p->ZP[KP]-p->F60);
                    zpos_w=eta(i,j+1)-(p->ZP[KP]-p->F60);
                    
                    b_quad_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    b_quad_num=zpos_w*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-zpos_e-zpos_p*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    b_quad=b_quad_num/b_quad_denom;
                    a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-b_quad/p->DYP[JP];
                    c_quad=zpos_p;
                    p_quad=b_quad/a_quad;
                    q_quad=c_quad/a_quad;
                    y1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    y2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                    b_fi_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    b_fi_num=Fifsf(i,j+1)*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j-1)-Fifsf(i,j)*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    b_fi=b_fi_num/b_fi_denom;
                    a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-b_fi/p->DYP[JP];
                    c_fi=Fifsf(i,j);
                    
                    if(y1_quad>=0.0 && y1_quad<=p->DYP[JP] && !(y2_quad>=0.0 && y2_quad<=p->DYP[JP]))
                        y_quad=y1_quad;
                    else if(y2_quad>=0.0 && y2_quad<=p->DYP[JP] && !(y1_quad>=0.0 && y1_quad<=p->DYP[JP]))
                        y_quad=y2_quad;
                    else if(y1_quad <= y2_quad+1.0e-06 && y1_quad >= y2_quad-1.0e-06)
                    {
                    //    cout<<"equal w"<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;;
                        y_quad=y1_quad;
                    }
                    else
                    {
                       // cout<<"mucho problemo w: "<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        y_quad=p->DYP[JP];
                    }
                    
                    teta=(y_quad)/p->DYP[JP];
                    Fival=a_fi*y_quad*y_quad+b_fi*y_quad+c_fi;
                    
                    teta = teta>1.0e-6?teta:1.0e20;
            
                    a->M.p[n] -= 1.0/(p->DYP[JP]*p->DYN[JP]);
                    a->M.p[n] += 1.0/(teta*p->DYP[JP]*p->DYN[JP]);
                
                    a->M.w[n] += 1.0/(p->DYP[JP]*p->DYN[JP]);
                    a->M.w[n] -= 1.0/(teta*p->DYP[JP]*p->DYN[JP]);
                
                    a->rhsvec.V[n] -= a->M.w[n]*Fival;
                    a->M.w[n] = 0.0;
                }
                
                if(p->A323==78)
                {
                    double zpos_p,zpos_w,zpos_ww,zpos_e;
                    double Fival,ypos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double a_quad,b_quad,b_quad_num,b_quad_denom,c_quad,y1_quad,y2_quad,y_quad,p_quad,q_quad;
                    double a_fi,b_fi_b,b_fi_denom,b_fi_num,c_fi,b_fi;
                    
                    if(p->flag4[IJm1K]<OUTFLOW)
                    {
                        zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                        zpos_w=eta(i,j+1)-(p->ZP[KP]-p->F60);
                        zpos_ww=eta(i,j+2)-(p->ZP[KP]-p->F60);
                        
                        b_quad_denom=p->DYP[JP]+p->DYP[JP1]-((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/p->DYP[JP];
                        b_quad_num=zpos_ww-zpos_p+zpos_p*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP])-zpos_w*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP]);
                        b_quad=b_quad_num/b_quad_denom;
                        c_quad=zpos_p;
                        a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-b_quad/p->DYP[JP];
                        p_quad=b_quad/a_quad;
                        q_quad=c_quad/a_quad;
                        y1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                        y2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                        
                        b_fi_denom=p->DYP[JP]+p->DYP[JP1]-((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/p->DYP[JP];
                        b_fi_num=Fifsf(i,j+2)-Fifsf(i,j)+Fifsf(i,j)*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j+1)*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP]);
                        b_fi=b_fi_num/b_fi_denom;
                        c_fi=Fifsf(i,j);
                        a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-b_fi/p->DYP[JP];
                    }
                    
                    else
                    {
                        zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                        zpos_e=eta(i,j-1)-(p->ZP[KP]-p->F60);
                        zpos_w=eta(i,j+1)-(p->ZP[KP]-p->F60);
                    
                        b_quad_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                        b_quad_num=zpos_w*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-zpos_e-zpos_p*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                        b_quad=b_quad_num/b_quad_denom;
                        a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-b_quad/p->DYP[JP];
                        c_quad=zpos_p;
                        p_quad=b_quad/a_quad;
                        q_quad=c_quad/a_quad;
                        y1_quad=(0.0-p_quad)/2.0+sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                        y2_quad=(0.0-p_quad)/2.0-sqrt(((p_quad/2.0)*(p_quad/2.0)-q_quad));
                    
                        b_fi_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                        b_fi_num=Fifsf(i,j+1)*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j-1)-Fifsf(i,j)*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                        b_fi=b_fi_num/b_fi_denom;
                        a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-b_fi/p->DYP[JP];
                        c_fi=Fifsf(i,j);
                    }
                    
                    if(y1_quad>=0.0 && y1_quad<=p->DYP[JP] && !(y2_quad>=0.0 && y2_quad<=p->DYP[JP]))
                        y_quad=y1_quad;
                    else if(y2_quad>=0.0 && y2_quad<=p->DYP[JP] && !(y1_quad>=0.0 && y1_quad<=p->DYP[JP]))
                        y_quad=y2_quad;
                    else if(y1_quad <= y2_quad+1.0e-06 && y1_quad >= y2_quad-1.0e-06)
                    {
                   //     cout<<"equal w"<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;;
                        y_quad=y1_quad;
                    }
                    else
                    {
                    //    cout<<"mucho problemo w: "<< " y1:"<<y1_quad<<" y2:"<<y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        y_quad=p->DYP[JP];
                    }
                    
                    teta=(y_quad)/p->DYP[JP];
                    Fival=a_fi*y_quad*y_quad+b_fi*y_quad+c_fi;
                    
                    teta = teta>1.0e-6?teta:1.0e20;
            
                    a->M.p[n] -= 1.0/(p->DYP[JP]*p->DYN[JP]);
                    a->M.p[n] += 1.0/(teta*p->DYP[JP]*p->DYN[JP]);
                
                    a->M.w[n] += 1.0/(p->DYP[JP]*p->DYN[JP]);
                    a->M.w[n] -= 1.0/(teta*p->DYP[JP]*p->DYN[JP]);
                
                    a->rhsvec.V[n] -= a->M.w[n]*Fival;
                    a->M.w[n] = 0.0;
                }
                
                if(p->A323==79)
                {
                    double zpos_p,zpos_w,zpos_ww,zpos_e;
                    double w_Fival,p_Fival,ypos_zero,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    double w_a_quad,w_b_quad,w_b_quad_num,w_b_quad_denom,w_c_quad,w_y1_quad,w_y2_quad,w_y_quad,w_p_quad,w_q_quad;
                    double w_a_fi,w_b_fi_denom,w_b_fi_num,w_c_fi,w_b_fi;
                    double p_a_quad,p_b_quad,p_b_quad_num,p_b_quad_denom,p_c_quad,p_y1_quad,p_y2_quad,p_y_quad,p_p_quad,p_q_quad;
                    double p_a_fi,p_b_fi_denom,p_b_fi_num,p_c_fi,p_b_fi;
                    double Fival, y_quad;
                    
                    zpos_e=eta(i,j-1)-(p->ZP[KP]-p->F60);
                    zpos_p=eta(i,j)-(p->ZP[KP]-p->F60);
                    zpos_w=eta(i,j+1)-(p->ZP[KP]-p->F60);
                    zpos_ww=eta(i,j+2)-(p->ZP[KP]-p->F60);
                        
                    w_b_quad_denom=p->DYP[JP]+p->DYP[JP1]-((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/p->DYP[JP];
                    w_b_quad_num=zpos_ww-zpos_p+zpos_p*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP])-zpos_w*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP]);
                    w_b_quad=w_b_quad_num/w_b_quad_denom;
                    w_c_quad=zpos_p;
                    w_a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-w_b_quad/p->DYP[JP];
                    w_p_quad=w_b_quad/w_a_quad;
                    w_q_quad=w_c_quad/w_a_quad;
                    w_y1_quad=(0.0-w_p_quad)/2.0+sqrt(((w_p_quad/2.0)*(w_p_quad/2.0)-w_q_quad));
                    w_y2_quad=(0.0-w_p_quad)/2.0-sqrt(((w_p_quad/2.0)*(w_p_quad/2.0)-w_q_quad));
                        
                    w_b_fi_denom=p->DYP[JP]+p->DYP[JP1]-((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/p->DYP[JP];
                    w_b_fi_num=Fifsf(i,j+2)-Fifsf(i,j)+Fifsf(i,j)*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j+1)*((p->DYP[JP]+p->DYP[JP1])*(p->DYP[JP]+p->DYP[JP1]))/(p->DYP[JP]*p->DYP[JP]);
                    w_b_fi=w_b_fi_num/w_b_fi_denom;
                    w_c_fi=Fifsf(i,j);
                    w_a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-w_b_fi/p->DYP[JP];
                    
                    p_b_quad_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    p_b_quad_num=zpos_w*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-zpos_e-zpos_p*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    p_b_quad=p_b_quad_num/p_b_quad_denom;
                    p_a_quad=zpos_w/(p->DYP[JP]*p->DYP[JP])-zpos_p/(p->DYP[JP]*p->DYP[JP])-p_b_quad/p->DYP[JP];
                    p_c_quad=zpos_p;
                    p_p_quad=p_b_quad/p_a_quad;
                    p_q_quad=p_c_quad/p_a_quad;
                    p_y1_quad=(0.0-p_p_quad)/2.0+sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    p_y2_quad=(0.0-p_p_quad)/2.0-sqrt(((p_p_quad/2.0)*(p_p_quad/2.0)-p_q_quad));
                    
                    p_b_fi_denom=p->DYP[JM1]*p->DYP[JM1]/p->DYP[JP]+p->DYP[JM1];
                    p_b_fi_num=Fifsf(i,j+1)*p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j-1)-Fifsf(i,j)*(p->DYP[JM1]*p->DYP[JM1]/(p->DYP[JP]*p->DYP[JP])-1.0);
                    p_b_fi=p_b_fi_num/p_b_fi_denom;
                    p_a_fi=Fifsf(i,j+1)/(p->DYP[JP]*p->DYP[JP])-Fifsf(i,j)/(p->DYP[JP]*p->DYP[JP])-p_b_fi/p->DYP[JP];
                    p_c_fi=Fifsf(i,j);
                    
                    
                    if(w_y1_quad>=0.0 && w_y1_quad<=p->DYP[JP] && !(w_y2_quad>=0.0 && w_y2_quad<=p->DYP[JP]))
                        w_y_quad=w_y1_quad;
                    else if(w_y2_quad>=0.0 && w_y2_quad<=p->DYP[JP] && !(w_y1_quad>=0.0 && w_y1_quad<=p->DYP[JP]))
                        w_y_quad=w_y2_quad;
                    else if(w_y1_quad <= w_y2_quad+1.0e-06 && w_y1_quad >= w_y2_quad-1.0e-06)
                    {
                        cout<<"equal w w"<< "w y1:"<<w_y1_quad<<"w y2:"<<w_y2_quad<<" zpos_ww:" << zpos_ww<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_ww:"<<p->DYP[JP1]<<" DY_w:"<<p->DYP[JP]<<endl;;
                        w_y_quad=w_y1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo w w: "<< "w y1:"<<w_y1_quad<<"w y2:"<<w_y2_quad<<" zpos_ww:" << zpos_ww<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_ww:"<<p->DYP[JP1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        w_y_quad=p->DYP[JP];
                    }
                    
                    if(p_y1_quad>=0.0 && p_y1_quad<=p->DYP[JP] && !(p_y2_quad>=0.0 && p_y2_quad<=p->DYP[JP]))
                        p_y_quad=p_y1_quad;
                    else if(p_y2_quad>=0.0 && p_y2_quad<=p->DYP[JP] && !(p_y1_quad>=0.0 && p_y1_quad<=p->DYP[JP]))
                        p_y_quad=p_y2_quad;
                    else if(p_y1_quad <= p_y2_quad+1.0e-06 && p_y1_quad >= p_y2_quad-1.0e-06)
                    {
                        cout<<"equal w p"<< "p y1:"<<p_y1_quad<<"p y2:"<<p_y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;;
                        p_y_quad=p_y1_quad;
                    }
                    else
                    {
                        cout<<"mucho problemo w p: "<< "p y1:"<<p_y1_quad<<"p y2:"<<p_y2_quad<<" zpos_e:" << zpos_e<<" zpos_p:"<<zpos_p<<" zpos_w:"<<zpos_w<<" DY_e:"<<p->DYP[JM1]<<" DY_w:"<<p->DYP[JP]<<endl;
                        p_y_quad=p->DYP[JP];
                    }
                    
                    y_quad=(w_y_quad+p_y_quad)/2;
                    teta=(y_quad)/p->DYP[JP];
                    
                    w_Fival=w_a_fi*y_quad*y_quad+w_b_fi*y_quad+w_c_fi;
                    p_Fival=p_a_fi*y_quad*y_quad+p_b_fi*y_quad+p_c_fi;
                    Fival=(w_Fival+p_Fival)/2;
                    
                    teta = teta>1.0e-6?teta:1.0e20;
            
                    a->M.p[n] -= 1.0/(p->DYP[JP]*p->DYN[JP]);
                    a->M.p[n] += 1.0/(teta*p->DYP[JP]*p->DYN[JP]);
                
                    a->M.w[n] += 1.0/(p->DYP[JP]*p->DYN[JP]);
                    a->M.w[n] -= 1.0/(teta*p->DYP[JP]*p->DYN[JP]);
                
                    a->rhsvec.V[n] -= a->M.w[n]*Fival;
                    a->M.w[n] = 0.0;
                }
            }

            // Free Surface BC
            if(p->flag4[IJKp1]==AIR)
            {
                // -----------
                if(p->A323==1)
                {
                a->rhsvec.V[n] -= a->M.t[n]*Fifsf(i,j);
                a->M.t[n] = 0.0;
                }
                
                // -----------
                if(p->A323==2)
                {
                double lsv0,lsv1;

                lsv0 = fabs(a->phi(i,j,k));
                lsv1 = fabs(a->phi(i,j,k+1));

                lsv0 = fabs(lsv0)>1.0e-6?lsv0:1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));

                a->rhsvec.V[n] -= a->M.t[n]*Fifsf(i,j)*(1.0 + lsv1/lsv0);
                a->M.p[n] -= a->M.t[n]*lsv1/lsv0;
                a->M.t[n] = 0.0;
                }
                
                
                // -----------
                if(p->A323==3)
                {
                double x0,x1,x2,y2;
                double x,y;
                double Lx0,Lx1,Lx2;
                double denom1,denom2,denom3,denom4,denom5,denom6;

                x0 = -fabs(a->phi(i,j,k-1));
                x1 = -fabs(a->phi(i,j,k));
                x2 = 0.0;

                y2 = Fifsf(i,j);
                
                denom1 = fabs(x0-x1)>1.0e-6?(x0-x1):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
                denom2 = fabs(x1-x0)>1.0e-6?(x1-x0):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
                denom3 = fabs(x2-x0)>1.0e-6?(x2-x0):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
                
                denom4 = fabs(x0-x2)>1.0e-6?(x0-x2):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
                denom5 = fabs(x1-x2)>1.0e-6?(x1-x2):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
                denom6 = fabs(x2-x1)>1.0e-6?(x2-x1):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
    

                x = fabs(a->phi(i,j,k+1));

                Lx0 = ((x-x1)/denom1) * ((x-x2)/denom4);
                Lx1 = ((x-x0)/denom2) * ((x-x2)/denom5);
                Lx2 = ((x-x0)/denom3) * ((x-x1)/denom6);

                a->rhsvec.V[n]  -= a->M.t[n]*Lx2*y2;
                a->M.p[n]       += a->M.t[n]*Lx1;
                a->M.b[n]       += a->M.t[n]*Lx0;
                a->M.t[n]       = 0.0;
                }
                
                // -----------
                if(p->A323==4)
                {
                    
                double teta;
                teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
                
                teta = teta>1.0e-6?teta:1.0e20;
                
                //cout<<" Teta: "<<teta<<" a->phi(i,j,k): "<<a->phi(i,j,k)<<" a->phi(i+1,j,k): "<<a->phi(i+1,j,k)<<endl;
                
                a->M.p[n] -= 1.0/(p->DZP[KP]*p->DZN[KP]);
                a->M.p[n] += 1.0/(teta*p->DZP[KP]*p->DZN[KP]);
                           
                a->M.t[n] += 1.0/(p->DZP[KP]*p->DZN[KP]);
                a->M.t[n] -= 1.0/(teta*p->DZP[KP]*p->DZN[KP]);
                
                a->rhsvec.V[n] -= a->M.t[n]*Fifsf(i,j);
                a->M.t[n] = 0.0;
                }
                
                
                if(p->A323==5 || p->A323==6)
                {
                    
                double teta;
                teta = (eta(i,j)-(p->ZP[KP]-p->F60))/p->DZP[KP];
                teta = teta>1.0e-6?teta:1.0e20;
                
                a->M.p[n] -= 1.0/(p->DZP[KP]*p->DZN[KP]);
                a->M.p[n] += 1.0/(teta*p->DZP[KP]*p->DZN[KP]);
                           
                a->M.t[n] += 1.0/(p->DZP[KP]*p->DZN[KP]);
                a->M.t[n] -= 1.0/(teta*p->DZP[KP]*p->DZN[KP]);
                
                a->rhsvec.V[n] -= a->M.t[n]*Fifsf(i,j);
                a->M.t[n] = 0.0;
                }
                
                if(p->A323==7 || p->A323==77 || p->A323==8 || p->A323==9 || p->A323==10 || p->A323==11 || p->A323==12 || p->A323==78 || p->A323==79)
                {
                    double Fival,teta;
                    double Z_t,Z_b,M_p_num,M_b_num,denom;
                    
                    teta = (eta(i,j)-(p->ZP[KP]-p->F60))/p->DZP[KP];
                    
                    if(teta<1.0e-7)
                        teta=1.0e-7;
                    Fival = Fifsf(i,j);
                    
                    Z_t=p->DZP[KP];
                    Z_b=p->DZP[KM1];
                    
                    denom=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)-(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b);
                        
                    M_p_num=(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b)+(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t);
                        
                    M_b_num=1.0-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_t+Z_t*Z_t)+(Z_b+teta*Z_t)*(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_b*Z_t+Z_b*Z_b*Z_t*Z_t)-(Z_b+teta*Z_t)*(Z_b+Z_t)*(Z_b+Z_t)/(Z_b*Z_b*Z_t+Z_b*Z_t*Z_t)+(Z_b+teta*Z_t)/(Z_t+Z_t*Z_t/Z_b)-(Z_b+teta*Z_t)*(Z_b+teta*Z_t)/(Z_b*Z_b);
                        
                    a->M.p[n]+=(M_p_num/denom)/(p->DZP[KP]*p->DZN[KP]);
                    a->M.b[n]+=(M_b_num/denom)/(p->DZP[KP]*p->DZN[KP]);
                    a->rhsvec.V[n]+=(Fival/denom)/(p->DZP[KP]*p->DZN[KP]);
                    a->M.t[n]=0.0;
                }
                
  
            }

            // KBEDBC
            if(p->flag4[IJKm1]<AIR)
            {
            a->M.p[n] += a->M.b[n];
            a->M.b[n] = 0.0;
            }
        }
	++n;
	}
    
    double starttime=pgc->timer();
    psolv->start(p,a,pgc,a->Fi,a->rhsvec,5);
    double endtime=pgc->timer();
    pgc->start4(p,a->Fi,250);

    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"Fi_iter: "<<p->solveriter<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
}