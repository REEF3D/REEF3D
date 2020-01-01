/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"density_vof.h"
#include"lexer.h"
#include"fdm.h"

density_vof::density_vof(lexer* p) : epsi(p->F45*p->dx), eps(2.1*p->dx)
{
}

density_vof::~density_vof()
{
}

double density_vof::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi,r,s;
	
	ii = aa-aa/(fabs(aa)>0?fabs(aa):1);
	jj = bb-bb/(fabs(bb)>0?fabs(bb):1);
	kk = cc-cc/(fabs(cc)>0?fabs(cc):1);
	
	
	if(p->D32==1)
	{
       
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

        psi = p->F45*p->DXM;
      
        if(phival>psi)
        H=1.0;

        if(phival<-psi)
        H=0.0;

        if(fabs(phival)<=psi)
        H=0.5*(1.0 + phival/psi + (1.0/PI)*sin((PI*phival)/psi));
        
            
        roval = p->W1*H + p->W3*(1.0-H);
	}
	
	

	if(p->D32==2)
	{
       
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
        
        if(p->j_dir==0)
        psi = p->F45*(1.0/2.0)*(p->DXN[IP]+p->DZN[KP]);
        
        if(p->j_dir==1)
        psi = p->F45*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
    
        if(phival>psi)
        H=1.0;

        if(phival<-psi)
        H=0.0;

        if(fabs(phival)<=psi)
        H=0.5*(1.0 + phival/psi + (1.0/PI)*sin((PI*phival)/psi));
        
            
        roval = p->W1*H + p->W3*(1.0-H);
	}
	
	// -----
	
	if(p->D32==3)
	roval = 0.5*(a->ro(i+ii,j+jj,k+kk) + a->ro(i+aa,j+bb,k+cc));
	
	// -----
	
	if(p->D32==4)
	{
        
        double fx,fy,fz,fn;

        
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
        
        if(aa<0)
        {
        fx = fabs((a->phi(i,j,k) - a->phi(i-1,j,k))/(p->DXP[IM1]));
        fy = fabs( 0.5*(a->phi(i,j+1,k) + a->phi(i-1,j+1,k)) - 0.5*(a->phi(i,j-1,k)+a->phi(i-1,j-1,k))/(p->DYP[JM1]+p->DYP[JP]));
        fz = fabs( 0.5*(a->phi(i,j,k+1) + a->phi(i-1,j,k+1)) - 0.5*(a->phi(i,j,k-1)+a->phi(i-1,j,k-1))/(p->DZP[KM1]+p->DZP[KP]));
        
        fn = sqrt(fx*fx + fy*fy + fz*fz);
        
        fx = fx/fn;
        fy = fy/fn;
        fz = fz/fn;
        
        psi = p->F45*sqrt(pow(fx*p->DXP[IP],2.0) + pow(fy*p->DYN[JP],2.0) + pow(fz*p->DZN[KP],2.0));
        }
        
        if(aa>0)
        {
        fx = fabs((a->phi(i+1,j,k) - a->phi(i,j,k))/(p->DXP[IP]));
        fy = fabs( 0.5*(a->phi(i,j+1,k) + a->phi(i+1,j+1,k)) - 0.5*(a->phi(i,j-1,k)+a->phi(i+1,j-1,k))/(p->DYP[JM1]+p->DYP[JP]));
        fz = fabs( 0.5*(a->phi(i,j,k+1) + a->phi(i+1,j,k+1)) - 0.5*(a->phi(i,j,k-1)+a->phi(i+1,j,k-1))/(p->DZP[KM1]+p->DZP[KP]));
        
        fn = sqrt(fx*fx + fy*fy + fz*fz);
        
        fx = fx/fn;
        fy = fy/fn;
        fz = fz/fn;
        
        psi = p->F45*sqrt(pow(fx*p->DXP[IP],2.0) + pow(fy*p->DYN[JP],2.0) + pow(fz*p->DZN[KP],2.0));
        }
        
        
        if(bb<0)
        {
        fx = fabs( 0.5*(a->phi(i+1,j,k) + a->phi(i+1,j-1,k)) - 0.5*(a->phi(i-1,j,k)+a->phi(i-1,j-1,k))/(p->DXP[IM1]+p->DXP[IP]));
        fy = fabs((a->phi(i,j,k) - a->phi(i,j-1,k))/p->DYP[JM1]);      
        fz = fabs( 0.5*(a->phi(i,j,k+1) + a->phi(i,j-1,k+1)) - 0.5*(a->phi(i,j,k-1)+a->phi(i,j-1,k-1))/(p->DZP[KM1]+p->DZP[KP]));
        
        fn = sqrt(fx*fx + fy*fy + fz*fz);
        
        fx = fx/fn;
        fy = fy/fn;
        fz = fz/fn;
        
        psi = p->F45*sqrt(pow(fx*p->DXP[IP],2.0) + pow(fy*p->DYN[JP],2.0) + pow(fz*p->DZN[KP],2.0));
        }
        
        if(bb>0)
        {
        fx = fabs( 0.5*(a->phi(i+1,j,k) + a->phi(i+1,j+1,k)) - 0.5*(a->phi(i-1,j,k)+a->phi(i-1,j+1,k))/(p->DXP[IM1]+p->DXP[IP]));
        fy = fabs((a->phi(i,j+1,k) - a->phi(i,j,k))/p->DYP[JP]);      
        fz = fabs( 0.5*(a->phi(i,j,k+1) + a->phi(i,j+1,k+1)) - 0.5*(a->phi(i,j,k-1)+a->phi(i,j+1,k-1))/(p->DZP[KM1]+p->DZP[KP]));
        
        fn = sqrt(fx*fx + fy*fy + fz*fz);
        
        fx = fx/fn;
        fy = fy/fn;
        fz = fz/fn;
        
        psi = p->F45*sqrt(pow(fx*p->DXP[IP],2.0) + pow(fy*p->DYN[JP],2.0) + pow(fz*p->DZN[KP],2.0));
        }
        
        
        if(cc<0)
        {
        fx = fabs( 0.5*(a->phi(i+1,j,k) + a->phi(i+1,j,k-1)) - 0.5*(a->phi(i-1,j,k)+a->phi(i-1,j,k-1))/(p->DXP[IM1]+p->DXP[IP]));
        fy = fabs( 0.5*(a->phi(i,j+1,k) + a->phi(i,j+1,k-1)) - 0.5*(a->phi(i,j-1,k)+a->phi(i,j-1,k-1))/(p->DYP[JM1]+p->DYP[JP]));      
        fz = fabs((a->phi(i,j,k) - a->phi(i,j,k-1))/p->DZP[KM1]);
        
        fn = sqrt(fx*fx + fy*fy + fz*fz);
        fn= fx+fy+fz;
        
        fx = fx/fn;
        fy = fy/fn;
        fz = fz/fn;
        
        psi = p->F45*sqrt(pow(fx*p->DXP[IP],2.0) + pow(fy*p->DYN[JP],2.0) + pow(fz*p->DZN[KP],2.0));
        }
        
        if(cc>0)
        {
        fx = fabs( 0.5*(a->phi(i+1,j,k) + a->phi(i+1,j,k+1)) - 0.5*(a->phi(i-1,j,k)+a->phi(i-1,j,k+1))/(p->DXP[IM1]+p->DXP[IP]));
        fy = fabs( 0.5*(a->phi(i,j+1,k) + a->phi(i,j+1,k+1)) - 0.5*(a->phi(i,j-1,k)+a->phi(i,j-1,k+1))/(p->DYP[JM1]+p->DYP[JP]));      
        fz = fabs((a->phi(i,j,k+1) - a->phi(i,j,k))/p->DZP[KP]);
        
        fn = sqrt(fx*fx + fy*fy + fz*fz);
        
        fx = fx/fn;
        fy = fy/fn;
        fz = fz/fn;
        
        psi = p->F45*sqrt(pow(fx*p->DXP[IP],2.0) + pow(fy*p->DYN[JP],2.0) + pow(fz*p->DZN[KP],2.0));
        }
        
        if(phival>psi)
        H=1.0;

        if(phival<-psi)
        H=0.0;

        if(fabs(phival)<=psi)
        H=0.5*(1.0 + phival/psi + (1.0/PI)*sin((PI*phival)/psi));
        
            
        roval = p->W1*H + p->W3*(1.0-H);
	}
	
	// -----
	
	if(p->D32==5)
	{
	
	roval = (fabs(a->phi(i,j,k))*a->ro(i,j,k) + fabs(a->phi(i+aa,j+bb,k+cc))*a->ro(i+aa,j+bb,k+cc))
		    /(fabs(a->phi(i,j,k)) + fabs(a->phi(i+aa,j+bb,k+cc)));
	}
	
	// -----
	
	if(p->D32==6)
	{
	double pval1,pval2,H1,H2;
	
	phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
	//phival = (1.0/16.0)*(- a->phi(i-aa,j-bb,k-cc) + 9.0*a->phi(i,j,k) + 9.0*a->phi(i+aa,j+bb,k+cc) - a->phi(i+2*aa,j+2*bb,k+2*cc));
	
	if(phival>0.5*p->dx)
	{
	H1=0.5;
	H2=0.5;
	}

	if(phival<-0.5*p->dx)
	{
	H1=0.5;
	H2=0.5;
	}
	
	if(fabs(phival)<0.5*p->dx)
	{
	H1 = fabs(a->phi(i,j,k))/p->dx;
	H2 = fabs(a->phi(i+aa,j+bb,k+cc))/p->dx;
	}
		
	roval = H1*a->ro(i,j,k) + H2*a->ro(i+aa,j+bb,k+cc);
	}
	
	if(p->D32==7)
	{
		H= 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

		H=MAX(H,0.0);
		H=MIN(H,1.0);

		roval = p->W1*H +   p->W3*(1.0-H);
	}
	
	if(p->D32==8)
	roval = 0.5*(a->ro(i+ii,j+jj,k+kk) + a->ro(i+aa,j+bb,k+cc));
	
	
	if(p->D32==9)
	{
	phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
	
	if(phival>epsi)
	H=1.0;

	if(phival<-epsi)
	H=0.0;

	if(fabs(phival)<=epsi)
	H=0.5*(1.0 + phival/eps + (1.0/PI)*sin((PI*phival)/eps));
	
		
	roval = p->W1*H + p->W3*(1.0-H);
	}
    
    // -----
	if(p->D32==10)
	{
	
	phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));
	
	if(phival>epsi)
	H=1.0;

	if(phival<-epsi)
	H=0.0;

	if(fabs(phival)<=epsi)
	H=0.5*(1.0 + phival/eps + (1.0/PI)*sin((PI*phival)/eps));
	
		
	roval = p->W1*H + p->W3*(1.0-H);
	}
	
	return roval;		
}




