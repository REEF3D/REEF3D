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

#include"iowave.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void iowave::active_beach2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
	double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M;
    double posx,posy,posz,uvel,vvel,uabs,fx,fy,pval,fp;
    double fxdir,fydir,dP,dQ;
    double x=0.0;
    double z=0.0;
    double fac1,fac,multiplier;
    int aa,bb,ii,jj;
        
        
		
    // UVEL
    for(n=0;n<p->gcslawa1_count;++n)
    {
		i=p->gcslawa1[n][0];
		j=p->gcslawa1[n][1];
		
		ii=0;
		if(p->gcslawa1[n][3]==4)
		ii=1;
			
		eta_T = 0.0;
		eta_M = eta(i+ii,j);
		eta_R = eta_M-eta_T;
        
		uvel = b->P(i,j);
		vvel = 0.5*(b->Q(i,j)+b->Q(i,j-1));
		uabs = sqrt(uvel*uvel + vvel*vvel);
		

		Uc=fx=fy=0.0;
		
        fx=(uvel/(uabs>1.0e-10?uabs:1.0e20));
        

		if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;
		
		aa=bb=0;
		if(p->gcslawa1[n][2]==1)
        {
		aa=-1;
		fx*=-1.0;
        }
        
		if(p->gcslawa1[n][2]==4)
		aa=1;
		
		if(p->gcslawa1[n][2]==3)
        {
		bb=-1;
        fx*=-0.0;
        }
		
		if(p->gcslawa1[n][2]==2)
        {
		bb=1;
        fx*=-1.0;
        }
        
        fx=1.0;
        
    if(p->B99==3)
    {    
        Uc = eta_R*sqrt(9.81/p->wd);

        P(i+1*aa,j+1*bb) = Uc*fx;
        P(i+2*aa,j+2*bb) = Uc*fx;
        P(i+3*aa,j+3*bb) = Uc*fx;
     }
    
    
  if(p->B99==4)
  {
    double dfx1,dfx4,dfy2,dfy3;
    
    dfx1 = (P(i+1,j)-P(i,j))/p->DXM;
    dfx4 = (P(i,j)-P(i-1,j))/p->DXM;
    dfy2 = (P(i,j)-P(i,j-1))/p->DXM;
    dfy3 = (P(i,j+1)-P(i,j))/p->DXM;

     /*       
	if(cs==1)
	for(q=0;q<margin;++q)
	P(i-q-1,j) = P(i,j) - p->dt*sqrt(9.81*p->wd)*dfx1;

	if(cs==2)
	for(q=0;q<margin;++q)
	P(i,j+q+1) = P(i,j) - p->dt*sqrt(9.81*p->wd)*dfy2;

	if(cs==3)
	for(q=0;q<margin;++q)
	P(i,j-q-1) = P(i,j) - p->dt*sqrt(9.81*p->wd)*dfy3;
*/
    
	for(int q=1;q<=3;++q)
	P(i+q*aa,j+q*bb) = P(i,j) - p->dt*sqrt(9.81*b->hp(i,j))*dfx4;
        
    }
    }
		
//-----------------------------------------------		
		// VVEL
		
		for(n=0;n<p->gcslawa2_count;++n)
		{
		i=p->gcslawa2[n][0];
		j=p->gcslawa2[n][1];
		
		
		jj=0;
		if(p->gcslawa2[n][3]==2)
		jj=1;
		

		eta_T = 0.0;
		eta_M = b->eta(i,j+jj);
		eta_R = eta_M-eta_T;
		
		
		uvel=0.5*(b->P(i,j)+b->P(i-1,j));
		vvel=b->Q(i,j);
		uabs = sqrt(uvel*uvel + vvel*vvel);
		
        
        fy=(vvel/(uabs>1.0e-10?uabs:1.0e20));
			
		if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;
		
		aa=bb=0;
		if(p->gcslawa2[n][2]==1)
        {
		aa=-1;
        fy*=1.0;
        }
		
		if(p->gcslawa2[n][2]==4)
		aa=1;
		
		if(p->gcslawa2[n][2]==3)
        {
		bb=-1;
        fy*=-1.0;
        }
		
		if(p->gcslawa2[n][2]==2)
        {
		bb=1;
        fy*=-1.0;
        }
        
        fy=0.0; // !

        
            
				Uc=eta_R*sqrt(9.81/p->wd);
				
				b->Q(i+1*aa,j+1*bb) = Uc*fy;
				b->Q(i+2*aa,j+2*bb) = Uc*fy;
				b->Q(i+3*aa,j+3*bb) = Uc*fy;
		}
		
}
