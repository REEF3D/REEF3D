/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void iowave::active_beach_fnpf(lexer *p, ghostcell* pgc, double *Fi, double *Uout, slice &Fifsf, slice &eta)
{
    /*
		double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
		double posx,posy,posz,uvel,vvel,uabs,fx,fy,pval,fp;
        double fxdir,fydir;
		double x=0.0;
		double z=0.0;
		double fac1,fac,multiplier;
		int aa,bb,ii,jj;
        
		
		// Fi
		for(n=0;n<gcawa4_count;n++)
		{
		i=gcawa4[n][0];
		j=gcawa4[n][1];
		
		wsf=wsfmax[i+ii][j];
		
		
		eta_T = 0.0;
		eta_M = wsf-p->wd; 
		eta_R = eta_M-eta_T;
        
        //cout<<p->mpirank<<" eta_R: "<<eta_R<<" eta_M: "<<eta_M<<"   wsf: "<<wsf<<endl;
		
		count=0;
		uvel=vvel=0.0;
		KLOOP
        PCHECK
        {
			posx=p->XN[IP1];
			posy=p->YP[JP];
			posz=p->ZP[KP];

             uvel+=p->ccipol1(a->u,posx,posy,posz);
			  vvel+=p->ccipol2(a->v,posx,posy,posz);
			
			++count;
        }
        
        // vertical sum uvel and vvel !
		
		uvel=uvel/double(count);
		vvel=vvel/double(count);
		uabs = sqrt(uvel*uvel + vvel*vvel);
		

		Uc=fx=fy=0.0;
		
        
        fx=(uvel/(uabs>1.0e-10?uabs:1.0e20));
        

		if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;
		
		aa=bb=0;
		if(gcawa1[n][2]==1)
        {
		aa=-1;
		fx*=-1.0;
        }
        
		if(gcawa1[n][2]==4)
		aa=1;
		
		if(gcawa1[n][2]==3)
        {
		bb=-1;
        fx*=-0.0;
        }
		
		if(gcawa1[n][2]==2)
        {
		bb=1;
        fx*=-1.0;
        }
        
        fx=1.0;

			if(wsf>-1.0e20)
			KLOOP 
			{
				if(p->pos_z()<=p->phimean)
				z=-(fabs(p->phimean-p->pos_z()));
				
				if(p->pos_z()>p->phimean)
				z=(fabs(p->phimean-p->pos_z()));
				
				if(p->B99==3)
				Uc=eta_R*sqrt(9.81/p->wd);
				
				if(p->B99==4)
				Uc=eta_R*p->ww*( cosh(p->wk*(p->wd+z))/sinh(p->wk*p->wd));
                
               if(p->B99==5)
               {
                   if(p->pos_z()>p->B123)
                   {
                   fac = (p->pos_z()-p->B123)/(wsf-p->B123);
                   multiplier = 2.0*((wsf)/(wsf-p->B123));
                   Uc =   multiplier*fac*eta_R*sqrt(9.81/p->wd);
                   }
                   
                   if(p->pos_z()<=p->B123)
                   Uc=0.0;
               }
                   

				if(z<=eta_M)
				{
				u(i+1*aa,j+1*bb,k)=Uc*fx;
				u(i+2*aa,j+2*bb,k)=Uc*fx;
				u(i+3*aa,j+3*bb,k)=Uc*fx;
				}

				if(z>=eta_M && z<eta_M+p->F45*p->DZP[KP])
				{
				fac=p->B122*(1.0 - fabs(a->phi(i-1,j,k))/p->F45*p->DZP[KP]);
				
				u(i+1*aa,j+1*bb,k)=Uc*fx*fac*fac1;
				u(i+2*aa,j+2*bb,k)=Uc*fx*fac*fac1;
				u(i+3*aa,j+3*bb,k)=Uc*fx*fac*fac1;
				}

				if(z>=eta_M+p->F45*p->DZP[KP])
				{
				u(i+1*aa,j+1*bb,k)=0.0;
				u(i+2*aa,j+2*bb,k)=0.0;
				u(i+3*aa,j+3*bb,k)=0.0;
				}
			}
		}

		
		//-----------------------------------------------		
		// PRESSURE
		for(n=0;n<gcawa4_count;n++)
		{
		i=gcawa4[n][0];
		j=gcawa4[n][1];
		
		aa=bb=0;
		
		if(gcawa4[n][2]==1)
		aa=-1;
		
		if(gcawa4[n][2]==4)
		aa=1;
		
		if(gcawa4[n][2]==3)
		bb=-1;
		
		if(gcawa4[n][2]==2)
		bb=1;
				
		wsf=wsfmax[i][j];
        
        
        
        eta_T = 0.0;
        eta_M = wsf-p->wd; 
        eta_R = fabs(eta_M-eta_T);
        
        double r=0.0;
        
        double wH=0.25*p->wH;
    
        if(p->B92>30)
        wH=0.25*p->wHs;
        
        wH=MAX(wH,0.5*p->DXM);

        x=fabs(eta_R/(fabs(wH)>1.0e-20?wH:1.0e20));
        x=MIN(x,1.0);
        
        r = -2.0*x*x*x + 3.0*x*x;
    
        //r=1.0;

			//if(wsf>-1.0e20)
			KLOOP 
			{
                            
                
			if(p->B78==1)
			pval=(wsf - p->pos_z()+0.5*p->DZP[KP])*a->ro(i,j,k)*fabs(p->W22);
			
			if(p->B78==2)
			pval=(p->wd - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
			
			if(p->B78==3)
			pval=a->press(i,j,k);
			
			a->press(i+1*aa,j+1*bb,k)=r*pval + (1.0-r)*a->press(i,j,k);
			a->press(i+2*aa,j+2*bb,k)=r*pval + (1.0-r)*a->press(i,j,k);
			a->press(i+3*aa,j+3*bb,k)=r*pval + (1.0-r)*a->press(i,j,k);
			
			a->w(i+1*aa,j+1*bb,k)=0.0;
			a->w(i+2*aa,j+2*bb,k)=0.0;
			a->w(i+3*aa,j+3*bb,k)=0.0;
			}		
		}
        
        
        */
      
}
