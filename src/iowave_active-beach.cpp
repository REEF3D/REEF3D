/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"fdm.h"
#include"ghostcell.h"

void iowave::active_beach(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
		double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
		double posx,posy,posz,uvel,vvel,uabs,fx,fy,pval,fp;
        double fxdir,fydir;
		double x=0.0;
		double z=0.0;
		double fac1,fac,multiplier;
		int aa,bb,ii,jj;
        double H,epsi;
        
        // get the fsf elevation
        LOOP
		wsfmax[i][j]=-1.0e20;

		LOOP
		if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
		wsfmax[i][j]=MAX(wsfmax[i][j],-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());

        for(int qn=0; qn<p->mz;++qn)
		pgc->verticalmax(p,a,wsfmax);
        
		
		// UVEL
		for(n=0;n<p->gcslawa1_count;++n)
		{
		i=p->gcslawa1[n][0];
		j=p->gcslawa1[n][1];
		
		ii=0;
		if(p->gcslawa1[n][2]==4)
		ii=1;
		
		wsf=wsfmax[i+ii][j];
        
		
		eta_T = 0.0;
		eta_M = wsf-p->wd; 
		eta_R = eta_M-eta_T;
        
        if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;
        //cout<<p->mpirank<<" eta_R: "<<eta_R<<" eta_M: "<<eta_M<<"   wsf: "<<wsf<<endl;
		
        aa=bb=0;
		if(p->gcslawa1[n][2]==1)
		aa=-1;

		if(p->gcslawa1[n][2]==4)
		aa=1;
		
		if(p->gcslawa1[n][2]==3)
		bb=-1;
 
		if(p->gcslawa1[n][2]==2)
		bb=1;


        
			if(wsf>-1.0e19)
			KLOOP 
           PCHECK
			{
                //cout<<p->mpirank<<" eta_R: "<<eta_R<<" eta_M: "<<eta_M<<"   wsf: "<<wsf<<endl;
                
				if(p->pos_z()<=p->phimean)
				z=-(fabs(p->phimean-p->pos_z()));
				
				if(p->pos_z()>p->phimean)
				z=(fabs(p->phimean-p->pos_z()));
				
				if(p->B99==3)
				Uc=eta_R*sqrt(9.81/p->wd);
				
				if(p->B99==4)
				Uc=eta_R*p->ww*(cosh(p->wk*(p->wd+z))/sinh(p->wk*p->wd));
                
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
               
               // inteface H
                epsi = p->F45*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
                if(a->phi(i,j,k)>epsi)
                H=1.0;

                if(a->phi(i,j,k)<-epsi)
                H=0.0;

                if(fabs(a->phi(i,j,k))<=epsi)
                H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));
                

				if(z<=eta_M)
				{
				u(i+1*aa,j+1*bb,k)=Uc;
				u(i+2*aa,j+2*bb,k)=Uc;
				u(i+3*aa,j+3*bb,k)=Uc;
				}

				if(z>=eta_M && z<eta_M+epsi)
				{
				u(i+1*aa,j+1*bb,k)=Uc*H*fac1;
				u(i+2*aa,j+2*bb,k)=Uc*H*fac1;
				u(i+3*aa,j+3*bb,k)=Uc*H*fac1;
				}

				if(z>=eta_M+epsi)
				{
				u(i+1*aa,j+1*bb,k)=0.0;
				u(i+2*aa,j+2*bb,k)=0.0;
				u(i+3*aa,j+3*bb,k)=0.0;
				}
			}
            
            if(wsf<-1.0e19)
			KLOOP
            PCHECK 
			{
            u(i+1*aa,j+1*bb,k)=0.0;
            u(i+2*aa,j+2*bb,k)=0.0;
            u(i+3*aa,j+3*bb,k)=0.0;
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
		
		wsf=wsfmax[i][j+jj];

		eta_T = 0.0;
		eta_M = wsf-p->wd; 
		eta_R = eta_M-eta_T;
        
        if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;
		
		
        /*
		count=0;
		uvel=vvel=0.0;
		KLOOP
        PCHECK
        {	
			posx=p->XP[IP];
			posy=p->YN[JP1];
			posz=p->ZP[KP];

            uvel=p->ccipol1(a->u,posx,posy,posz);
			 vvel=p->ccipol2(a->v,posx,posy,posz);
			
			++count;
        }
		
		uvel=uvel/double(count);
		vvel=vvel/double(count);
		uabs = sqrt(uvel*uvel + vvel*vvel);
		
        
        fy=(vvel/(uabs>1.0e-10?uabs:1.0e20));
			
		if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;
		
		aa=bb=0;
		if(p->gcslawa2[n][3]==1)
        {
		aa=-1;
        fy*=1.0;
        }
		
		if(p->gcslawa2[n][3]==4)
		aa=1;
		
		if(p->gcslawa2[n][3]==3)
        {
		bb=-1;
        fy*=-1.0;
        }
		
		if(p->gcslawa2[n][3]==2)
        {
		bb=1;
        fy*=-1.0;
        }*/
        
        aa=bb=0;
		if(p->gcslawa2[n][2]==1)
		aa=-1;
		
		if(p->gcslawa2[n][2]==4)
		aa=1;
		
		if(p->gcslawa2[n][2]==3)
		bb=-1;
		
		if(p->gcslawa2[n][2]==2)
		bb=1;
 
        
        fy=0.0; // !
 


			if(wsf>-1.0e19)
			KLOOP
            PCHECK
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
               
               // inteface H
                epsi = p->F45*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
                if(a->phi(i,j,k)>epsi)
                H=1.0;

                if(a->phi(i,j,k)<-epsi)
                H=0.0;

                if(fabs(a->phi(i,j,k))<=epsi)
                H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));
                
			
				if(z<=eta_M)
				{
				v(i+1*aa,j+1*bb,k)=Uc*fy;
				v(i+2*aa,j+2*bb,k)=Uc*fy;
				v(i+3*aa,j+3*bb,k)=Uc*fy;
				}

				if(z>=eta_M && z<eta_M+epsi)
				{
				v(i+1*aa,j+1*bb,k)=Uc*fy*H*fac1;
				v(i+2*aa,j+2*bb,k)=Uc*fy*H*fac1;
				v(i+3*aa,j+3*bb,k)=Uc*fy*H*fac1;
				}

				if(z>=eta_M+epsi)
				{
				v(i+1*aa,j+1*bb,k)=0.0;
				v(i+2*aa,j+2*bb,k)=0.0;
				v(i+3*aa,j+3*bb,k)=0.0;
				}
			}
            
            if(wsf>-1.0e19)
			KLOOP 
			{
            v(i+1*aa,j+1*bb,k)=0.0;
            v(i+2*aa,j+2*bb,k)=0.0;
            v(i+3*aa,j+3*bb,k)=0.0;
            }
		}
		
		//-----------------------------------------------		
		// PRESSURE
    
		for(n=0;n<p->gcslout_count;n++)
		{
		i=p->gcslout[n][0];
		j=p->gcslout[n][1];
		
		aa=bb=0;
		
		if(p->gcslout[n][3]==1)
		aa=-1;
		
		if(p->gcslout[n][3]==4)
		aa=1;
		
		if(p->gcslout[n][3]==3)
		bb=-1;
		
		if(p->gcslout[n][3]==2)
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
        
        H = -2.0*x*x*x + 3.0*x*x;
    
        /*
                epsi = p->F45*p->DZN[KP];
        
                if(a->phi(i,j,k)>epsi)
                H=1.0;

                if(a->phi(i,j,k)<-epsi)
                H=0.0;

                if(fabs(a->phi(i,j,k))<=epsi)
                H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));*/

		//cout<<p->mpirank<<" eta_R: "<<eta_R<<" x: "<<x<<" r: "<<r<<endl;
		
            if(wsf>-1.0e19)
			KLOOP 
            PCHECK
			{
			pval=(wsf - p->pos_z()+0.5*p->DZP[KP])*a->ro(i,j,k)*fabs(p->W22);
			

			a->press(i+1*aa,j+1*bb,k)=H*pval + (1.0-H)*a->press(i,j,k);
			a->press(i+2*aa,j+2*bb,k)=H*pval + (1.0-H)*a->press(i,j,k);
			a->press(i+3*aa,j+3*bb,k)=H*pval + (1.0-H)*a->press(i,j,k);
            
            /*
            if(p->B30==2)
            {
            a->press(i+1*aa,j+1*bb,k) = a->press(i,j,k);
            a->press(i+2*aa,j+2*bb,k) = a->press(i,j,k);
            a->press(i+3*aa,j+3*bb,k) = a->press(i,j,k);
            }*/
			
			w(i+1*aa,j+1*bb,k)=0.0;
			w(i+2*aa,j+2*bb,k)=0.0;
			w(i+3*aa,j+3*bb,k)=0.0;
			}	
		}
        
        
        
        // NSEWAVE
        if(p->A10==55)
        for(n=0;n<p->gcslawa1_count;++n)
		{
		i=p->gcslawa1[n][0];
		j=p->gcslawa1[n][1];
            
            a->P(i+1,j)=0.0;
            double d=0.0;
            double epsi=p->A440*p->DXM;
            KULOOP
            {
            phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
            
                if(phival>epsi)
                H=1.0;

                if(phival<-epsi)
                H=0.0;

                if(fabs(phival)<=epsi)
                H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
                
                a->P(i+1,j) += u(i+1,j,k)*p->DZP[KP]*H;
                d+=p->DXM*H;
            }
            a->P(i+1,j)/=d;
        }
}
