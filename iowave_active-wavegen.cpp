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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void iowave::active_wavegen(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{		
		double eta_R,Uc,Un,Ucorr,Vc,Wc,eta_T,eta_M,wsf;
		double posx,posy,posz,uvel,vvel,uabs,fx,fy,pval;
		x=0.0;
		z=0.0;
		double fac1;
		int aa,bb,ii,jj;
		
		LOOP
		wsfmax[i][j]=-1.0e20;

		LOOP
		if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
		wsfmax[i][j]=MAX(wsfmax[i][j],-(a->phi(i,j,k)*p->DXM)/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());

        for(int qn=0; qn<p->mz;++qn)
		pgc->verticalmax(p,a,wsfmax);
	
	// UVEL
		for(n=0;n<gcgen1_count;++n)
		{
		i=gcgen1[n][0];
		j=gcgen1[n][1];
				
		
		ii=0;
		if(gcgen1[n][3]==4)
		ii=1;
		
		wsf=wsfmax[i+ii][j];
        
        //cout<<p->mpirank<<" WSF: "<<setprecision(5)<<wsf<<endl;


		eta_T = wave_h(p,pgc,x,0.0,0.0)-p->wd;
		eta_M = wsf-p->wd; 
		eta_R = -eta_M+eta_T;
		
		
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
		
		
		uvel=uvel/double(count);
		vvel=vvel/double(count);
		uabs = sqrt(uvel*uvel + vvel*vvel);
		
		
		Ucorr = sqrt(9.81/p->wd)*eta_R; 
		
		Uc=fx=fy=0.0;
		
		if(gcgen1[n][2]==1)
		{
		Uc=-Ucorr*(uvel/(uabs>1.0e-10?uabs:1.0e20));
		fx=1.0;
		}
		
		if(gcgen1[n][2]==4)
		{
		Uc=Ucorr*(uvel/(uabs>1.0e-10?uabs:1.0e20));
		fx=1.0;
		}
		
		
		Uc=Ucorr;

		if(Uc>=0.0)
		fac1=1.0;
		
		if(Uc<0.0)
		fac1=0.0;
		
		aa=bb=0;
		if(gcgen1[n][2]==1)
		aa=-1;
		
		if(gcgen1[n][2]==4)
		aa=1;
		
		if(gcgen1[n][2]==3)
		bb=-1;
		
		if(gcgen1[n][2]==2)
		bb=1;


           KLOOP
           PCHECK 
			{
				if(a->phi(i-1,j,k)>=0.0)
				{
					if(p->pos_z()<=p->phimean)
                    {
                    z=-(fabs(p->phimean-p->pos_z()));
                    z3=z+0.5*p->DXM;
                    }
                    
                    if(p->pos_z()>p->phimean)
                    {
                    z=(fabs(p->phimean-p->pos_z()));
                    z3=z;
                    }

					uvel=wave_u(p,pgc,x,0.0,z);
					wvel=wave_w(p,pgc,x,0.0,z3);
				}
				
				if(a->phi(i-1,j,k)<0.0)
				{
				z = wave_h(p,pgc,xg,0.0,0.0)-p->pos_z();	
                z3 = z;
					
					uvel=wave_u(p,pgc,x,0.0,z);
					wvel=wave_w(p,pgc,x,0.0,z3);
				}
					
				if(p->B98==4)
				Uc=eta_R*sqrt(9.81/p->wd);//*fabs((uvel/(uabs>1.0e-10?uabs:1.0e20)));
				
				if(p->B98==5)
				Uc=eta_R*p->ww*( cosh(p->wk*(p->wd+z))/sinh(p->wk*p->wd));//*fabs((uvel/(uabs>1.0e-10?uabs:1.0e20)));
				
			
				if(z<=eta_M)
				{
				u(i+1*aa,j+1*bb,k)=uvel + p->Ui + Uc*fx;
				u(i+2*aa,j+2*bb,k)=uvel + p->Ui + Uc*fx;
				u(i+3*aa,j+3*bb,k)=uvel + p->Ui + Uc*fx;
				
				w(i-1,j,k)=wvel;
				w(i-2,j,k)=wvel;
				w(i-3,j,k)=wvel;
				}

				if(z>=eta_M && z<eta_M+p->F45*p->DXM)
				{
				fac=1.0 - fabs(a->phi(i-1,j,k))/p->F45*p->DXM;
				
				u(i+1*aa,j+1*bb,k)=uvel*fac + p->Ui*fac + Uc*fx*fac*fac1;
				u(i+2*aa,j+2*bb,k)=uvel*fac + p->Ui*fac + Uc*fx*fac*fac1;
				u(i+3*aa,j+3*bb,k)=uvel*fac + p->Ui*fac + Uc*fx*fac*fac1;
				
				w(i-1,j,k)=wvel*fac;
				w(i-2,j,k)=wvel*fac;
				w(i-3,j,k)=wvel*fac;
				}
             
    
				if(z>=eta_M+p->F45*p->DXM)
				{
				u(i+1*aa,j+1*bb,k)=0.0;
				u(i+2*aa,j+2*bb,k)=0.0;
				u(i+3*aa,j+3*bb,k)=0.0;
				
				w(i-1,j,k)=0.0;
				w(i-2,j,k)=0.0;
				w(i-3,j,k)=0.0;
				}
			}
		}
	
        // NSEWAVE
        if(p->A10==5)
        for(n=0;n<gcgen1_count;++n)
		{
		i=gcgen1[n][0];
		j=gcgen1[n][1];
        
        x=xgen(p);
        y=ygen(p);
        a->eta(i-1,j)=wave_eta(p,pgc,x,y);
        a->eta(i-2,j)=wave_eta(p,pgc,x,y);
        a->eta(i-3,j)=wave_eta(p,pgc,x,y);
        }
}

void iowave::active_wavegen2(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{	

		int q;
		double eta_M, eta_T;
		double x=0.0;
		double z=0.0;
		double fac1,fac2,fac3;
		double z1,z2;
		uvel = 0.0;
		wvel = 0.0;
		
		eta_T = wave_h(p,pgc,x,0.0,0.0)-p->phimean;
		
		eta_M = p->phiin-p->wd;
		
		if(eta_M>=eta_T)
		{
		fac1 = 1.0;
		fac2 = 0.0;
		fac3 = 0.0;
		
		z1 = eta_T;
		z2 = eta_M;
		}
		
		if(eta_T>eta_M)
		{
		fac1 = 1.0;
		fac2 = 1.0;
		fac3 = 1.0;

		z1 = eta_M;
		z2 = eta_T;
		}
				
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];

                if(a->phi(i-1,j,k)>=0.0)
				{
					if(p->pos_z()<=p->phimean)
                    {
                    z=-(fabs(p->phimean-p->pos_z()));
                    z3=z+0.5*p->DXM;
                    }
                    
                    if(p->pos_z()>p->phimean)
                    {
                    z=(fabs(p->phimean-p->pos_z()));
                    z3=z;
                    }

					uvel=wave_u(p,pgc,x,0.0,z);
					wvel=wave_w(p,pgc,x,0.0,z3);
				}
				
				if(a->phi(i-1,j,k)<0.0)
				{
				z = wave_h(p,pgc,xg,0.0,0.0)-p->pos_z();	
                z3 = z;
					
					uvel=wave_u(p,pgc,x,0.0,z);
					wvel=wave_w(p,pgc,x,0.0,z3);
				}
		
			if(z<z1)
			{
			u(i-1,j,k)=uvel+p->Ui;
			u(i-2,j,k)=uvel+p->Ui;
			u(i-3,j,k)=uvel+p->Ui;
			
			w(i-1,j,k)=wvel;
			w(i-2,j,k)=wvel;
			w(i-3,j,k)=wvel;
			}
			
			if(z>=z1 && z<z2)
			{
			u(i-1,j,k)=uvel*fac2 + p->Ui;
			u(i-2,j,k)=uvel*fac2 + p->Ui;
			u(i-3,j,k)=uvel*fac2 + p->Ui;
			
			w(i-1,j,k)=wvel*fac2;
			w(i-2,j,k)=wvel*fac2;
			w(i-3,j,k)=wvel*fac2;
			}

			if(z>=z2 && z<z2+p->F45*p->DXM)
			{
			fac=p->B122*(1.0 - fabs(a->phi(i-1,j,k))/p->F45*p->DXM);
			u(i-1,j,k)=uvel*fac*fac3 + p->Ui;
			u(i-2,j,k)=uvel*fac*fac3 + p->Ui;
			u(i-3,j,k)=uvel*fac*fac3 + p->Ui;

			w(i-1,j,k)=wvel*fac*fac3;
			w(i-2,j,k)=wvel*fac*fac3;
			w(i-3,j,k)=wvel*fac*fac3;
			}

			if(z>=z2+p->F45*p->DXM)
			{
			pgc->dirichlet_ortho(p,u,p->DXM,10,1,1);

			w(i-1,j,k)=0.0;
			w(i-2,j,k)=0.0;
			w(i-3,j,k)=0.0;
			}
		}	
		
		
		if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
		if(p->B64==1)
		{
		for(int q=0;q<4;++q)
		for(n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0]+q;
		j=p->gcin[n][1];
		k=p->gcin[n][2];

		if(a->phi(i,j,k)<0.0)
		a->eddyv(i,j,k)=MIN(a->eddyv(i,j,k),1.0e-4);
		}
		pgc->start4(p,a->eddyv,24);
		}
        
        // NSEWAVE
        if(p->A10==5)
        for(n=0;n<gcgen1_count;++n)
		{
		i=gcgen1[n][0];
		j=gcgen1[n][1];
        
        x=xgen(p);
        y=ygen(p);
        a->eta(i-1,j)=wave_eta(p,pgc,x,y);
        a->eta(i-2,j)=wave_eta(p,pgc,x,y);
        a->eta(i-3,j)=wave_eta(p,pgc,x,y);
        }

}
