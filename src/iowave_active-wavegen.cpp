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

void iowave::active_wavegen(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    int ii;
    double fac,fac1,epsi,H;
    double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
    
    
        // get the fsf elevation
        LOOP
		wsfmax[i][j]=-1.0e20;

		LOOP
		if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
		wsfmax[i][j]=MAX(wsfmax[i][j],-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());

        for(int qn=0; qn<p->mz;++qn)
		pgc->verticalmax(p,a,wsfmax);
        
        // wavegen
        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];		

        uvel=uval[count]*ramp(p);
        vvel=vval[count]*ramp(p);
        wvel=wval[count]*ramp(p);
        
        

			if(a->phi(i-1,j,k)>=0.0)
			{
            
			u(i-1,j,k)=uvel+p->Ui;
			u(i-2,j,k)=uvel+p->Ui;
			u(i-3,j,k)=uvel+p->Ui;
            
             v(i-1,j,k)=vvel;
			v(i-2,j,k)=vvel;
			v(i-3,j,k)=vvel;
			
			w(i-1,j,k)=wvel;
			w(i-2,j,k)=wvel;
			w(i-3,j,k)=wvel;
			}

			if(a->phi(i-1,j,k)<0.0 && a->phi(i-1,j,k)>=-p->F45*p->DZP[KP])
			{
			fac= p->B122*(1.0 - fabs(a->phi(i-1,j,k))/(p->F45*p->DZP[KP]));
            
			u(i-1,j,k)=uvel*fac + p->Ui;
			u(i-2,j,k)=uvel*fac + p->Ui;
			u(i-3,j,k)=uvel*fac + p->Ui;
            
             v(i-1,j,k)=vvel*fac;
			v(i-2,j,k)=vvel*fac;
			v(i-3,j,k)=vvel*fac;
            
			w(i-1,j,k)=wvel*fac;
			w(i-2,j,k)=wvel*fac;
			w(i-3,j,k)=wvel*fac;
			}

			if(a->phi(i-1,j,k)<-p->F45*p->DXM)
			{
			//pgc->dirichlet_ortho(p,u,p->DXM,10,1,1);
			u(i-1,j,k)=0.0 + p->Ui;
			u(i-2,j,k)=0.0 + p->Ui;
			u(i-3,j,k)=0.0 + p->Ui;
            
            v(i-1,j,k)=0.0;
			v(i-2,j,k)=0.0;
			v(i-3,j,k)=0.0;

			w(i-1,j,k)=0.0;
			w(i-2,j,k)=0.0;
			w(i-3,j,k)=0.0;
			}
            
                
                // fsf deviation
                
                wsf=wsfmax[i][j];
                
                eta_T = wave_eta(p,pgc,x,0.0);
                eta_M = wsf-p->wd; 
                eta_R = eta_T-eta_M;
                
                
                if(eta_R>=0.0)
                fac1=1.0;
                
                if(eta_R<0.0)
                fac1=0.0;
        
        
                if(p->pos_z()<=p->phimean)
				z=-(fabs(p->phimean-p->pos_z()));
				
				if(p->pos_z()>p->phimean)
				z=(fabs(p->phimean-p->pos_z()));
				
				if(p->B98==4)
				Uc=eta_R*sqrt(9.81/p->wd);
                
                
                
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
				u(i-1,j,k)+=Uc;
				u(i-2,j,k)+=Uc;
				u(i-3,j,k)+=Uc;
				}

				if(z>=eta_M && z<eta_M+epsi)
				{
				u(i-1,j,k)+=Uc*H*fac1;
				u(i-2,j,k)+=Uc*H*fac1;
				u(i-3,j,k)+=Uc*H*fac1;
				}

				if(z>=eta_M+epsi)
				{
				u(i-1,j,k)+=0.0;
				u(i-2,j,k)+=0.0;
				u(i-3,j,k)+=0.0;
				}
            
            
            
        ++count;
		}
        
        
        if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
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
    if(p->A10==55)
    {
        for(n=0;n<gcgen4_count;++n)
		{
		i=gcgen4[n][0];
		j=gcgen4[n][1];
        
        x=xgen(p);
        y=ygen(p);
        
        a->eta(i-1,j)=a->eta(i-2,j)=a->eta(i-3,j)=wave_eta(p,pgc,x,y);        
            
            KLOOP
            {
            a->phi(i-1,j,k) = a->eta(i-1,j) + p->phimean - p->pos_z();
            a->phi(i-2,j,k) = a->eta(i-2,j) + p->phimean - p->pos_z();
            a->phi(i-3,j,k) = a->eta(i-3,j) + p->phimean - p->pos_z();
            }
        }
        
        

        for(n=0;n<gcgen1_count;++n)
		{
		i=gcgen1[n][0];
		j=gcgen1[n][1];
            
            a->P(i-1,j)=0.0;
            double d=0.0;
            double epsi=p->A440*p->DXM;
            
            KULOOP
            {
            phival = 0.5*(a->phi(i,j,k)+a->phi(i-1,j,k));
            
                if(phival>epsi)
                H=1.0;

                if(phival<-epsi)
                H=0.0;

                if(fabs(phival)<=epsi)
                H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
                
                a->P(i-1,j) += a->u(i-1,j,k)*p->DXM*H;
                d+=p->DXM*H;
            }
            a->P(i-1,j)/=d;
            a->P(i-3,j)=a->P(i-2,j)=a->P(i-1,j);

        }
        
    }
}

