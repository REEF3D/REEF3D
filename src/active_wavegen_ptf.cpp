#include"iowave.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"

void iowave::active_wavegen_ptf(lexer *p, fdm_ptf* e, ghostcell* pgc, field& u, field& v, field& w)
{
    int ii;
    double fac,fac1,epsi,H;
    double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
    
    
        // get the fsf elevation
        LOOP
		wsfmax[i][j]=-1.0e20;

		LOOP
		if(e->phi(i,j,k)>=0.0 && e->phi(i,j,k+1)<0.0)
		wsfmax[i][j]=MAX(wsfmax[i][j],-(e->phi(i,j,k)*p->DZP[KP])/(e->phi(i,j,k+1)-e->phi(i,j,k)) + p->pos_z());

        for(int qn=0; qn<p->mz;++qn)
		pgc->verticalmax_ptf(p,e,wsfmax);
        
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
        
        

			if(e->phi(i-1,j,k)>=0.0)
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

			if(e->phi(i-1,j,k)<0.0 && e->phi(i-1,j,k)>=-p->F45*p->DZP[KP])
			{
			fac= p->B122*(1.0 - fabs(e->phi(i-1,j,k))/(p->F45*p->DZP[KP]));
            
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

			if(e->phi(i-1,j,k)<-p->F45*p->DXM)
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
        
                if(e->phi(i,j,k)>epsi)
                H=1.0;

                if(e->phi(i,j,k)<-epsi)
                H=0.0;

                if(fabs(e->phi(i,j,k))<=epsi)
                H=0.5*(1.0 + e->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*e->phi(i,j,k))/epsi));
                
                
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
        }
        }
        
}

