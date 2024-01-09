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
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_active_beach(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
		double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
		double posx,posy,posz,uvel,vvel,uabs,fx,fy,pval,fp;
        double fxdir,fydir;
		double x=0.0;
		double z=0.0;
		double fac1,fac,multiplier;
        int aa,bb;
    
		
		// UVEL
		for(n=0;n<p->gcslout_count;++n)
		{
		i=p->gcslout[n][0];
		j=p->gcslout[n][1];
		
		eta_T = 0.0;
		eta_M = d->eta(i,j); 
		eta_R = eta_M-eta_T;
        
        aa=bb=0;
		if(p->gcslout[n][3]==1)
		aa=-1;
		
		if(p->gcslout[n][3]==4)
		aa=1;
		
		if(p->gcslout[n][3]==3)
		bb=-1;
		
		if(p->gcslout[n][3]==2)
		bb=1;
        
        if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;

        fx=1.0;
            
            if(p->wet[IJ]==1)
			KLOOP
			{
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
                   if(p->ZP[KP]>p->B123)
                   {
                   fac = (p->pos_z()-p->B123)/(wsf-p->B123);
                   multiplier = 2.0*((wsf)/(wsf-p->B123));
                   Uc =   multiplier*fac*eta_R*sqrt(9.81/p->wd);
                   }
                   
                   if(p->ZP[KP]<=p->B123)
                   Uc=0.0;
               }
               
                U[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = Uc*fx;
                U[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = Uc*fx;
                U[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = Uc*fx;
                
                UH[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = (d->eta(i,j)+d->depth(i,j))*Uc*fx;
                UH[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = (d->eta(i,j)+d->depth(i,j))*Uc*fx;
                UH[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = (d->eta(i,j)+d->depth(i,j))*Uc*fx;
                
			}
            
             if(p->wet[IJ]==0)
			KLOOP 
			{
             U[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = 0.0;
             U[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = 0.0;
             U[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = 0.0;
             
             UH[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = 0.0;
             UH[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = 0.0;
             UH[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = 0.0;
            }
		}
//-----------------------------------------------		
		// VVEL
		
		for(n=0;n<p->gcslout_count;++n)
		{
		i=p->gcslout[n][0];
		j=p->gcslout[n][1];
		
		eta_T = 0.0;
		eta_M = d->eta(i,j); 
		eta_R = eta_M-eta_T;
        
        if(eta_R>=0.0)
		fac1=1.0;
		
		if(eta_R<0.0)
		fac1=0.0;
		
		
        aa=bb=0;
		if(p->gcslout[n][3]==1)
		aa=-1;
		
		if(p->gcslout[n][3]==4)
		aa=1;
		
		if(p->gcslout[n][3]==3)
		bb=-1;
		
		if(p->gcslout[n][3]==2)
		bb=1;
 
        
        fy=0.0; // !
 


			if(p->wet[IJ]==1)
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
               
                V[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = Uc*fy;
                V[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = Uc*fy;
                V[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = Uc*fy;
                
                VH[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = (d->eta(i,j)+d->depth(i,j))*Uc*fy;
                VH[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = (d->eta(i,j)+d->depth(i,j))*Uc*fy;
                VH[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = (d->eta(i,j)+d->depth(i,j))*Uc*fy;
				
			}
            
            if(p->wet[IJ]==0)
			KLOOP 
			{
             V[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = 0.0;
             V[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = 0.0;
             V[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = 0.0;
             
             VH[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = 0.0;
             VH[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = 0.0;
             VH[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = 0.0;
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
                    
        
        eta_T = 0.0;
        eta_M = d->eta(i,j); 
        eta_R = fabs(eta_M-eta_T);
        

            if(p->wet[IJ]==1)
			KLOOP 
			{
			d->P[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = 0.0;
             d->P[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = 0.0;
             d->P[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = 0.0;
			
			W[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = 0.0;
             W[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = 0.0;
             W[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = 0.0;
             
             WH[(i-p->imin+1*aa)*p->jmax*p->kmax + (j-p->jmin+1*bb)*p->kmax + k-p->kmin] = 0.0;
             WH[(i-p->imin+2*aa)*p->jmax*p->kmax + (j-p->jmin+2*bb)*p->kmax + k-p->kmin] = 0.0;
             WH[(i-p->imin+3*aa)*p->jmax*p->kmax + (j-p->jmin+3*bb)*p->kmax + k-p->kmin] = 0.0;
            }
        }
}
