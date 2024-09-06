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
Author: Fabian Knoblauch
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

void VOF_PLIC::redistancePhiByPlane_Bonn
(
    fdm* a,
    lexer* p
)
{
    int bandWidth=1;
    int jjbandWidth;
    if(p->j_dir>0)
        jjbandWidth=1;
    else
        jjbandWidth=0;
    double phitemp;
    double dd=-alpha(i,j,k);
    int ip, jp, kp;
    
    for(int ii=-bandWidth;ii<bandWidth+1;ii++)
    {
        for(int jj=-jjbandWidth;jj<jjbandWidth+1;jj++)
        {
            for(int kk=-bandWidth;kk<bandWidth+1;kk++)
            {
                        ip = i + ii;
						jp = j + jj;
						kp = k + kk;
  
						if (ip >= 0 && jp >= 0 && kp >= 0 && ip < p->knox && jp < p->knoy && kp < p->knoz)
						{
                                if(ii==0 && jj==0 && kk==0)
                                {
                                    if(a->vof(i,j,k)>0.0001 && a->vof(i,j,k)<0.9999)
                                    {
                                        if(alpha(i,j,k)*phistep(i,j,k)>=0.0)
                                        {
                                        phitemp=alpha(i,j,k);
                                        if(fabs(phitemp)<fabs(phiaux(i,j,k)))
                                            phiaux(i,j,k)=phitemp;
                                        }
                                    }
                                }
                                else
                                {   
                                    phitemp=ShortestDistanceOnBoundaryCandidate(a,p,ii,jj,kk,dd);
                                    if(fabs(phitemp)<1E05)
                                    {
                                        if(fabs(phitemp)<fabs(phiaux(ip,jp,kp)))
                                            phiaux(ip,jp,kp)=copysign(phitemp,phistep(i,j,k));
                                    }
                                    else
                                    {   
                                        phitemp=ProjectionPointCandidate(a,p,ii,jj,kk,dd);
                                        if(fabs(phitemp)<1E05)
                                        {
                                            if(fabs(phitemp)<fabs(phiaux(ip,jp,kp)))
                                                phiaux(ip,jp,kp)=copysign(phitemp,phistep(i,j,k));
                                        }   
                                        else
                                        {
                                            phitemp=IntersectionPointCandidate(a,p,ii,jj,kk,dd);
                                            if(fabs(phitemp)<1E05)
                                            {
                                                if(fabs(phitemp)<fabs(phiaux(ip,jp,kp)))
                                                    phiaux(ip,jp,kp)=copysign(phitemp,phistep(i,j,k));
                                            }
                                            else
                                            {
                                                cout<<"no nearest phase pint found"<<endl;
                                                phitemp=phistep(ip,jp,kp);
                                                if(fabs(phitemp)<fabs(phiaux(ip,jp,kp)))
                                                    phiaux(ip,jp,kp)=copysign(phitemp,phistep(i,j,k));
                                            }
                                        }
                                    }
                                }
                        }       
            }          
        }
    }
}
    


double VOF_PLIC::ShortestDistanceOnBoundaryCandidate
(
    fdm* a,
    lexer* p,
    int ii,
    int jj,
    int kk,
    double dd
)
{
    // NOT ADAPTED FOR STRETCHED GRID YET !!!
    
    double xvx, xvy,xvz,dist,ret;
    xvx=0.5*p->DXN[IP]*max(-1,min(1,ii));
    xvy=0.5*p->DYN[JP]*max(-1,min(1,jj));
    xvz=0.5*p->DZN[KP]*max(-1,min(1,kk));
    dist=nx(i,j,k)*xvx+ny(i,j,k)*xvy+nz(i,j,k)*xvz+dd;
    if(-dist*(phistep(i+ii,j+jj,k+kk))<0.0)
        ret=copysign(sqrt(xvx*xvx+xvy*xvy+xvz*xvz),(vofstep(i+ii,j+jj,k+kk)-0.5));
    else 
        ret=1E06;
        
    return ret;
        
}

double VOF_PLIC::ProjectionPointCandidate
(
    fdm* a,
    lexer* p,
    int ii,
    int jj,
    int kk,
    double dd
)
{
    // NOT ADAPTED FOR STRETCHED GRID YET !!!
    
    double xpx, xpy,xpz, dist, ret;
    dist=nx(i,j,k)*(p->DXN[IP]*ii)+ny(i,j,k)*(p->DYN[JP]*jj)+nz(i,j,k)*(p->DZN[KP]*kk)+dd;
    xpx=(p->DXN[IP]*ii)-dist*nx(i,j,k);
    xpy=(p->DYN[JP]*jj)-dist*ny(i,j,k);
    xpz=(p->DZN[KP]*kk)-dist*nz(i,j,k);
    if(fabs(xpx)<=p->DXN[IP]*0.5+1E-04*p->DXN[IP] && fabs(xpy)<=0.5*p->DYN[JP]+1E-04*p->DYN[JP] && fabs(xpz)<=0.5*p->DZN[KP]+1E-04*p->DZN[KP])
    {
        double dx, dy, dz;
        dx=(p->DXN[IP]*ii)-xpx;
        dy=(p->DYN[JP]*jj)-xpy;
        dz=(p->DZN[KP]*kk)-xpz;
        ret=copysign(sqrt(dx*dx+dy*dy+dz*dz),(vofstep(i+ii,j+jj,k+kk)-0.5));
    }
    else
        ret=1E06;
        
    return ret;
}

double VOF_PLIC::IntersectionPointCandidate
(
    fdm* a,
    lexer* p,
    int ii,
    int jj,
    int kk,
    double dd
)// NOT ADAPTED FOR STRETCHED GRID YET !!!
{
    double dxn, dyn, dzn, x_x,x_y,x_z;
    double p1x,p1y,p1z, p2x,p2y,p2z,Dp1,Dp2;
    double pisx,pisy,pisz;
    double ppr1x,ppr1y,ppr1z,ppr2x,ppr2y,ppr2z;
    double distl1, distl2;
    double r0, n_x,n_y,n_z;
    double nxl1,nyl1,nzl1,r0l1;
    double nxl2,nyl2,nzl2,r0l2;
    double mindist=1E06;
    int dimchange;
    dxn=p->DXN[IP];
    dyn=p->DYN[JP];
    dzn=p->DZN[KP];
    x_x=dxn*ii;
    x_y=dyn*jj;
    x_z=dzn*kk;
    r0=alpha(i,j,k);
    n_x=nx(i,j,k);
    n_y=ny(i,j,k);
    n_z=nz(i,j,k);
    
    //interate through all 12 lines of a cube
    for(int cubeline=1; cubeline<13; cubeline++)
    {
        if(cubeline<=4)
        {
            p1x=-0.5*dxn;
            p2x=-0.5*dxn;
            
        }
        else if(cubeline <=8)
        {
            p1x=-0.5*dxn;
            p2x=0.5*dxn;
            dimchange=1;
        }
        else
        {
            p1x=0.5*dxn;
            p2x=0.5*dxn;
        }
        
        if(cubeline==1 || cubeline==3 || cubeline==9 || cubeline==11)
        {
            p1y=-0.5*dyn;
            p2y=0.5*dyn;
            dimchange=2;
        }
        else if(cubeline==2 || cubeline==5 || cubeline==6 || cubeline==10)
        {
            p1y=-0.5*dyn;
            p2y=-0.5*dyn;
        }
        else
        {
            p1y=0.5*dyn;
            p2y=0.5*dyn;
        }
        
        if(cubeline==2 || cubeline==4 || cubeline==10 || cubeline==12)
        {
            p1z=-0.5*dzn;
            p2z=0.5*dzn;
            dimchange=3;
        }
        else if(cubeline==1 || cubeline==5 || cubeline==8 || cubeline==9)
        {
            p1z=-0.5*dzn;
            p2z=-0.5*dzn;
        }
        else
        {
            p1z=0.5*dzn;
            p2z=0.5*dzn;
        }
        
        Dp1=n_x*p1x+n_y*p1y+n_z*p1z-r0;
        Dp2=n_x*p2x+n_y*p2y+n_z*p2z-r0;
        
        if(Dp1*Dp2<0.0)
        {
            double dx,dy,dz;
            if(dimchange==1)
            {   
                //würfelkante entlang x
                pisy=p1y;
                pisz=p1z;
                pisx=(r0-n_y*pisy-n_z*pisz)/n_x;
            }
            else if (dimchange==2)
            {
                //würfelkante entlang y
                pisx=p1x;
                pisz=p1z;
                pisy=(r0-n_x*pisx-n_z*pisz)/n_y;
            }
            else
            {
                //würfelkante entlang z
                pisx=p1x;
                pisy=p1y;
                pisz=(r0-n_x*pisx-n_y*pisy)/n_z;
            }
            dx=x_x-pisx;
            dy=x_y-pisy;
            dz=x_z-pisz;
            mindist=min(sqrt(dx*dx+dy*dy+dz*dz),mindist);
            
            if(dimchange==1)
            {   // x=change
            
                //1. y festhalten
                nxl1=n_x/sqrt(n_x*n_x+n_z*n_z);
                nzl1=n_z/sqrt(n_x*n_x+n_z*n_z);
                r0l1=(r0-n_y*pisy)/sqrt(n_x*n_x+n_z*n_z);
                distl1=nxl1*x_x+nzl1*x_z-r0l1;
                ppr1x=x_x-distl1*nxl1;
                ppr1y=pisy;
                ppr1z=x_z-distl1*nzl1;
                
                //2, z festhalten
                nxl2=n_x/sqrt(n_x*n_x+n_y*n_y);
                nyl2=n_y/sqrt(n_x*n_x+n_y*n_y);
                r0l2=(r0-n_z*pisz)/sqrt(n_x*n_x+n_z*n_z);
                distl2=nxl2*x_x+nyl2*x_y-r0l2;
                ppr2x=x_x-distl2*nxl2;
                ppr2y=x_y-distl2*nyl2;
                ppr2z=pisz;
            }
            else if(dimchange==2)
            {
                // y=change
                
                // 1. x festhalten
                nyl1=n_y/sqrt(n_y*n_y+n_z*n_z);
                nzl1=n_z/sqrt(n_y*n_y+n_z*n_z);
                r0l1=(r0-n_x*pisx)/sqrt(n_y*n_y+n_z*n_z);
                distl1=nyl1*x_y+nzl1*x_z-r0l1;
                ppr1x=pisx;
                ppr1y=x_y-distl1*nyl1;
                ppr1z=x_z-distl1*nzl1;
                
                // 2. z festhalten
                nxl2=n_x/sqrt(n_x*n_x+n_y*n_y);
                nyl2=n_y/sqrt(n_x*n_x+n_y*n_y);
                r0l2=(r0-n_z*pisz)/sqrt(n_x*n_x+n_y*n_y);
                distl2=nxl2*x_x+nyl2*x_y-r0l2;
                ppr2x=x_x-distl2*nxl2;
                ppr2y=x_y-distl2*nyl2;
                ppr2z=pisz;
            }
            else
            {
                //z=change
                
                //1. x festhalten
                nyl1=n_y/sqrt(n_y*n_y+n_z*n_z);
                nzl1=n_z/sqrt(n_y*n_y+n_z*n_z);
                r0l1=(r0-n_x*pisx)/sqrt(n_y*n_y+n_z*n_z);
                distl1=nyl1*x_y+nzl1*x_z-r0l1;
                ppr1x=pisx;
                ppr1y=x_y-distl1*nyl1;
                ppr1z=x_z-distl1*nzl1;
                
                //2. y festhalten
                nxl2=n_x/sqrt(n_x*n_x+n_z*n_z);
                nzl2=n_z/sqrt(n_x*n_x+n_z*n_z);
                r0l2=(r0-n_y*pisy)/sqrt(n_x*n_x+n_z*n_z);
                distl2=nxl2*x_x+nzl2*x_z-r0l2;
                ppr2x=x_x-distl2*nxl2;
                ppr2y=pisy;
                ppr2z=x_z-distl2*nzl2;
            }
            
            if(fabs(ppr1x)<=0.5*dxn+1E-04*dxn && fabs(ppr1y)<=0.5*dyn+1E-04*dyn && fabs(ppr1z)<=0.5*dzn+1E-04*dzn)
            {
                dx=x_x-ppr1x;
                dy=x_y-ppr1y;
                dz=x_z-ppr1z;
                mindist=min(sqrt(dx*dx+dy*dy+dz*dz),mindist);
            }
            
            if(fabs(ppr2x)<=0.5*dxn+1E-04*dxn && fabs(ppr2y)<=0.5*dyn+1E-04*dyn && fabs(ppr2z)<=0.5*dzn+1E-04*dzn)
            {
                dx=x_x-ppr2x;
                dy=x_y-ppr2y;
                dz=x_z-ppr2z;
                mindist=min(sqrt(dx*dx+dy*dy+dz*dz),mindist);
            }
            
            
        }
    }
    
    return mindist;
}

