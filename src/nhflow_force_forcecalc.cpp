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

#include"nhflow_force.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include <math.h>

void nhflow_force::force_calc(lexer* p, fdm_nhf *d, ghostcell *pgc)
{  
    Ax=0.0;
    Ay=0.0;
    Az=0.0;
    
    Fx=Fy=Fz=0.0;
    A_tot=0.0;
    
    for(n=0;n<polygon_num;++n)
    { 
            // triangle
            if(numpt[n]==3)
            {
            // tri1
            x1 = ccpt[facet[n][0]][0];
            y1 = ccpt[facet[n][0]][1];
            z1 = ccpt[facet[n][0]][2];
            
            x2 = ccpt[facet[n][1]][0];
            y2 = ccpt[facet[n][1]][1];
            z2 = ccpt[facet[n][1]][2];
            
            x3 = ccpt[facet[n][2]][0];
            y3 = ccpt[facet[n][2]][1];
            z3 = ccpt[facet[n][2]][2];
            
            xc = (1.0/3.0)*(x1 + x2 + x3);
            yc = (1.0/3.0)*(y1 + y2 + y3);
            zc = (1.0/3.0)*(z1 + z2 + z3);
            
            at = sqrt(pow(x2-x1,2.0) + pow(y2-y1,2.0) + pow(z2-z1,2.0));
            bt = sqrt(pow(x2-x3,2.0) + pow(y2-y3,2.0) + pow(z2-z3,2.0));
            ct = sqrt(pow(x3-x1,2.0) + pow(y3-y1,2.0) + pow(z3-z1,2.0));
            
            st = 0.5*(at+bt+ct);
            
            A = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
            }
            
            //quadrilidral
            if(numpt[n]==4)
            {
            x1 = ccpt[facet[n][0]][0];
            y1 = ccpt[facet[n][0]][1];
            z1 = ccpt[facet[n][0]][2];
            
            x2 = ccpt[facet[n][1]][0];
            y2 = ccpt[facet[n][1]][1];
            z2 = ccpt[facet[n][1]][2];
            
            x3 = ccpt[facet[n][3]][0];
            y3 = ccpt[facet[n][3]][1];
            z3 = ccpt[facet[n][3]][2];
            
            x4 = ccpt[facet[n][2]][0];
            y4 = ccpt[facet[n][2]][1];
            z4 = ccpt[facet[n][2]][2];
            
            xc = (1.0/4.0)*(x1 + x2 + x3 + x4);
            yc = (1.0/4.0)*(y1 + y2 + y3 + y4);
            zc = (1.0/4.0)*(z1 + z2 + z3 + z4);
            
            //tri1
            at = sqrt(pow(x2-x1,2.0) + pow(y2-y1,2.0) + pow(z2-z1,2.0));
            bt = sqrt(pow(x2-x3,2.0) + pow(y2-y3,2.0) + pow(z2-z3,2.0));
            ct = sqrt(pow(x3-x1,2.0) + pow(y3-y1,2.0) + pow(z3-z1,2.0));
            
            st = 0.5*(at+bt+ct);

            A = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
            
            //tri2
            at = sqrt(pow(x3-x1,2.0) + pow(y3-y1,2.0) + pow(z3-z1,2.0));
            bt = sqrt(pow(x4-x3,2.0) + pow(y4-y3,2.0) + pow(z4-z3,2.0));
            ct = sqrt(pow(x4-x1,2.0) + pow(y4-y1,2.0) + pow(z4-z1,2.0));
            
            st = 0.5*(at+bt+ct);

            A += sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
            }
            
            xp1 = x2-x1;
            yp1 = y2-y1;
            zp1 = z2-z1;
            
            xp2 = x3-x1;
            yp2 = y3-y1;
            zp2 = z3-z1;
            
            nx=fabs(yp1*zp2 - zp1*yp2);
            ny=fabs(zp1*xp2 - xp1*zp2);
            nz=fabs(xp1*yp2 - yp1*xp2);
            
            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            nx/=norm>1.0e-20?norm:1.0e20;
            ny/=norm>1.0e-20?norm:1.0e20;
            nz/=norm>1.0e-20?norm:1.0e20;
            
            //----
            
            i = p->posc_i(xc);
            j = p->posc_j(yc);
            k = p->posc_k(zc);
            
            sgnx = (d->SOLID[Ip1JK] - d->SOLID[Im1JK])/(p->DXP[IM1] + p->DXP[IP]);
            sgny = (d->SOLID[IJp1K] - d->SOLID[IJm1K])/(p->DYP[JM1] + p->DYP[JP]);
            sgnz = (d->SOLID[IJKp1] - d->SOLID[IJKm1])/(p->DZP[KM1] + p->DZP[KP])*p->sigz[IJ];
            
            nx = nx*sgnx/fabs(fabs(sgnx)>1.0e-20?sgnx:1.0e20);
            ny = ny*sgny/fabs(fabs(sgny)>1.0e-20?sgny:1.0e20);
            nz = nz*sgnz/fabs(fabs(sgnz)>1.0e-20?sgnz:1.0e20);
            
            xloc = xc - nx*p->DXP[IP]*p->P91;
            yloc = yc - ny*p->DYP[JP]*p->P91;
            zloc = zc - nz*p->DZP[KP]*d->WL(i,j)*p->P91;
            
            xlocvel = xc + nx*p->DXP[IP];
            ylocvel = yc + ny*p->DYP[JP];
            zlocvel = zc + nz*p->DZP[KP]*d->WL(i,j);
            
            /*uval = p->ccipol4V(d->U, d->WL, d->bed,xlocvel,ylocvel,zlocvel);
            vval = p->ccipol4V(d->V, d->WL, d->bed,xlocvel,ylocvel,zlocvel);
            wval = p->ccipol4V(d->W, d->WL, d->bed,xlocvel,ylocvel,zlocvel);
            
            du = uval/p->DXN[IP];
            dv = vval/p->DYN[JP];
            dw = wval/(p->DZN[KP]*d->WL(i,j));*/
            
            density =   p->ccipol4V(d->RO, d->WL, d->bed,xloc,yloc,zloc);
            viscosity = p->ccipol4V(d->VISC, d->WL, d->bed,xloc,yloc,zloc);
            if(p->P82==1)
            viscosity += p->ccipol4V(d->EV, d->WL, d->bed,xloc,yloc,zloc);
            
            // pressure
            pval   = p->ccipol4V(d->P, d->WL, d->bed,xloc,yloc,zloc);// - p->pressgage;
            etaval = p->ccslipol4(d->eta,xloc,yloc);  
            hspval = (p->wd + etaval - zloc)*p->W1*fabs(p->W22);
            
            /*pval   = p->ccipol4V(d->P, d->WL, d->bed,xc,yc,zc);// - p->pressgage;
            etaval = p->ccslipol4(d->eta,xc,yc);    
            hspval = (p->wd + etaval - zc)*p->W1*fabs(p->W22);*/
            
            // Force
            Fx += -(pval + hspval)*A*nx;
                       //+ 0.0*density*viscosity*A*(du*ny+du*nz);
                       
            Fy += -(pval + hspval)*A*ny;
                      // + 0.0*density*viscosity*A*(dv*nx+dv*nz);
                    
            Fz += -(pval + hspval)*A*nz;
                      // + 0.0*density*viscosity*A*(dw*nx+dw*ny); 
                      
    Ax+=A*nx;    
    Ay+=A*ny;
    Az+=A*nz;
    
    A_tot+=A;
    }
    
    Fx = pgc->globalsum(Fx);
    Fy = pgc->globalsum(Fy);
    Fz = pgc->globalsum(Fz);
    
    Ax = pgc->globalsum(Ax);
    Ay = pgc->globalsum(Ay);
    Az = pgc->globalsum(Az);
    A_tot = pgc->globalsum(A_tot);
}





