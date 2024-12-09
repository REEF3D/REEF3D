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

#include"6DOF_obj.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

void sixdof_obj::hydrodynamic_forces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, bool finalize)
{
	// forcecalc
    if(p->X60==1)
    force_calc_stl(p,d,pgc,finalize);
        
    
    if(p->X60==2)
    {
    triangulation(p,d,pgc);
	reconstruct(p,d);
    force_calc_lsm(p,d,pgc);
        
    deallocate(p,d,pgc);
    }
} 

void sixdof_obj::force_calc_lsm(lexer* p, fdm_nhf *d, ghostcell *pgc)
{  
    Ax=0.0;
    Ay=0.0;
    Az=0.0;
    
    Fx=Fy=Fz=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    

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
            k = p->posc_sig(i,j,zc);
            
            sgnx = (d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IM1] + p->DXP[IP]);
            sgny = (d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JM1] + p->DYP[JP]);
            sgnz = (d->FB[IJKp1] - d->FB[IJKm1])/((p->DZP[KM1] + p->DZP[KP])*d->WL(i,j));
            
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
            pval   = p->ccipol7P(d->P, d->WL, d->bed, xc, yc, zc);// - p->pressgage;
            etaval = p->ccslipol4(d->eta,xc,yc);  
            hspval = (p->wd + etaval - zc)*p->W1*fabs(p->W22);
    
            
            // Force
            Fx = -(pval + hspval)*A*nx;
                       //+ 0.0*density*viscosity*A*(du*ny+du*nz);
                       
            if(p->j_dir==1)
            Fy = -(pval + hspval)*A*ny;
                      // + 0.0*density*viscosity*A*(dv*nx+dv*nz);
                    
            Fz = -(pval + hspval)*A*nz;
                      // + 0.0*density*viscosity*A*(dw*nx+dw*ny); 
                      
                      
            Ax+=A*nx;    
            Ay+=A*ny;
            Az+=A*nz;
    
    
            Xe += Fx;
            Ye += Fy;
            Ze += Fz;

            Ke += (yc - c_(1))*Fz - (zc - c_(2))*Fy;
            Me += (zc - c_(2))*Fx - (xc - c_(0))*Fz;
            Ne += (xc - c_(0))*Fy - (yc - c_(1))*Fx;
    }
    
    Xe = pgc->globalsum(Xe);
	Ye = pgc->globalsum(Ye);
	Ze = pgc->globalsum(Ze);
	Ke = pgc->globalsum(Ke);
	Me = pgc->globalsum(Me);
	Ne = pgc->globalsum(Ne);
    
    Fx = Xe;
    Fy = Ye;
    Fz = Ze;
    
    Ax = pgc->globalsum(Ax);
    Ay = pgc->globalsum(Ay);
    Az = pgc->globalsum(Az);
    
    Xe += p->W20*Mass_fb;
	Ye += p->W21*Mass_fb;
	Ze += p->W22*Mass_fb;
    
    if(p->mpirank==0)
    {
    cout<<"Mass_fb: "<<Mass_fb<<" G_fb: "<<9.81*Mass_fb<<endl;
    cout<<"Fx: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<endl;
    cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
    }
    
    /*
    // Print results	
    if (p->mpirank==0 && finalize==1) 
    {
        printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   
    }*/
}

void sixdof_obj::allocate(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    p->Iarray(tri,numtri,4);
    p->Darray(pt,numvert,3);
    p->Darray(ls,numvert);
    p->Iarray(facet,numtri,4);
    p->Iarray(confac,numtri);
    p->Iarray(numfac,numtri);
	p->Iarray(numpt,numtri);
    p->Darray(ccpt,numtri*4,3);
    
    // ini
    int n,q;
    
    for(n=0; n<numtri; ++n)
    {
    for(q=0;q<4;++q)
    tri[n][q]=0;
    
    for(q=0;q<4;++q)
    facet[n][q]=0;
    
    confac[n]=0;
    numfac[n]=0;
    numpt[n]=0;
    }
    
    for(n=0; n<numvert; ++n)
    {
    for(q=0;q<3;++q)
    pt[n][q]=0.0;
    
    ls[n]=0.0;
    }
    
    for(n=0; n<numtri*4; ++n)
    {
    for(q=0;q<3;++q)
    ccpt[n][q]=0;
    }
}

void sixdof_obj::deallocate(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    p->del_Iarray(tri,numtri,4);
    p->del_Darray(pt,numvert,3);
    p->del_Darray(ls,numvert);
    p->del_Iarray(facet,numtri,4);
    p->del_Iarray(confac,numtri);
    p->del_Iarray(numfac,numtri);
	p->del_Iarray(numpt,numtri);
    p->del_Darray(ccpt,numtri*4,3);
}

