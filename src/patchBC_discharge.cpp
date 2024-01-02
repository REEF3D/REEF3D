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

#include"patchBC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC::patchBC_discharge(lexer *p, fdm* a, ghostcell *pgc)
{
    double area=0.0;
    double hval=0.0;
    int hcount=0;
    double Ai=0.0;
    double Qi=0.0;
    double Ui=0.0;
    double Hi=0.0;
    double zval=0.0;
    
    // hydrograph interpolation
    // discharge
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->hydroQ_flag==1)
    {
    patch[qq]->Q = patchBC_hydrograph_Q_ipol(p,pgc,qq,patch[qq]->ID);
    }

    
    // Q calc
    for(qq=0;qq<obj_count;++qq)
    {
        Ai=area=hval=Qi=Hi=Ui=zval=0.0;
        hcount=0;
        
        for(n=0;n<patch[qq]->gcb_count;++n)
        {
        i=patch[qq]->gcb[n][0];
        j=patch[qq]->gcb[n][1];
        k=patch[qq]->gcb[n][2];
            
            
            // sides 1
            if(patch[qq]->gcb[n][3]==1)
            {            
            if(a->phi(i,j,k)>=0.5*p->DZN[KP])
            area=p->DYN[JP]*p->DZN[KP];

            if(a->phi(i,j,k)<0.5*p->DZN[KP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 + a->phi(i,j,k));
			
            if(a->phi(i,j,k)>=-0.5*p->DZN[KP] -1.0e-20 && a->phi(i,j,k)<=0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 - fabs(a->phi(i,j,k)));
            
            Ai+=area;
            Qi+=area*(patch[qq]->cosalpha*a->v(i-1,j,k) + patch[qq]->sinalpha*a->u(i-1,j,k));
            
                if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                {
                zval+=-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z();
                ++hcount;
                }
            }
            
            // side 4
            if(patch[qq]->gcb[n][3]==4)
            {
            if(a->phi(i,j,k)>=0.5*p->DZN[KP])
            area=p->DYN[JP]*p->DZN[KP];

            if(a->phi(i,j,k)<0.5*p->DZN[KP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 + a->phi(i,j,k));
			
			if(a->phi(i,j,k)>=-0.5*p->DZN[KP] -1.0e-20 && a->phi(i,j,k)<=0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 - fabs(a->phi(i,j,k)));
            
            Ai+=area;
            Qi+=area*(patch[qq]->cosalpha*a->v(i+1,j,k) + patch[qq]->sinalpha*a->u(i,j,k));
            
                if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                {
                zval+=-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z();
                ++hcount;
                }
            }
            
            // side 3 
            if(patch[qq]->gcb[n][3]==3)
            {
            if(a->phi(i,j,k)>=0.5*p->DZN[KP])
            area=p->DYN[JP]*p->DZN[KP];

            if(a->phi(i,j,k)<0.5*p->DZN[KP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 + a->phi(i,j,k));
			
			if(a->phi(i,j,k)>=-0.5*p->DZN[KP] -1.0e-20 && a->phi(i,j,k)<=0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 - fabs(a->phi(i,j,k)));
            
            Ai+=area;
            Qi+=area*(patch[qq]->sinalpha*a->v(i,j-1,k) + patch[qq]->cosalpha*a->u(i,j-1,k));
            
                if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                {
                zval+=-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z();
                ++hcount;
                }
            }
            
            // side 2
            if(patch[qq]->gcb[n][3]==2)
            {
            if(a->phi(i,j,k)>=0.5*p->DZN[KP])
            area=p->DYN[JP]*p->DZN[KP];

            if(a->phi(i,j,k)<0.5*p->DZN[KP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 + a->phi(i,j,k));
			
			if(a->phi(i,j,k)>=-0.5*p->DZN[KP] -1.0e-20 && a->phi(i,j,k)<=0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 - fabs(a->phi(i,j,k)));
            
            Ai+=area;
            Qi+=area*(patch[qq]->sinalpha*a->v(i,j,k) + patch[qq]->cosalpha*a->u(i,j+1,k));
            
                if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                {
                zval+=-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z();
                ++hcount;
                }
            }
            
            // side 5 
            if(patch[qq]->gcb[n][3]==5)
            {
            if(a->phi(i,j,k)>=0.5*p->DXN[IP])
            area=p->DYN[JP]*p->DXN[IP];

            if(a->phi(i,j,k)<0.5*p->DXN[IP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DXN[IP]*0.5 + a->phi(i,j,k));
			
			if(a->phi(i,j,k)>=-0.5*p->DXN[IP] -1.0e-20 && a->phi(i,j,k)<=0.0)
            area=p->DYN[JP]*(p->DXN[IP]*0.5 - fabs(a->phi(i,j,k)));
            
            Ai+=area;
            Qi+=area*a->w(i,j,k-1);
            }
            
            // side 6
            if(patch[qq]->gcb[n][3]==6)
            {
            if(a->phi(i,j,k)>=0.5*p->DXN[IP])
            area=p->DYN[JP]*p->DXN[IP];

            if(a->phi(i,j,k)<0.5*p->DXN[IP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DXN[IP]*0.5 + a->phi(i,j,k));
			
			if(a->phi(i,j,k)>=-0.5*p->DXN[IP] -1.0e-20 && a->phi(i,j,k)<=0.0)
            area=p->DYN[JP]*(p->DXN[IP]*0.5 - fabs(a->phi(i,j,k)));
            
            Ai+=area;
            Qi+=area*a->w(i,j,k+1);
            }
            
        }
            
            Ai=pgc->globalsum(Ai);
            Qi=pgc->globalsum(Qi);
            zval=pgc->globalsum(zval);
            hcount=pgc->globalisum(hcount);
    
            patch[qq]->Uq = patch[qq]->Q/(Ai>1.0e-20?Ai:1.0e20); 
            Ui = Qi/(Ai>1.0e-20?Ai:1.0e20); 
            
            if(hcount>0)
            {
            Hi=zval/double(hcount);
            }
            
            Hi=pgc->globalmax(Hi);
            

        if(p->mpirank==0)
        {
        cout<<"PatchBC Discharge | ID: "<<patch[qq]->ID<<" Qq: "<<patch[qq]->Q<<" Uq: "<<patch[qq]->Uq<<" Qi: "<<setprecision(5)<<Qi<<" Ui: "<<Ui<<" Hi: "<<Hi<<endl;
        }

        
    }
    
}

void patchBC::patchBC_discharge2D(lexer *p, fdm2D*, ghostcell *pgc, slice &P, slice &Q, slice &eta, slice &bed)
{
    
}