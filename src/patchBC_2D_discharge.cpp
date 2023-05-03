/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"patchBC_2D.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC_2D::patchBC_discharge2D(lexer *p, fdm2D* b, ghostcell *pgc, slice &P, slice &Q, slice &eta, slice &bed)
{
    double area=0.0;
    double hval=0.0;
    int hcount=0;
    double Ai=0.0;
    double Qi=0.0;
    double Ui=0.0;
    double Hi=0.0;
    
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
        Ai=area=hval=Qi=Hi=Ui=0.0;
        hcount=0;
        
        for(n=0;n<patch[qq]->gcb_count;++n)
        {
        i=patch[qq]->gcb[n][0];
        j=patch[qq]->gcb[n][1];
            
            
            // sides 1 & 4
            if(patch[qq]->gcb[n][3]==1 && p->wet[IJ]==1)
            {
            area = p->DYN[JP]*b->hp(i-1,j);
            
            Ai+=area;
            Qi+=area*(patch[qq]->cosalpha*b->Q(i-1,j) + patch[qq]->sinalpha*b->P(i-1,j));
            
            hval += b->hp(i-1,j);
            ++hcount;
            }
            
            if(patch[qq]->gcb[n][3]==4 && p->wet[IJ]==1)
            {
            area = p->DYN[JP]*b->hp(i+1,j);
            
            Ai+=area;
            Qi+=area*(patch[qq]->cosalpha*b->Q(i+1,j) + patch[qq]->sinalpha*b->P(i+1,j));
            
            hval += b->hp(i+1,j);
            ++hcount;
            }
            
            // sides 3 & 2
            if(patch[qq]->gcb[n][3]==3 && p->wet[IJ]==1)
            {
            area = p->DXN[IP]*b->hp(i,j-1);
            
            Ai+=area;
            Qi+=area*(patch[qq]->sinalpha*b->Q(i,j-1) + patch[qq]->cosalpha*b->P(i,j+1));
            
            hval += b->hp(i,j-1);
            ++hcount;
            }
            
            if(patch[qq]->gcb[n][3]==2  && p->wet[IJ]==1)
            {
            area = p->DXN[IP]*b->hp(i,j+1);
            
            Ai+=area;
            Qi+=area*(patch[qq]->sinalpha*b->Q(i,j+1) + patch[qq]->cosalpha*b->P(i,j+1));
            hval += b->hp(i,j+1);
            ++hcount;
            }
            
        }
            
            Ai=pgc->globalsum(Ai);
            Qi=pgc->globalsum(Qi);
            Hi=pgc->globalsum(hval);
            hcount=pgc->globalisum(hcount);
            
            Hi = Hi/(hcount>1.0e-20?hcount:1.0e20); 
    
            patch[qq]->Uq = patch[qq]->Q/(Ai>1.0e-20?Ai:1.0e20); 
            Ui = Qi/(Ai>1.0e-20?Ai:1.0e20); 
            

        if(p->mpirank==0)
        {
        cout<<"PatchBC Discharge | ID: "<<patch[qq]->ID<<" Qq: "<<patch[qq]->Q<<" Uq: "<<patch[qq]->Uq<<" Qi: "<<setprecision(5)<<Qi<<" Ui: "<<Ui<<" Ai: "<<Ai<<" Hi: "<<Hi<<endl;
        }

    }

}


void patchBC_2D::patchBC_discharge(lexer *p, fdm* a, ghostcell *pgc)
{
    
}