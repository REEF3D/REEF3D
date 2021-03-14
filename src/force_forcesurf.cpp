/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include "force.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"


void force::force_surface(lexer* p, fdm *a, ghostcell *pgc)
{
    double ux,vy,wz,vel,pressure,density,viscosity;
    double du,dv,dw;
    double xloc,yloc,zloc;
	double xlocvel,ylocvel,zlocvel;


    FDs=FLs=FTs=Fvert=0.0;
	FD_Press=FD_Shear=0.0;

    for(n=0;n<fnum;++n)
    {
        i=fid[n][0];
        j=fid[n][1];
        k=fid[n][2];

        if(a->phi(i,j,k)>fsfmargin*p->DXM || p->P92==1)
        {
        xloc = fx[n] + fn[n][0]*p->DXP[IP]*p->P91;
        yloc = fy[n] + fn[n][1]*p->DYP[JP]*p->P91;
        zloc = fz[n] + fn[n][2]*p->DZP[KP]*p->P91;
		
		xlocvel = fx[n] + fn[n][0]*p->DXP[IP];
        ylocvel = fy[n] + fn[n][1]*p->DYP[JP];
        zlocvel = fz[n] + fn[n][2]*p->DZP[KP];

        ux = p->ccipol1(a->u,xlocvel,ylocvel,zlocvel);
        vy = p->ccipol2(a->v,xlocvel,ylocvel,zlocvel);
        wz = p->ccipol3(a->w,xlocvel,ylocvel,zlocvel);
		
		pressure = p->ccipol4(a->press,xloc,yloc,zloc);
        
        //cout<<"PRESSURE: "<<pressure<<" x: "<<xloc<<" y: "<<yloc<<" z: "<<zloc<<endl;
		density = p->ccipol4(a->ro,xloc,yloc,zloc);
		viscosity = p->ccipol4(a->visc,xloc,yloc,zloc);

        du = ux/p->DXN[IP];
        dv = vy/p->DYN[JP];
        dw = wz/p->DZN[KP];

        vel = sqrt(ux*ux + vy*vy + wz*wz);

        FTs += -pressure*farea[n]*(fn[n][0] + fn[n][1])
               + density*viscosity*farea[n]*(du*fn[n][1]+dv*fn[n][0]);

        FDs += -(pressure)*farea[n] * fn[n][0]
               + density*viscosity*farea[n]*(du*fn[n][1]);
			
		FD_Press += -(pressure)*farea[n] * fn[n][0];

		FD_Shear += density*viscosity*farea[n]*(du*fn[n][1]);

        FLs += -(pressure)*farea[n] * fn[n][1]
               + density*viscosity*farea[n]*(dv*fn[n][0]);
			   
		Fvert += -(pressure)*farea[n] * fn[n][2]
               + density*viscosity*farea[n]*(du*fn[n][1]+dv*fn[n][0]);
        }

    }
	
    FDs  = pgc->globalsum(FDs);
    FLs  = pgc->globalsum(FLs);
    FTs  = pgc->globalsum(FTs);
	Fvert  = pgc->globalsum(Fvert);
	
	FD_Press  = pgc->globalsum(FD_Press);
	FD_Shear  = pgc->globalsum(FD_Shear);
	
	if(p->count==2)
	Fvert0 = Fvert;

    FDs_norm = FDs/(p->W1*fabs(p->W22)*p->wH*p->wH);
	Fvert_norm =(Fvert-Fvert0)/(p->W1*fabs(p->W22)*PI*p->DXM*0.25);
}


