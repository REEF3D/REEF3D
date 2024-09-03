/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "partres.h"
#include "particles_obj.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"

void partres::timestep(lexer *p, ghostcell &pgc, particles_obj &PP)
{
        double maxVelU=0,maxVelV=0,maxVelW=0;
        for(size_t n=0;n<PP.loopindex;n++)
        {
            if(PP.Flag[n]>0)
            {
                maxVelU=max(maxVelU,fabs(PP.U[n]));
                maxVelV=max(maxVelV,fabs(PP.V[n]));
                maxVelW=max(maxVelW,fabs(PP.W[n]));
            }
        }

        dx = pgc.globalmin(dx);
        if(p->mpirank==0)
        cout<<"TimeStep: dx: "<<dx<<" vel: "<<sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW)<<endl;
        double meanVel=sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW);
        if(meanVel==0)
        meanVel=dx*p->S14/p->S13;
        p->dtsed = p->S14 * (dx/meanVel);
        p->dtsed = min(p->dtsed,p->S13);
        p->dtsed = min(p->dtsed,p->dt);
        p->dtsed = pgc.globalmin(p->dtsed);
}