/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"

// particle boundary conditions
// walls: inelastic reflection (no penetration)
// open boundaries (inflow, outflow): parcels leaving the domain are removed, their volume is logged
// sides: 1 -x, 2 +y, 3 -y, 4 +x, 5 -z, 6 +z
void CPM::boundcheck(lexer *p, int mode)
{
    double *PX,*PY,*PZ,*PU,*PV,*PW;
    
    if(mode==1)
    {
        PX=P.XRK1; PY=P.YRK1; PZ=P.ZRK1;
        PU=P.URK1; PV=P.VRK1; PW=P.WRK1;
    }
    else
    {
        PX=P.X; PY=P.Y; PZ=P.Z;
        PU=P.U; PV=P.V; PW=P.W;
    }
    
    double eps = 1.0e-3*hmin;
    const double vpar = P.ParcelFactor*Vp;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        // open boundaries
        if((PX[n]<p->global_xmin && open_side[0]==1) || (PX[n]>p->global_xmax && open_side[3]==1)
        || (p->j_dir==1 && ((PY[n]<p->global_ymin && open_side[2]==1) || (PY[n]>p->global_ymax && open_side[1]==1))))
        {
            P.remove(n);
            outvol += vpar;
            continue;
        }
        
        // walls in x
        if(perx==0)
        {
        if(PX[n]<p->global_xmin)
        {
            PX[n] = MIN(2.0*p->global_xmin - PX[n], p->global_xmin + 0.5*hmin);
            PU[n] = MAX(PU[n],0.0);
        }
        
        if(PX[n]>p->global_xmax)
        {
            PX[n] = MAX(2.0*p->global_xmax - PX[n], p->global_xmax - 0.5*hmin);
            PU[n] = MIN(PU[n],0.0);
        }
        
        PX[n] = MAX(PX[n],p->global_xmin+eps);
        PX[n] = MIN(PX[n],p->global_xmax-eps);
        }
        
        // walls in y
        if(p->j_dir==1 && pery==0)
        {
            if(PY[n]<p->global_ymin)
            {
                PY[n] = MIN(2.0*p->global_ymin - PY[n], p->global_ymin + 0.5*hmin);
                PV[n] = MAX(PV[n],0.0);
            }
            
            if(PY[n]>p->global_ymax)
            {
                PY[n] = MAX(2.0*p->global_ymax - PY[n], p->global_ymax - 0.5*hmin);
                PV[n] = MIN(PV[n],0.0);
            }
            
            PY[n] = MAX(PY[n],p->global_ymin+eps);
            PY[n] = MIN(PY[n],p->global_ymax-eps);
        }
        
        // walls in z
        if(PZ[n]<p->global_zmin)
        {
            PZ[n] = MIN(2.0*p->global_zmin - PZ[n], p->global_zmin + 0.5*hmin);
            PW[n] = MAX(PW[n],0.0);
        }
        
        if(PZ[n]>p->global_zmax)
        {
            PZ[n] = MAX(2.0*p->global_zmax - PZ[n], p->global_zmax - 0.5*hmin);
            PW[n] = MIN(PW[n],0.0);
        }
        
        PZ[n] = MAX(PZ[n],p->global_zmin+eps);
        PZ[n] = MIN(PZ[n],p->global_zmax-eps);
    }
}
