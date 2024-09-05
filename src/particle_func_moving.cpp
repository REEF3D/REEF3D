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

#include "particle_func.h"

#include "lexer.h"
#include "fdm.h"

#define PARTICLELOOP for(size_t n=0;n<PP->loopindex;n++) if(PP->Flag[n]>INT32_MIN)

/// @brief Applies advection to positions of particles in @param PP
/// @param p partition object
/// @param a fdm object contains flow field
/// @param PP tracers_obj contains tracer information
/// @param minflag PP.Flag[n] needs to be bigger than minflag for PP[n] to be affected by avection
void particle_func::advect(lexer* p, fdm* a, particles_obj* PP, int minflag, double source_u, double source_v, double source_w)
{
    double coord1, coord2, coord3, u1, u2, v1, v2, w1, w2;
    PARTICLELOOP
        if(PP->Flag[n]>minflag)
        {
            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);
            cellSum[IJK]-=PP->ParcelFactor[n];

            source_w -=settling_vel(p,a,PP,n);
            u1=p->dt*(p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n])+source_u);
            coord1=PP->X[n]+u1;
            
            v1=p->dt*(p->ccipol2(a->v,PP->X[n],PP->Y[n],PP->Z[n])+source_v);
            coord2=PP->Y[n]+v1;
            
            w1=p->dt*(p->ccipol3(a->w,PP->X[n],PP->Y[n],PP->Z[n])+source_w);
            coord3=PP->Z[n]+w1;
            
            
            u2=0.25*u1 + 0.25*p->dt*(p->ccipol1(a->u,coord1,coord2,coord3)+source_u);
            coord1=PP->X[n]+u2;
            
            v2=0.25*v1 + 0.25*p->dt*(p->ccipol2(a->v,coord1,coord2,coord3)+source_v);
            coord2=PP->Y[n]+v2;
            
            w2=0.25*w1 + 0.25*p->dt*(p->ccipol3(a->w,coord1,coord2,coord3)+source_w);
            coord3=PP->Z[n]+w2;
            
            
            PP->X[n] += (2.0/3.0)*u2 + (2.0/3.0)*p->dt*(p->ccipol1(a->u,coord1,coord2,coord3)+source_u);

            PP->Y[n] += (2.0/3.0)*v2 + (2.0/3.0)*p->dt*(p->ccipol2(a->v,coord1,coord2,coord3)+source_v);
            
            PP->Z[n] += (2.0/3.0)*w2 + (2.0/3.0)*p->dt*(p->ccipol3(a->w,coord1,coord2,coord3)+source_w);

            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);
            cellSum[IJK]+=PP->ParcelFactor[n];
        }
}
