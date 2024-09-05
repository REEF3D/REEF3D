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

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void partres::move_RK2(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb)
{
    double F,G,H;
    particlePerCell(p,pgc,PP);
    particleStressTensor(p,a,pgc,PP);
    timestep(p,pgc,PP);
    
    // RK step 1 
    for(size_t n=0;n<PP.loopindex;n++)
    if(PP.Flag[n]>=0)
    {
        if(p->Q11==1)
        advec_plain(p, a, PP, n, s, pturb, 
                        PP.X, PP.Y, PP.Z, PP.U, PP.V, PP.W,
                        F, G, H, 1.0);
        
        if(p->Q11==2)
        advec_pic(p, a, PP, n, s, pturb, 
                        PP.X, PP.Y, PP.Z, PP.U, PP.V, PP.W,
                        F, G, H, 1.0);
                        
        //cout<<"F: "<<F<<" G: "<<G<<" H: "<<H<<endl;
                            
        // Velocity update
        PP.URK1[n] = PP.U[n] + p->dtsed*F;
        PP.VRK1[n] = PP.V[n] + p->dtsed*G;
        PP.WRK1[n] = 0.0; //PP.W[n] + p->dtsed*H;
        
        // Position update
        PP.XRK1[n] = PP.X[n] + p->dtsed*PP.URK1[n];
        PP.YRK1[n] = PP.Y[n] + p->dtsed*PP.VRK1[n];
        PP.ZRK1[n] = PP.Z[n] + p->dtsed*PP.WRK1[n];

        // Particel sum update
        cellSum[IJK]-=PP.ParcelFactor[n];
        bedChange[IJ]-=PP.ParcelFactor[n];
        i=p->posc_i(PP.XRK1[n]);
        j=p->posc_j(PP.YRK1[n]);
        k=p->posc_k(PP.ZRK1[n]);
        cellSum[IJK]+=PP.ParcelFactor[n];
        bedChange[IJ]+=PP.ParcelFactor[n];
        particleStressTensorUpdateIJK(p,a,PP);
    }
    
    particleStressTensor(p,a,pgc,PP);
    
    
    // RK step 2
    for(size_t n=0;n<PP.loopindex;n++)
    if(PP.Flag[n]>=0)
    {
        if(p->Q11==1)
        advec_plain(p, a, PP, n, s, pturb, 
                        PP.XRK1, PP.YRK1, PP.ZRK1, PP.URK1, PP.VRK1, PP.WRK1,
                        F, G, H, 0.5);
                        
        if(p->Q11==2)
        advec_pic(p, a, PP, n, s, pturb, 
                        PP.XRK1, PP.YRK1, PP.ZRK1, PP.URK1, PP.VRK1, PP.WRK1,
                        F, G, H, 0.5);
                        
        // Velocity update
        PP.U[n] = 0.5*PP.U[n] + 0.5*PP.URK1[n] + 0.5*p->dtsed*F;
        PP.V[n] = 0.5*PP.V[n] + 0.5*PP.VRK1[n] + 0.5*p->dtsed*G;
        PP.W[n] = 0.0;//0.5*PP.W[n] + 0.5*PP.URK1[n] + 0.5*p->dtsed*H;
        
        // Position update
        PP.X[n] = 0.5*PP.X[n] + 0.5*PP.XRK1[n] + 0.5*p->dtsed*PP.U[n];
        PP.Y[n] = 0.5*PP.Y[n] + 0.5*PP.YRK1[n] + 0.5*p->dtsed*PP.V[n];
        PP.Z[n] = 0.5*PP.Z[n] + 0.5*PP.ZRK1[n] + 0.5*p->dtsed*PP.W[n];

        // Particel sum update
        cellSum[IJK]-=PP.ParcelFactor[n];
        bedChange[IJ]-=PP.ParcelFactor[n];
        i=p->posc_i(PP.X[n]);
        j=p->posc_j(PP.Y[n]);
        k=p->posc_k(PP.Z[n]);
        cellSum[IJK]+=PP.ParcelFactor[n];
        bedChange[IJ]+=PP.ParcelFactor[n];
        particleStressTensorUpdateIJK(p,a,PP);
    }

    if(p->mpirank==0)
    {
        p->sedtime += p->dtsed;
        cout<<"Sediment time: "<<p->sedtime<<" time step: "<<p->dtsed<<endl;
    }
    
}