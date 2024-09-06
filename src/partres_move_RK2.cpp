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
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void partres::move_RK2_step1(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb, int &xchanged, int &removed)
{
    particles_obj Send[6]={particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),
    particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1)};
    particles_obj Recv[6]={particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),
    particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1)};
    
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

        addParticleForTransfer(p,PP,n,Send,xchanged);
    }
    

    {
        pgc.para_tracersobj(p,Send,Recv);

        size_t sum=PP.size;
        for(int n=0;n<6;n++)
            sum += Recv[n].size;
        if(sum>PP.capacity)
            PP.reserve(sum);

        for(int n=0;n<6;n++)
        {
            for(size_t m=0;m<Recv[n].loopindex;m++)
            {
                i = p->posc_i(Recv[n].XRK1[m]);
                j = p->posc_j(Recv[n].YRK1[m]);
                k = p->posc_k(Recv[n].ZRK1[m]);
                transfer(p,Recv[n],m);
                cellSum[IJK]+=Recv[n].ParcelFactor[n];
            }
            PP.add_obj(&Recv[n]);
        }
    }


    {
        bool inBounds=false;
        int i,j,k;
        boundarycheck bounderies;

        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]>0)
            {
                i = p->posc_i(PP.XRK1[n]);
                j = p->posc_j(PP.YRK1[n]);
                k = p->posc_k(PP.ZRK1[n]);

                inBounds=bounderies.minboundcheck(p,i,j,k,1);
                if (inBounds)
                    inBounds=bounderies.maxboundcheck(p,i,j,k,1);

                // remove out of bounds particles
                if(!inBounds)
                {
                    remove(p,PP,n);
                    PP.erase(n);
                    removed++;
                }
            }
    }
    
    particleStressTensor(p,a,pgc,PP);
}
    
void partres::move_RK2_step2(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb, int &xchanged, int &removed)
{
    particles_obj Send[6]={particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),
    particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1)};
    particles_obj Recv[6]={particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),
    particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1),particles_obj(maxcount,PP.d50,PP.density,1)};

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

        addParticleForTransfer(p,PP,n,Send,xchanged);
    }

    {
        pgc.para_tracersobj(p,Send,Recv);

        size_t sum=PP.size;
        for(int n=0;n<6;n++)
            sum += Recv[n].size;
        if(sum>PP.capacity)
            PP.reserve(sum);

        for(int n=0;n<6;n++)
        {
            for(size_t m=0;m<Recv[n].loopindex;m++)
            {
                i = p->posc_i(Recv[n].X[m]);
                j = p->posc_j(Recv[n].Y[m]);
                k = p->posc_k(Recv[n].Z[m]);
                transfer(p,Recv[n],m);
                cellSum[IJK]+=Recv[n].ParcelFactor[n];
            }
            PP.add_obj(&Recv[n]);
        }
    }


    {
        bool inBounds=false;
        int removed=0;
        int i,j,k;
        boundarycheck bounderies;

        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]>0)
            {
                i = p->posc_i(PP.X[n]);
                j = p->posc_j(PP.Y[n]);
                k = p->posc_k(PP.Z[n]);

                inBounds=bounderies.minboundcheck(p,i,j,k,1);
                if (inBounds)
                    inBounds=bounderies.maxboundcheck(p,i,j,k,1);

                // remove out of bounds particles
                if(!inBounds)
                {
                    remove(p,PP,n);
                    PP.erase(n);
                    removed++;
                }
            }
    }

    if(p->mpirank==0)
    {
        p->sedtime += p->dtsed;
        cout<<"Sediment time: "<<p->sedtime<<" time step: "<<p->dtsed<<endl;
    }
}

void partres::move_RK2(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb)
{
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