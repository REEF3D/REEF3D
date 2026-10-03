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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include<fstream>
#include<cstdio>
#include<sys/stat.h>

// parcel state for the hotstart, written next to the CFD state files:
// ./REEF3D_CFD_STATE/REEF3D-CFD-CPM-State-<num>-<rank>.r3d
void CPM::state_write(lexer *p, int num)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_CFD_STATE/REEF3D-CFD-CPM-State-%08i-%06i.r3d",num,p->mpirank+1);
    
    ofstream result;
    result.open(name, ios::binary);
    
    if(!result.is_open())
    return;
    
    int numpt=0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    ++numpt;
    
    int version=1;
    result.write((char*)&version, sizeof(int));
    result.write((char*)&numpt, sizeof(int));
    result.write((char*)&P.ParcelFactor, sizeof(double));
    result.write((char*)&outvol, sizeof(double));
    
    double val[8];
    int flag;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        val[0]=P.X[n]; val[1]=P.Y[n]; val[2]=P.Z[n];
        val[3]=P.U[n]; val[4]=P.V[n]; val[5]=P.W[n];
        val[6]=P.D[n]; val[7]=P.RO[n];
        flag=P.Flag[n];
        
        result.write((char*)val, 8*sizeof(double));
        result.write((char*)&flag, sizeof(int));
    }
    
    result.close();
}

// replaces the seeded parcels by the parcels of the state file, if it exists
void CPM::state_read(lexer *p, ghostcell *pgc, int num)
{
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_CFD_STATE/REEF3D-CFD-CPM-State-%08i-%06i.r3d",num,p->mpirank+1);
    
    ifstream result;
    result.open(name, ios::binary);
    
    int ok = result.is_open() ? 1 : 0;
    ok = pgc->globalimin(ok);
    
    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"CPM hotstart: no parcel state file found, the parcels are seeded from the bed"<<endl;
        
        return;
    }
    
    int version,numpt;
    result.read((char*)&version, sizeof(int));
    result.read((char*)&numpt, sizeof(int));
    result.read((char*)&P.ParcelFactor, sizeof(double));
    result.read((char*)&outvol, sizeof(double));
    
    // remove all parcels
    for(n=0;n<P.index;++n)
    P.Flag[n]=EMPTY;
    
    P.resize(p,numpt+100);
    
    P.index_empty=0;
    for(n=P.index-1;n>=0;--n)
    {
        P.Flag[n]=EMPTY;
        P.Empty[P.index_empty]=n;
        ++P.index_empty;
    }
    
    double val[8];
    int flag;
    
    for(int q=0;q<numpt;++q)
    {
        result.read((char*)val, 8*sizeof(double));
        result.read((char*)&flag, sizeof(int));
        
        --P.index_empty;
        n=P.Empty[P.index_empty];
        
        P.X[n]=P.XRK1[n]=val[0];
        P.Y[n]=P.YRK1[n]=val[1];
        P.Z[n]=P.ZRK1[n]=val[2];
        P.U[n]=P.URK1[n]=val[3];
        P.V[n]=P.VRK1[n]=val[4];
        P.W[n]=P.WRK1[n]=val[5];
        P.D[n]=val[6];
        P.RO[n]=val[7];
        P.Flag[n]=flag;
    }
    
    result.close();
    
    restored=1;
    
    int total = pgc->globalisum(numpt);
    
    if(p->mpirank==0)
    cout<<"CPM hotstart: "<<total<<" parcels read from the state files "<<num<<endl;
}

// sediment log: parcels, sediment volume, outflow and the volumetric transport rate per unit width
//   q_x = sum(V_p u_p)/(L_x L_y)
void CPM::sedlog(lexer *p, ghostcell *pgc)
{
    const double vpar = P.ParcelFactor*Vp;
    
    double vol=0.0, qx=0.0, qy=0.0, vmov=0.0;
    int np=0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        ++np;
        vol += vpar;
        qx += vpar*P.U[n];
        qy += vpar*P.V[n];
        
        if(P.U[n]*P.U[n] + P.V[n]*P.V[n] + P.W[n]*P.W[n] > 1.0e-6)
        vmov += vpar;
    }
    
    np = pgc->globalisum(np);
    vol = pgc->globalsum(vol);
    qx = pgc->globalsum(qx);
    qy = pgc->globalsum(qy);
    vmov = pgc->globalsum(vmov);
    double vout = pgc->globalsum(outvol);
    
    double Lx = p->global_xmax-p->global_xmin;
    double Ly = p->j_dir==1 ? p->global_ymax-p->global_ymin : p->DYN[marge];
    
    if(p->mpirank==0)
    {
        if(logini==0)
        {
            mkdir("./REEF3D_CFD_CPM_Particle",0777);
            logout.open("./REEF3D_CFD_CPM_Particle/REEF3D-CFD-CPM-Log.dat");
            logout<<"time\tparcels\tsediment_volume[m3]\toutflow_volume[m3]\tmoving_volume[m3]\tqx[m2/s]\tqy[m2/s]"<<endl;
            logini=1;
        }
        
        logout<<p->simtime<<"\t"<<np<<"\t"<<vol<<"\t"<<vout<<"\t"<<vmov<<"\t"<<qx/(Lx*Ly)<<"\t"<<qy/(Lx*Ly)<<endl;
    }
}
