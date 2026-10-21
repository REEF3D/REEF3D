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
    
    // version 2: the hop length of the bedload layer (Q 58) follows the 8 values of a parcel;
    // version 3: and the parcel weight (merge and split, Q 78), written only if a parcel of the file has a weight
    // other than 1 (the files of runs without merged parcels stay as they were)
    int wtd=0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE && P.Wt[n]!=1.0)
    wtd=1;
    
    int version = wtd ? 3 : 2;
    result.write((char*)&version, sizeof(int));
    result.write((char*)&numpt, sizeof(int));
    result.write((char*)&P.ParcelFactor, sizeof(double));
    result.write((char*)&outvol, sizeof(double));
    
    double val[10];
    int flag;
    const int nwr = version>=3 ? 10 : 9;

    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        val[0]=P.X[n]; val[1]=P.Y[n]; val[2]=P.Z[n];
        val[3]=P.U[n]; val[4]=P.V[n]; val[5]=P.W[n];
        val[6]=P.D[n]; val[7]=P.RO[n]; val[8]=P.Hop[n];
        val[9]=P.Wt[n];
        flag=P.Flag[n];

        result.write((char*)val, nwr*sizeof(double));
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
        cout<<"CPM hotstart: no parcel state file found, the parcels are seeded from the bed of the state"<<endl;

        // the parcels seeded from the grid bed before the state was read are removed: ini_cfd seeds again
        // from the bed of the state (it added a second set of parcels to the first before)
        for(n=0;n<P.index;++n)
        P.Flag[n]=EMPTY;

        P.index_empty=0;
        for(n=P.index-1;n>=0;--n)
        {
            P.Empty[P.index_empty]=n;
            ++P.index_empty;
        }

        // the bed level of the layer starts again from the new parcels (the relaxed level seen by the fluid
        // kept the grid bed and pulled the fluid bed back to it)
        zbf_ini=0;
        zbl_ok=0;

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
    
    double val[10];
    int flag;
    const int nval = version>=3 ? 10 : (version>=2 ? 9 : 8);
    val[8]=0.0;
    val[9]=1.0;

    for(int q=0;q<numpt;++q)
    {
        result.read((char*)val, nval*sizeof(double));
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
        P.Hop[n]=val[8];
        P.Wt[n]=val[9];
        
        // fixed parcels are released unless all parcels stay fixed (Q 44 1),
        // the open boundaries fix theirs again below
        if(flag==BEDBC)
        P.Flag[n] = p->Q44==1 ? BEDBC : ACTIVE;
    }
    
    result.close();
    
    if(p->Q44!=1)
    wallbc(p,pgc,nullptr);
    
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
    
    double vol=0.0, qx=0.0, qy=0.0, vmov=0.0, vbl=0.0, qbx=0.0;
    int np=0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        const double vp = vpar*P.Wt[n];
        
        ++np;
        vol += vp;
        qx += vp*P.U[n];
        qy += vp*P.V[n];
        
        if(P.U[n]*P.U[n] + P.V[n]*P.V[n] + P.W[n]*P.W[n] > 1.0e-6)
        vmov += vp;
        
        if(P.Hop[n]>0.0)
        {
            vbl += vp;
            qbx += vp*P.U[n];
        }
    }
    
    // mean bed shear stress of the bedload layer (Q 58)
    double tb=0.0, nb=0.0;

    if(p->S10!=2)
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        tb += sqrt(blTx(i,j)*blTx(i,j) + blTy(i,j)*blTy(i,j));
        nb += 1.0;
    }

    tb = pgc->globalsum(tb);
    nb = pgc->globalsum(nb);
    tb = nb>0.0 ? tb/nb : 0.0;

    np = pgc->globalisum(np);
    vol = pgc->globalsum(vol);
    qx = pgc->globalsum(qx);
    qy = pgc->globalsum(qy);
    vmov = pgc->globalsum(vmov);
    vbl = pgc->globalsum(vbl);
    qbx = pgc->globalsum(qbx);
    int npick = pgc->globalisum(bl_npick);
    int ndep = pgc->globalisum(bl_ndep);
    int nsus = pgc->globalisum(bl_nsus);
    int nset = pgc->globalisum(bl_nset);
    int nsl = pgc->globalisum(bl_nslide);
    double nslp = pgc->globalsum(double(bl_nsleep));
    double nupd = pgc->globalsum(double(bl_nupd));
    bl_npick = bl_ndep = bl_nsus = bl_nset = bl_nslide = bl_nsleep = bl_nupd = 0;
    
    // sleeping grains of the jammed bed (bed_asleep): share of the parcel updates (parcels not in a hop)
    // skipped since the last line
    if(nupd>0.0 && p->mpirank==0)
    {
        std::streamsize prec = cout.precision();
        cout<<"CPM sleeping grains: "<<setprecision(3)<<100.0*nslp/nupd<<" % of the parcel updates"<<endl;
        cout.precision(prec);
    }
    
    // parcel merge and split (Q 78): since the last line, on the screen
    if(p->Q78>0)
    {
        int nmg = pgc->globalisum(rs_nmerge);
        int nsp = pgc->globalisum(rs_nsplit);
        rs_nmerge = rs_nsplit = 0;
        
        if(p->mpirank==0)
        cout<<"CPM merge and split (Q 78): "<<nmg<<" merged, "<<nsp<<" split, "<<np<<" parcels"<<endl;
    }
    
    // sand slide (Q 75): moved grains since the last line, on the screen (the log keeps its columns)
    if(p->Q75>0.0 && nsl>0 && p->mpirank==0)
    cout<<"CPM sand slide: "<<nsl<<" parcels"<<endl;
    
    // hybrid suspension (Q 58 3): sediment volume in the concentration
    double vsus = afdm!=nullptr ? susp_volume(p,afdm,pgc) : 0.0;
    double vout = pgc->globalsum(outvol);
    
    double Lx = p->global_xmax-p->global_xmin;
    double Ly = p->j_dir==1 ? p->global_ymax-p->global_ymin : p->DYN[marge];
    
    if(p->mpirank==0)
    {
        if(logini==0)
        {
            mkdir("./REEF3D_CFD_CPM_Particle",0777);
            logout.open("./REEF3D_CFD_CPM_Particle/REEF3D-CFD-CPM-Log.dat");
            logout<<"time\tparcels\tsediment_volume[m3]\toutflow_volume[m3]\tmoving_volume[m3]\tqx[m2/s]\tqy[m2/s]\tbedload_volume[m3]\tpickups\tdeposits\treleases\ttau_b[Pa]\tqx_bedload[m2/s]\tsuspended_volume[m3]\tsettled_from_suspension"<<endl;
            logini=1;
        }
        
        logout<<p->simtime<<"\t"<<np<<"\t"<<vol<<"\t"<<vout<<"\t"<<vmov<<"\t"<<qx/(Lx*Ly)<<"\t"<<qy/(Lx*Ly)<<"\t"<<vbl<<"\t"<<npick<<"\t"<<ndep<<"\t"<<nsus<<"\t"<<tb<<"\t"<<qbx/(Lx*Ly)<<"\t"<<vsus<<"\t"<<nset<<endl;
    }
}
