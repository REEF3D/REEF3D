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
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

// seed parcels in the cells below the initial bed (topo<0) and in the boxes Q 110 (suspension)
// Q 29 1: regular lattice, 2: lattice with jitter, 3: random
// the parcel factor is set such that the seeded bed has the solid fraction 1-S 24
void CPM::seed_particles(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    const int irand(100000);
    const double drand(100000.0);
    
    epsi = p->psi;
    
    // hotstart: the parcels come from the state file
    if(restored==1)
    return;
    
    int dim = p->j_dir==1 ? 3 : 2;
    int ppd = MAX(1, int(pow(double(MAX(p->Q24,1)), 1.0/double(dim)) + 0.5));
    int ppdy = p->j_dir==1 ? ppd : 1;
    int nppc = p->Q29==3 ? MAX(p->Q24,1) : ppd*ppd*ppdy;
    
    if(p->mpirank==0)
    {
        if(p->Q24<1)
        cout<<"CPM warning: Q 24 (parcels per cell) not set, using 1"<<endl;
        
        if(p->Q29!=3 && nppc!=p->Q24)
        cout<<"CPM: lattice seeding with "<<ppd<<" parcels per direction, "<<nppc<<" parcels per cell (Q 24 "<<p->Q24<<")"<<endl;
    }
    
    // cells to seed: below the initial bed or with the centre in a box Q 110
    auto seedcell = [&]()
    {
        // not inside solid bodies (their parcels are removed below; they must not count for the parcel volume)
        if(p->solidread>0 && a->solid(i,j,k)<0.0)
        return false;
        
        if(a->topo(i,j,k)<=0.0)
        return true;
        
        for(int qn=0;qn<p->Q110;++qn)
        if(p->XP[IP]>=p->Q110_xs[qn] && p->XP[IP]<p->Q110_xe[qn]
        && (p->j_dir==0 || (p->YP[JP]>=p->Q110_ys[qn] && p->YP[JP]<p->Q110_ye[qn]))
        && p->ZP[KP]>=p->Q110_zs[qn] && p->ZP[KP]<p->Q110_ze[qn])
        return true;
        
        return false;
    };
    
    // non-uniform grids: all parcels have the same volume, so a cell gets parcels in proportion to its
    // size; Q 24 holds for the smallest seeded cell, larger cells get ppd*dx/dx_min parcels per
    // direction (rounded; stretched grids should use cell sizes that are multiples of the smallest)
    double dxr=1.0e20, dyr=1.0e20, dzr=1.0e20;
    
    BASELOOP
    if(seedcell())
    {
        dxr = MIN(dxr, p->DXN[IP]);
        dyr = MIN(dyr, p->DYN[JP]);
        dzr = MIN(dzr, p->DZN[KP]);
    }
    
    dxr = pgc->globalmin(dxr);
    dyr = pgc->globalmin(dyr);
    dzr = pgc->globalmin(dzr);
    
    auto npx = [&](){return MAX(1, int(double(ppd)*p->DXN[IP]/dxr + 0.5));};
    auto npy = [&](){return p->j_dir==1 ? MAX(1, int(double(ppdy)*p->DYN[JP]/dyr + 0.5)) : 1;};
    auto npz = [&](){return MAX(1, int(double(ppd)*p->DZN[KP]/dzr + 0.5));};
    auto npc = [&](){return p->Q29==3 ? MAX(1, int(double(nppc)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]/(dxr*(p->j_dir==1?dyr:p->DYN[JP])*dzr) + 0.5)) : npx()*npy()*npz();};
    
    // number of parcels, bed volume, largest deviation of the packing from uniform
    int count=0, npar=0;
    double volsum=0.0, dev=0.0;
    
    BASELOOP
    if(seedcell())
    {
        ++count;
        npar += npc();
        volsum += p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        dev = MAX(dev, fabs(double(npc())*dxr*(p->j_dir==1?dyr:p->DYN[JP])*dzr/(double(nppc)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]) - 1.0));
    }
    
    double cellsum = double(pgc->globalsum(count));
    double parsum = double(pgc->globalsum(npar));
    volsum = pgc->globalsum(volsum);
    dev = pgc->globalmax(dev);
    
    if(cellsum>0.0)
    P.ParcelFactor = (p->Q45>0.0 ? p->Q45 : 1.0-p->S24)*volsum/(parsum*Vp);
    
    if(p->mpirank==0)
    {
    cout<<"CPM ParcelFactor: "<<P.ParcelFactor<<" particles per parcel, parcels per cell: "<<nppc<<endl;
    
    if(dev>0.02)
    cout<<"CPM warning: non-uniform grid, the packing of the seeded cells deviates by up to "<<int(100.0*dev+0.5)<<"% (cell sizes not multiples of the smallest)"<<endl;
    }

    P.resize(p,int(double(npar+100*nppc)*MAX(p->Q25,1.0)));
    
    double xs,ys,zs;
    
    BASELOOP
    if(seedcell())
    {
        if(p->Q29==1 || p->Q29==2)
        {
            const int nx=npx(), ny=npy(), nz=npz();
            
            for(int ii = 0; ii < nx; ++ii)
            for(int jj = 0; jj < ny; ++jj)
            for(int kk = 0; kk < nz; ++kk)
            {
                xs = p->XN[IP] + (double(ii)+0.5)*p->DXN[IP]/double(nx);
                ys = p->j_dir==1 ? p->YN[JP] + (double(jj)+0.5)*p->DYN[JP]/double(ny) : p->YP[JP];
                zs = p->ZN[KP] + (double(kk)+0.5)*p->DZN[KP]/double(nz);
                
                if(p->Q29==2)
                {
                    xs += 0.5*(double(rand() % irand)/drand-0.5)*p->DXN[IP]/double(nx);
                    zs += 0.5*(double(rand() % irand)/drand-0.5)*p->DZN[KP]/double(nz);
                    
                    if(p->j_dir==1)
                    ys += 0.5*(double(rand() % irand)/drand-0.5)*p->DYN[JP]/double(ny);
                }
                
                if(P.index_empty<=0)
                break;
                
                --P.index_empty;
                n=P.Empty[P.index_empty];
                
                P.X[n] = xs;
                P.Y[n] = ys;
                P.Z[n] = zs;
                P.U[n] = P.V[n] = P.W[n] = 0.0;

                P.D[n] = seed_diameter(p,n);
                P.RO[n] = p->S22;

                P.Flag[n] = ACTIVE;
                P.Hop[n] = 0.0;
            }
        }
        
        if(p->Q29==3)
        for(int qn=0;qn<npc();++qn)
        {
            if(P.index_empty<=0)
            break;
            
            --P.index_empty;
            n=P.Empty[P.index_empty];
            
            P.X[n] = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
            P.Y[n] = p->j_dir==1 ? p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand : p->YP[JP];
            P.Z[n] = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
            P.U[n] = P.V[n] = P.W[n] = 0.0;

            P.D[n] = seed_diameter(p,n);
            P.RO[n] = p->S22;

            P.Flag[n] = ACTIVE;
            P.Hop[n] = 0.0;
        }
    }
    
    // remove above bed and inside solids
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        double topoval  = p->ccipol4_b(a->topo,P.X[n],P.Y[n],P.Z[n]);
        double solidval = p->ccipol4_b(a->solid,P.X[n],P.Y[n],P.Z[n]);

        bool box=false;
        
        for(int qn=0;qn<p->Q110;++qn)
        if(P.X[n]>=p->Q110_xs[qn] && P.X[n]<p->Q110_xe[qn]
        && (p->j_dir==0 || (P.Y[n]>=p->Q110_ys[qn] && P.Y[n]<p->Q110_ye[qn]))
        && P.Z[n]>=p->Q110_zs[qn] && P.Z[n]<p->Q110_ze[qn])
        box=true;

        if((topoval>0.0 && !box) || solidval<0.0)
        P.remove(n);
    }
    
    for(n=0;n<P.index;++n)
    {
        P.XRK1[n] = P.X[n];
        P.YRK1[n] = P.Y[n];
        P.ZRK1[n] = P.Z[n];
        P.URK1[n] = P.VRK1[n] = P.WRK1[n] = 0.0;
    }
    
    wallbc(p,pgc,s);
    
    // all parcels fixed (Q 44 1)
    if(p->Q44==1)
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    P.Flag[n]=BEDBC;
}

void CPM::ini_fields(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    grid_update(p,a,pgc,s,P.X,P.Y,P.Z,P.U,P.V,P.W);
    count_particles(p,a,pgc,s);
    
    // bed seen by the parcels with the bedload layer (Q 58)
    if(p->Q58>0 && p->S10==1 && zsplit==0)
    topo_iso(p,pgc,Tiso);
}

// parcel diameter: S 20, or for a mixture (S 51 d fa fs) a fraction drawn with the volume
// fractions fa of the bed; all parcels have the same volume, so the number of parcels of a
// fraction follows its volume fraction. A golden ratio sequence gives an even mixture.
double CPM::seed_diameter(lexer *p, int nn)
{
    if(p->S51<=0)
    return p->S20;
    
    double sum=0.0;
    for(int q=0;q<p->S51;++q)
    sum += MAX(p->S51_fa[q],0.0);
    
    if(sum<=0.0)
    return p->S20;
    
    static long seq=0;
    double u = fmod(0.6180339887498949*double(seq + 7919*p->mpirank) + 0.5, 1.0);
    ++seq;
    
    double cum=0.0;
    for(int q=0;q<p->S51;++q)
    {
        cum += MAX(p->S51_fa[q],0.0)/sum;
        
        if(u<cum)
        return p->S51_d[q];
    }
    
    return p->S51_d[p->S51-1];
}
