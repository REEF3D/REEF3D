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

#include"lexer.h"
#include"fdm_nhf.h"
#include"sliceint.h"
#include"ghostcell.h"

// iterate over the precomputed boundary-cell list L (same order as the full sweep)
// entries are (i,j,k,h): h=1 if the cell has a flagged neighbour in x or y,
// h=0 if only the bottom/top neighbour is flagged. For h=0 cells only the
// vertical BC statements (the tail of each loop body) can fire.
#define GCBL_LOOP(T) const std::vector<int> &gcbl_L_ = gcbl_get(p,T,i,j,k); int gcbl_h=0; \
    for(size_t qq_=0; qq_<gcbl_L_.size(); qq_+=4) \
    if((i=gcbl_L_[qq_], j=gcbl_L_[qq_+1], k=gcbl_L_[qq_+2], gcbl_h=gcbl_L_[qq_+3], true))

// boundary-cell lists for the V-type BC sweeps (NHFLOW/FNPF):
// cells of a ULOOP/VLOOP/WLOOP/LOOP/FLOOP that have at least one face
// neighbour flagged <0. Only these cells can satisfy any of the BC branches,
// so iterating over them in the original order gives identical results.
// Rebuilt once per time step (p->count) and whenever flags are rebuilt.
// The lists live in the lexer they were built from (lexer::gcbl_ijk), so every grid - the rank
// grid and each mesh refinement patch - has its own, and a patch's lists go with its lexer.
#include<vector>

void gcbl_reset_all(lexer *p)
{
    for(int n=0; n<8; ++n)
    p->gcbl_count[n]=-2;
}

static void gcbl_build_impl(lexer *p, int type, int &i, int &j, int &k)
{
    // during initialisation (count==0) flags may still change: always rebuild
    if(p->gcbl_count[type]==p->count && p->count>0)
    return;

    std::vector<int> &L = p->gcbl_ijk[type];
    L.clear();

    if(type==1)
    ULOOP
    {
    const int h = (p->flag1[Im1JK]<0 || p->flag1[Ip1JK]<0 || p->flag1[IJm1K]<0 || p->flag1[IJp1K]<0) ? 1 : 0;
    if(h==1 || p->flag1[IJKm1]<0 || p->flag1[IJKp1]<0)
    {L.push_back(i); L.push_back(j); L.push_back(k); L.push_back(h);}
    }

    if(type==2)
    VLOOP
    {
    const int h = (p->flag2[Im1JK]<0 || p->flag2[Ip1JK]<0 || p->flag2[IJm1K]<0 || p->flag2[IJp1K]<0) ? 1 : 0;
    if(h==1 || p->flag2[IJKm1]<0 || p->flag2[IJKp1]<0)
    {L.push_back(i); L.push_back(j); L.push_back(k); L.push_back(h);}
    }

    if(type==3)
    WLOOP
    {
    const int h = (p->flag3[Im1JK]<0 || p->flag3[Ip1JK]<0 || p->flag3[IJm1K]<0 || p->flag3[IJp1K]<0) ? 1 : 0;
    if(h==1 || p->flag3[IJKm1]<0 || p->flag3[IJKp1]<0)
    {L.push_back(i); L.push_back(j); L.push_back(k); L.push_back(h);}
    }

    if(type==4)
    LOOP
    {
    const int h = (p->flag4[Im1JK]<0 || p->flag4[Ip1JK]<0 || p->flag4[IJm1K]<0 || p->flag4[IJp1K]<0) ? 1 : 0;
    if(h==1 || p->flag4[IJKm1]<0 || p->flag4[IJKp1]<0)
    {L.push_back(i); L.push_back(j); L.push_back(k); L.push_back(h);}
    }

    if(type==7)
    FLOOP
    {
    const int h = (p->flag7[FIm1JK]<0 || p->flag7[FIp1JK]<0 || p->flag7[FIJm1K]<0 || p->flag7[FIJp1K]<0) ? 1 : 0;
    if(h==1 || p->flag7[FIJKm1]<0 || p->flag7[FIJKp1]<0)
    {L.push_back(i); L.push_back(j); L.push_back(k); L.push_back(h);}
    }

    p->gcbl_count[type]=p->count;
}

static const std::vector<int>& gcbl_get(lexer *p, int type, int &i, int &j, int &k)
{
    gcbl_build_impl(p,type,i,j,k);
    return p->gcbl_ijk[type];
}

void ghostcell::start1V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    if(do_comms)
    gcparaxV1(p, f, gcv);
    if(do_comms)
    gcparacoxV1(p, f, gcv);
    

    // inflow outflow logic

    int inflow=0;
    int outflow=0;

    if(p->B60>=1)
        inflow=1;

    if(p->B98>=3)
        inflow=1;

    if(p->B99>=3)
        outflow=1;

    if(p->B60>=1)
        outflow=1;

    // iowave Riemann / Flather edges (B 520 method 3, 4)
    if(p->open_xm==1)
        inflow=1;

    if(p->open_xp==1)
        outflow=1;

    // 10 U
    // 11 V
    // 12 W
    // 14 ETA
    starttime=timer();
    GCBL_LOOP(1)
    {
    if(gcbl_h==1)
    {
        // s
        // U
        if(p->flag1[Im1JK]<0 && gcv==10 && inflow==0)
            f[Im1JK] = 0.5*fabs(p->W22)*d->eta(i-1,j)*d->eta(i-1,j) + fabs(p->W22)*d->eta(i-1,j)*d->dfx(i,j);

        if(p->flag1[Im1JK]<0 && gcv==10 && inflow==1)
            f[Im1JK] = d->UH[Im1JK]*d->U[Im1JK] + 0.5*fabs(p->W22)*d->eta(i-1,j)*d->eta(i-1,j) + fabs(p->W22)*d->eta(i-1,j)*d->dfx(i,j);

        // V
        if(p->flag1[Im1JK]<0 && gcv==11 && inflow==0)
            f[Im1JK] = 0.0;

        if(p->flag1[Im1JK]<0 && gcv==11 && inflow>=1)
            f[Im1JK] = d->VH[Im1JK]*d->U[Im1JK];

        // W
        if(p->flag1[Im1JK]<0 && gcv==12 && inflow==0)
            f[Im1JK] = 0.0;

        if(p->flag1[Im1JK]<0 && gcv==12 && inflow>=1)
            f[Im1JK] = d->WH[Im1JK]*d->U[Im1JK];

        // ETA
        if(p->flag1[Im1JK]<0 && gcv==14 && inflow==0)
            f[Im1JK] = 0.0;

        if(p->flag1[Im1JK]<0 && gcv==14 && inflow>=1)
            f[Im1JK] = d->UH[Im1JK];

        // n
        // U
        if(p->flag1[Ip1JK]<0 && gcv==10 && outflow==0)
            f[Ip1JK] = 0.5*fabs(p->W22)*d->eta(i+2,j)*d->eta(i+2,j) + fabs(p->W22)*d->eta(i+2,j)*d->dfx(i+1,j);

        if(p->flag1[Ip1JK]<0 && gcv==10 && outflow==1)
            f[Ip1JK] = d->UH[Ip2JK]*d->U[Ip2JK] + 0.5*fabs(p->W22)*d->eta(i+2,j)*d->eta(i+2,j) + fabs(p->W22)*d->eta(i+2,j)*d->dfx(i+1,j);

        // V
        if(p->flag1[Ip1JK]<0 && gcv==11 && outflow==0)
            f[Ip1JK] = 0.0;

        if(p->flag1[Ip1JK]<0 && gcv==11 && outflow==1)
            f[Ip1JK] = d->VH[Ip2JK]*d->U[Ip2JK];

        // W
        if(p->flag1[Ip1JK]<0 && gcv==12 && outflow==0)
            f[Ip1JK] = 0.0;

        if(p->flag1[Ip1JK]<0 && gcv==12 && outflow==1)
            f[Ip1JK] = d->WH[Ip2JK]*d->U[Ip2JK];

        // ETA
        if(p->flag1[Ip1JK]<0 && gcv==14 && outflow==0)
            f[Ip1JK] = 0.0;

        if(p->flag1[Ip1JK]<0 && gcv==14 && outflow==1)
            f[Ip1JK] = d->UH[Ip2JK];

        // e
        if(p->flag1[IJm1K]<0 && p->j_dir==1)
            f[IJm1K] = 0.0;

        // w
        if(p->flag1[IJp1K]<0 && p->j_dir==1)
        f[IJp1K] = 0.0;

        // b
        if(p->flag1[IJKm1]<0)
            f[IJKm1] = 0.0;

        // t
        if(p->flag1[IJKp1]<0)
            f[IJKp1] = 0.0;
    }
    else
    {
        if(p->flag1[IJKm1]<0)
            f[IJKm1] = 0.0;

        // t
        if(p->flag1[IJKp1]<0)
            f[IJKp1] = 0.0;
    
    }
    }
    p->gctime+=timer()-starttime;
}

void ghostcell::start2V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    if(do_comms)
    gcparaxV1(p, f, gcv);
    if(do_comms)
    gcparacoxV1(p, f, gcv);
    
    int inflow=0;
    int outflow=0;

    if(p->B98>=3 || p->B60>=1)
    inflow=1;

    if(p->B99>=3)
    outflow=1;

    if(p->B60>=1)
    outflow=1;

    // 10 U
    // 11 V
    // 12 W
    // 14 ETA
    starttime=timer();
    GCBL_LOOP(2)
    {
    if(gcbl_h==1)
    {
    // s
        if(p->flag2[Im1JK]<0)
        {
        f[Im1JK] = 0.0;
        }

    // n
        if(p->flag2[Ip1JK]<0)
        {
        f[Ip1JK] = 0.0;
        }

    // e
        // iowave Riemann / Flather edge on y- (ghost cells set by iowave): fluxes from the ghost cell
        const int openm = (p->open_ym==1 && p->origin_j+j==0) ? 1 : 0;
        const int openp = (p->open_yp==1 && p->origin_j+j+1==p->gknoy-1) ? 1 : 0;
        
        if(p->flag2[IJm1K]<0 && p->j_dir==1 && openm==1)
        {
        if(gcv==11)
        f[IJm1K] = d->VH[IJm1K]*d->V[IJm1K] + 0.5*fabs(p->W22)*d->eta(i,j-1)*d->eta(i,j-1) + fabs(p->W22)*d->eta(i,j-1)*d->dfy(i,j);
        
        if(gcv==14)
        f[IJm1K] = d->VH[IJm1K];
        
        if(gcv==10)
        f[IJm1K] = d->UH[IJm1K]*d->V[IJm1K];
        
        if(gcv==12)
        f[IJm1K] = d->WH[IJm1K]*d->V[IJm1K];
        }
        
        // V
        if(p->flag2[IJm1K]<0 &&  gcv==11 && p->j_dir==1 && openm==0)
        {
        f[IJm1K] = 0.5*fabs(p->W22)*d->eta(i,j-1)*d->eta(i,j-1) + fabs(p->W22)*d->eta(i,j-1)*d->dfy(i,j);
        }

        // ETA
        if(p->flag2[IJm1K]<0 &&  gcv==14 && p->j_dir==1 && openm==0)
        {
        f[IJm1K] = 0.0;
        }

        // U,W
        if(p->flag2[IJm1K]<0 && gcv!=11 && gcv!=14 && p->j_dir==1 && openm==0)
        {
        f[IJm1K] = 0.0;
        }

    // w
        // iowave Riemann / Flather edge on y+: fluxes from the first ghost cell
        if(p->flag2[IJp1K]<0 && p->j_dir==1 && openp==1)
        {
        if(gcv==11)
        f[IJp1K] = d->VH[IJp2K]*d->V[IJp2K] + 0.5*fabs(p->W22)*d->eta(i,j+2)*d->eta(i,j+2) + fabs(p->W22)*d->eta(i,j+2)*d->dfy(i,j+1);
        
        if(gcv==14)
        f[IJp1K] = d->VH[IJp2K];
        
        if(gcv==10)
        f[IJp1K] = d->UH[IJp2K]*d->V[IJp2K];
        
        if(gcv==12)
        f[IJp1K] = d->WH[IJp2K]*d->V[IJp2K];
        }
        
        // V
        if(p->flag2[IJp1K]<0 &&  gcv==11 && p->j_dir==1 && openp==0)
        {
        f[IJp1K] = 0.5*fabs(p->W22)*d->eta(i,j+2)*d->eta(i,j+2) + fabs(p->W22)*d->eta(i,j+2)*d->dfy(i,j+1);
        }

        // ETA
        if(p->flag2[IJp1K]<0 &&  gcv==14 && p->j_dir==1 && openp==0)
        {
        f[IJp1K] = 0.0;
        }

        // U,W
        if(p->flag2[IJp1K]<0 && gcv!=11 && gcv!=14 && p->j_dir==1 && openp==0)
        {
        f[IJp1K] = 0.0;
        }

    // b
        if(p->flag2[IJKm1]<0 && p->j_dir==1)
        {
        f[IJKm1] = 0.0;
        }

    // t
        if(p->flag2[IJKp1]<0 && p->j_dir==1)
        {
        f[IJKp1] = 0.0;
        }
    }
    else
    {
        if(p->flag2[IJKm1]<0 && p->j_dir==1)
        {
        f[IJKm1] = 0.0;
        }

    // t
        if(p->flag2[IJKp1]<0 && p->j_dir==1)
        {
        f[IJKp1] = 0.0;
        }
    
    }
    }
    p->gctime+=timer()-starttime;
}

void ghostcell::start3V(lexer *p, double *f, int gcv)
{
    if(do_comms)
    gcparaxV1(p, f, gcv);
    if(do_comms)
    gcparacoxV1(p, f, gcv);

    int inflow=0;
    int outflow=0;

    if(p->B98>=3 || p->B60>=1)
    inflow=1;

    if(p->B99>=3)
    outflow=1;

    if(p->B60>=1)
    outflow=1;

    // 10 U
    // 11 V
    // 12 W
    // 14 ETA
    starttime=timer();
    GCBL_LOOP(3)
    {
    if(gcbl_h==1)
    {
        if(p->flag3[Im1JK]<0)
        {
        f[Im1JK] = 0.0;
        }

        if(p->flag3[Ip1JK]<0)
        {
        f[Ip1JK] = 0.0;
        }

        if(p->flag3[IJm1K]<0)
        {
        f[IJm1K] = 0.0;
        }

        if(p->flag3[IJp1K]<0)
        {
        f[IJp1K] = 0.0;
        }

        if(p->flag3[IJKm1]<0)
        {
        f[IJKm1] = 0.0;
        }

        if(p->flag3[IJKp1]<0 && gcv!=10)
        {
        f[IJKp1] = 0.0;
        }

        if(p->flag3[IJKp1]<0 && gcv==10)
        {
        f[IJKp1] = 0.0;
        }
    }
    else
    {
        if(p->flag3[IJKm1]<0)
        {
        f[IJKm1] = 0.0;
        }

        if(p->flag3[IJKp1]<0 && gcv!=10)
        {
        f[IJKp1] = 0.0;
        }

        if(p->flag3[IJKp1]<0 && gcv==10)
        {
        f[IJKp1] = 0.0;
        }
    
    }
    }
    p->gctime+=timer()-starttime;
}

void ghostcell::start4V_par(lexer *p, double *f, int gcv)
{
    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);
}

void ghostcell::start4V(lexer *p, double *f, int gcv)
{
    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);

    int inflow=0;
    int outflow=0;

    if(p->B98>=3 || p->B60>=1)
        inflow=1;

    if(p->B99>=3)
        outflow=1;

    if(p->B60>=1)
        outflow=2;

    // waves on a current with relaxation wave generation (B 98 2, B 60 1): the relaxation zone
    // prescribes the velocity, the inflow ghost cells take it over (zero gradient) as without current
    if(p->B98==2 && p->B60>=1)
        inflow=0;

    // iowave Riemann / Flather edges (B 520 method 3, 4): ghost cells set by iowave
    if(p->open_xm==1)
        inflow=1;

    if(p->open_xp==1)
        outflow=1;

    starttime=timer();
    GCBL_LOOP(4)
    {
    if(gcbl_h==1)
    {
        // xxxxxxx
        // s
        if(p->flag4[Im1JK]<0 && (gcv==10 || gcv==14) && (inflow==0 && p->B98!=2))
        {
            f[Im1JK] = 0.0;
            f[Im2JK] = 0.0;
            f[Im3JK] = 0.0;
        }
        
        if(p->flag4[Im1JK]<0 && (gcv==10 || gcv==14) && (inflow==0 && p->B98==2))
        {
            f[Im1JK] = f[IJK];
            f[Im2JK] = f[IJK];
            f[Im3JK] = f[IJK];
        }

        if(p->flag4[Im1JK]<0 && (gcv!=10 && gcv!=14) && inflow==0)
        {
            f[Im1JK] = 0.0;
            f[Im2JK] = 0.0;
            f[Im3JK] = 0.0;
        }

        // n
        if(p->flag4[Ip1JK]<0 && (gcv==10 || gcv==14) && outflow==0)
        {
            f[Ip1JK] = 0.0;
            f[Ip2JK] = 0.0;
            f[Ip3JK] = 0.0;
        }

        if(p->flag4[Ip1JK]<0 && (gcv==10 || gcv==14) && outflow==2)
        {
            f[Ip1JK] = MAX(0.0,f[IJK] - p->dt/p->DXP[IP]*p->Uo*(f[IJK]-f[Im1JK]));
            f[Ip2JK] = MAX(0.0,f[IJK]);
            f[Ip3JK] = MAX(0.0,f[IJK]);
        }

        if(p->flag4[Ip1JK]<0 && (gcv!=10 && gcv!=14) && outflow==0)
        {
            f[Ip1JK] = 0.0;
            f[Ip2JK] = 0.0;
            f[Ip3JK] = 0.0;
        }

        if(p->flag4[Ip1JK]<0 && (gcv!=10 && gcv!=14) && outflow==2)
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        // yyyyy
        // side walls are slip walls: the normal velocity (V, VH: gcv 11, 15) is mirrored
        // antisymmetrically, all other fields (U, W, UH, WH, scalars) symmetrically.
        // (Setting the tangential components to 0 acted like a no-slip wall and, through
        // the reconstruction and the HLL dissipation, damped waves in 3D.)
        if(p->flag4[IJm1K]<0 && p->j_dir==1 && (gcv==11 || gcv==15) && (p->open_ym==0 || p->origin_j+j>0))
        {
            f[IJm1K] = -f[IJK];
            f[IJm2K] = -f[IJp1K];
            f[IJm3K] = -f[IJp2K];
        }

        if(p->flag4[IJm1K]<0 && p->j_dir==1 && (gcv!=11 && gcv!=15) && (p->open_ym==0 || p->origin_j+j>0))
        {
            f[IJm1K] = f[IJK];
            f[IJm2K] = f[IJp1K];
            f[IJm3K] = f[IJp2K];
        }

        if(p->flag4[IJp1K]<0 && p->j_dir==1 && (gcv==11 || gcv==15) && (p->open_yp==0 || p->origin_j+j<p->gknoy-1))
        {
            f[IJp1K] = -f[IJK];
            f[IJp2K] = -f[IJm1K];
            f[IJp3K] = -f[IJm2K];
        }

        if(p->flag4[IJp1K]<0 && p->j_dir==1 && (gcv!=11 && gcv!=15) && (p->open_yp==0 || p->origin_j+j<p->gknoy-1))
        {
            f[IJp1K] = f[IJK];
            f[IJp2K] = f[IJm1K];
            f[IJp3K] = f[IJm2K];
        }

        // zzzzz
        if(p->flag4[IJKp1]<0 && (gcv==14))
        {
            f[IJKp1] = 0.0;
            f[IJKp2] = 0.0;
            f[IJKp3] = 0.0;
        }

        if(p->flag4[IJKp1]<0 && (gcv==10||gcv==11||gcv==14||gcv==15))
        {
            f[IJKp1] = f[IJK];
            f[IJKp2] = f[IJK];
            f[IJKp3] = f[IJK];
        }

        if(p->flag4[IJKp1]<0 && gcv==12)
        {
            f[IJKp1] = f[IJK];
            f[IJKp2] = f[IJK];
            f[IJKp3] = f[IJK];
        }

        // bed
        if(p->flag4[IJKm1]<0 && (gcv==10||gcv==11||gcv==14||gcv==15))
        {
            if(p->A518==1)
            {
                f[IJKm1] = f[IJK];
                f[IJKm2] = f[IJK];
                f[IJKm3] = f[IJK];
            }
             if(p->A518==2)
            {
                f[IJKm1] = 0.0;
                f[IJKm2] = 0.0;
                f[IJKm3] = 0.0;
            }
        }
    }
    else
    {
        if(p->flag4[IJKp1]<0 && (gcv==14))
        {
            f[IJKp1] = 0.0;
            f[IJKp2] = 0.0;
            f[IJKp3] = 0.0;
        }

        if(p->flag4[IJKp1]<0 && (gcv==10||gcv==11||gcv==14||gcv==15))
        {
            f[IJKp1] = f[IJK];
            f[IJKp2] = f[IJK];
            f[IJKp3] = f[IJK];
        }

        if(p->flag4[IJKp1]<0 && gcv==12)
        {
            f[IJKp1] = f[IJK];
            f[IJKp2] = f[IJK];
            f[IJKp3] = f[IJK];
        }

        // bed
        if(p->flag4[IJKm1]<0 && (gcv==10||gcv==11||gcv==14||gcv==15))
        {
            if(p->A518==1)
            {
                f[IJKm1] = f[IJK];
                f[IJKm2] = f[IJK];
                f[IJKm3] = f[IJK];
            }
             if(p->A518==2)
            {
                f[IJKm1] = 0.0;
                f[IJKm2] = 0.0;
                f[IJKm3] = 0.0;
            }
        }
    
    }
    }

    p->gctime+=timer()-starttime;
}

void ghostcell::start5V(lexer *p, double *f, int gcv)
{
    GCBL_LOOP(4)
    {
    if(gcbl_h==1)
    {
        if(p->flag4[Im1JK]<0)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im3JK] = f[IJK];
        }
        
        //
        if(p->flag4[Ip1JK]<0)
        {
        f[Ip1JK] = f[IJK];
        f[Ip2JK] = f[IJK];
        f[Ip3JK] = f[IJK];
        }
        
        //
        if(p->flag4[IJm1K]<0 && p->j_dir==1)
        {
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        }
        
        //
        if(p->flag4[IJp1K]<0 && p->j_dir==1)
        {
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
        }

        //
        if(p->flag4[IJKm1]<0)
        {
        f[IJKm1] = f[IJK];
        f[IJKm2] = f[IJK];
        f[IJKm3] = f[IJK];
        }

        //
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
    }
    else
    {
        if(p->flag4[IJKm1]<0)
        {
        f[IJKm1] = f[IJK];
        f[IJKm2] = f[IJK];
        f[IJKm3] = f[IJK];
        }

        //
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
    
    }
    }

    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);
}

void ghostcell::start5Vfull(lexer *p, double *f, int gcv)
{
    if(p->j_dir==0)
    LOOP
    {
        
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
    }
        
    LOOP
    {
        if(p->flag4[Im1JK]<0)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im2JK] = f[IJK];
        }
        
        
        if(p->flag4[Ip1JK]<0)
        {
        f[Ip1JK] = f[IJK];
        f[Ip2JK] = f[IJK];
        f[Ip3JK] = f[IJK];
        }
        
        
        if(p->flag4[IJm1K]<0 && p->j_dir==1)
        {
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        }
        
        
        if(p->flag4[IJp1K]<0 && p->j_dir==1)
        {
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
        }

        
        if(p->flag4[IJKm1]<0)
        {
        f[IJKm1] = f[IJK];
        f[IJKm2] = f[IJK];
        f[IJKm3] = f[IJK];
        }

        
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
    }

    if(do_comms)
    gcparaxV(p, f, gcv);
    //gcparacoxV(p, f, gcv);
}

// x- inflow ghost cells of k, eps/omega and nu_t keep their values only where nhflow_rans_io::inflow writes the
// equilibrium profile (RANS with a discharge inflow, B 60 >= 1, IO 1); all other inflows (LES, wave generation
// B 98 >= 3) get zero gradient (the ghosts were never written and stayed 0)
static inline bool nhflow_turb_profile_ghost(lexer *p, int ijk)
{
    return p->B60>=1 && (p->A560==1 || p->A560==21 || p->A560==2 || p->A560==22) && p->IO[ijk]==1;
}

void ghostcell::start20V(lexer *p, double *f, int gcv) //KIN
{
    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);

    int inflow=0;
    int outflow=0;

    if(p->B98>=3 || p->B60>=1)
    inflow=1;

    if(p->B99>=3)
    outflow=1;

    if(p->B60>=1)
    outflow=2;

    starttime=timer();
    LOOP
    if(p->DF[IJK]>0)
    {

        // xxxxxxx
        // s
        if((p->flag4[Im1JK]<0 && (inflow==0 || !nhflow_turb_profile_ghost(p,Im1JK))) || (p->DF[Im1JK]<0))
        {
            if(p->B11==1)
            {
                f[Im1JK] = f[IJK];
                f[Im2JK] = f[IJK];
                f[Im3JK] = f[IJK];
            }
            else if(p->B11==0)
            {
                f[Im1JK] = 0.0;
                f[Im2JK] = 0.0;
                f[Im3JK] = 0.0;
            }
        }

        /*if(p->flag4[Im1JK]<0 && inflow==1)
        {
            f[Im1JK] = 0.0;
            f[Im2JK] = 0.0;
            f[Im3JK] = 0.0;
        }*/

        // n
        if((p->flag4[Ip1JK]<0 && outflow==0) || (p->DF[Ip1JK]<0))
        {
            if(p->B11==1)
            {
                f[Ip1JK] = f[IJK];
                f[Ip2JK] = f[IJK];
                f[Ip3JK] = f[IJK];
            }
            else if(p->B11==0)
            {
                f[Ip1JK] = 0.0;
                f[Ip2JK] = 0.0;
                f[Ip3JK] = 0.0;
            }
        }

        if(p->flag4[Ip1JK]<0 && outflow==2)
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        // yyyyy
        if((p->flag4[IJm1K]<0  || (p->DF[IJm1K]<0)) && p->j_dir==1)
        {
            if(p->B11==1)
            {
                f[IJm1K] = f[IJK];
                f[IJm2K] = f[IJK];
                f[IJm3K] = f[IJK];
            }
            else if(p->B11==0)
            {
                f[IJm1K] = 0.0;
                f[IJm2K] = 0.0;
                f[IJm3K] = 0.0;
            }
        }

        if((p->flag4[IJp1K]<0  || (p->DF[IJp1K]<0)) && p->j_dir==1)
        {
            if(p->B11==1)
            {
                f[IJp1K] = f[IJK];
                f[IJp2K] = f[IJK];
                f[IJp3K] = f[IJK];
            }
            else if(p->B11==0)
            {
                f[IJp1K] = 0.0;
                f[IJp2K] = 0.0;
                f[IJp3K] = 0.0;
            }
        }

        // zzzzz
        if(p->flag4[IJKp1]<0  || (p->DF[IJKp1]<0) || k==p->knoz-1)
        {

                f[IJKp1] = f[IJK];
                f[IJKp2] = f[IJK];
                f[IJKp3] = f[IJK];
        }

        // bed
        if(p->flag4[IJKm1]<0  || (p->DF[IJKm1]<0 && p->B11==1)   || (k==0 && p->B11==1))
        {
            if(p->B11==1)
            {
                f[IJKm1] = f[IJK];
                f[IJKm2] = f[IJK];
                f[IJKm3] = f[IJK];
            }
            else if(p->B11==0)
            {
                f[IJKm1] = 0.0;
                f[IJKm2] = 0.0;
                f[IJKm3] = 0.0;
            }
        }
    }

    if(do_comms)
    gcparacoxV(p, f, gcv);

    p->gctime+=timer()-starttime;
}

void ghostcell::start24V(lexer *p, double *f, int gcv) //EDDYV
{
    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);

    int inflow=0;
    int outflow=0;

    if(p->B98>=3 || p->B60>=1)
        inflow=1;

    if(p->B99>=3)
        outflow=1;

    if(p->B60>=1)
        outflow=2;

    starttime=timer();
    LOOP
    if(p->DF[IJK]>0)
    {
        // xxxxxxx
        // s
        if(p->flag4[Im1JK]<0 && (inflow==0 || !nhflow_turb_profile_ghost(p,Im1JK)))
        {
            f[Im1JK] = f[IJK];
            f[Im2JK] = f[IJK];
            f[Im3JK] = f[IJK];
        }

        /*if(p->flag4[Im1JK]<0 && inflow==1)
        {
            f[Im1JK] = 0.0;
            f[Im2JK] = 0.0;
            f[Im3JK] = 0.0;
        }*/

        // n
        if(p->flag4[Ip1JK]<0 && outflow==0)
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        if(p->flag4[Ip1JK]<0 && outflow==1)
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        if(p->flag4[Ip1JK]<0 && outflow==2)
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        // yyyyy
        if(p->flag4[IJm1K]<0 && p->j_dir==1)
        {
            f[IJm1K] = f[IJK];
            f[IJm2K] = f[IJK];
            f[IJm3K] = f[IJK];
        }

        if(p->flag4[IJp1K]<0 && p->j_dir==1)
        {
            f[IJp1K] = f[IJK];
            f[IJp2K] = f[IJK];
            f[IJp3K] = f[IJK];
        }

        // zzzzz
        if(p->flag4[IJKp1]<0 || k==p->knoz-1)
        {
            f[IJKp1] = f[IJK];
            f[IJKp2] = f[IJK];
            f[IJKp3] = f[IJK];
        }

        // bed
        if(p->flag4[IJKm1]<0)
        {
            f[IJKm1] = f[IJK];
            f[IJKm2] = f[IJK];
            f[IJKm3] = f[IJK];
        }
    }

    if(do_comms)
    gcparacoxV(p, f, gcv);

    p->gctime+=timer()-starttime;
}

void ghostcell::start30V(lexer *p, double *f, int gcv) // EPS
{
    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);

    int inflow=0;
    int outflow=0;

    if(p->B98>=3 || p->B60>=1)
        inflow=1;

    if(p->B99>=3)
        outflow=1;

    if(p->B60>=1)
        outflow=2;

    starttime=timer();
    LOOP
    if(p->DF[IJK]>0)
    {
        // xxxxxxx
        // s
        if((p->flag4[Im1JK]<0 && (inflow==0 || !nhflow_turb_profile_ghost(p,Im1JK))) || (p->DF[Im1JK]<0))
        {
            f[Im1JK] = f[IJK];
            f[Im2JK] = f[IJK];
            f[Im3JK] = f[IJK];
        }

        /*if(p->flag4[Im1JK]<0 && inflow==1)
        {
            f[Im1JK] = f[IJK];
            f[Im2JK] = f[IJK];
            f[Im3JK] = f[IJK];
        }*/

        // n
        if((p->flag4[Ip1JK]<0 && outflow==0) || (p->DF[Ip1JK]<0))
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        if(p->flag4[Ip1JK]<0 && outflow==2)
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        // yyyyy
        if((p->flag4[IJm1K]<0  || (p->DF[IJm1K]<0)) && p->j_dir==1)
        {
            f[IJm1K] = f[IJK];
            f[IJm2K] = f[IJK];
            f[IJm3K] = f[IJK];
        }

        if((p->flag4[IJp1K]<0  || (p->DF[IJp1K]<0)) && p->j_dir==1)
        {
            f[IJp1K] = f[IJK];
            f[IJp2K] = f[IJK];
            f[IJp3K] = f[IJK];
        }

        // zzzzz
        if(p->flag4[IJKp1]<0  || (p->DF[IJKp1]<0) || k==p->knoz-1)
        {
            f[IJKp1] = f[IJK];
            f[IJKp2] = f[IJK];
            f[IJKp3] = f[IJK];
        }

        // bed
        if(p->flag4[IJKm1]<0  || (p->DF[IJKm1]<0) || k==0)
        {
            f[IJKm1] = f[IJK];
            f[IJKm2] = f[IJK];
            f[IJKm3] = f[IJK];
        }
    }

    if(do_comms)
    gcparacoxV(p, f, gcv);
}

void ghostcell::start49V(lexer *p, double *f, int gcv)
{
    LOOP
    {
        // if(p->flag4[Im1JK]<0 && p->IO[Im1JK]!=1)
        //     f[Im1JK] = f[IJK];

        if(p->flag4[Im1JK]<0)
            f[Im1JK] = f[IJK] - p->Ui*p->DXP[IM1];   // dpsi/dx = Ui at the inflow, same sign as the Laplace rhs in nhflow_potential_f

        // if(p->flag4[Ip1JK]<0 && p->IO[Ip1JK]!=2)
        //     f[Ip1JK] = f[IJK];

        if(p->flag4[Ip1JK]<0)
            f[Ip1JK] = p->Uo*p->DXP[IP] + f[IJK];

        if(p->flag4[IJm1K]<0)
            f[IJm1K] = f[IJK];

        if(p->flag4[IJp1K]<0)
            f[IJp1K] = f[IJK];

        if(p->flag4[IJKm1]<0)
            f[IJKm1] = f[IJK];

        if(p->flag4[IJKp1]<0)
            f[IJKp1] = f[IJK];
    }

    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);
}

void ghostcell::start60V(lexer *p, double *f, int gcv) // EPS
{
    if(do_comms)
    gcparaxV(p, f, gcv);
    if(do_comms)
    gcparacoxV(p, f, gcv);

    int inflow=0;
    int outflow=0;

    if(p->B98>=3 || p->B60>=1)
        inflow=1;

    if(p->B99>=3)
        outflow=1;

    if(p->B60>=1)
        outflow=2;

    starttime=timer();
    LOOP
    if(p->DF[IJK]>0)
    {
        // xxxxxxx
        // s
        if((p->flag4[Im1JK]<0 && inflow==0) || (p->DF[Im1JK]<0))
        {
            f[Im1JK] = f[IJK];
            f[Im2JK] = f[IJK];
            f[Im3JK] = f[IJK];
        }

        /*if(p->flag4[Im1JK]<0 && inflow==1)
        {
            f[Im1JK] = f[IJK];
            f[Im2JK] = f[IJK];
            f[Im3JK] = f[IJK];
        }*/

        // n
        if((p->flag4[Ip1JK]<0 && outflow==0) || (p->DF[Ip1JK]<0))
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        if(p->flag4[Ip1JK]<0 && outflow==2)
        {
            f[Ip1JK] = f[IJK];
            f[Ip2JK] = f[IJK];
            f[Ip3JK] = f[IJK];
        }

        // yyyyy
        if((p->flag4[IJm1K]<0  || (p->DF[IJm1K]<0)) && p->j_dir==1)
        {
            f[IJm1K] = f[IJK];
            f[IJm2K] = f[IJK];
            f[IJm3K] = f[IJK];
        }

        if((p->flag4[IJp1K]<0  || (p->DF[IJp1K]<0)) && p->j_dir==1)
        {
            f[IJp1K] = f[IJK];
            f[IJp2K] = f[IJK];
            f[IJp3K] = f[IJK];
        }

        // zzzzz
        if(p->flag4[IJKp1]<0  || p->DF[IJKp1]<0 || k==p->knoz-1)
        {
            f[IJKp1] = 0.0;
            f[IJKp2] = 0.0;
            f[IJKp3] = 0.0;
        }

        // bed
        if(p->flag4[IJKm1]<0  || p->DF[IJKm1]<0 || k==0)
        {
            f[IJKm1] = f[IJK];
            f[IJKm2] = f[IJK];
            f[IJKm3] = f[IJK];
        }
    }

    if(do_comms)
    gcparacoxV(p, f, gcv);
}

void ghostcell::startintV(lexer *p, int *f, int gcv)
{
    LOOP
    {
        if(p->flag4[Im1JK]<0)
        f[Im1JK] = f[IJK];

        if(p->flag4[Ip1JK]<0)
        f[Ip1JK] = f[IJK];

        if(p->flag4[IJm1K]<0)
        f[IJm1K] = f[IJK];

        if(p->flag4[IJp1K]<0)
        f[IJp1K] = f[IJK];

        if(p->flag4[IJKm1]<0)
        f[IJKm1] = f[IJK];

        if(p->flag4[IJKp1]<0)
        f[IJKp1] = f[IJK];
    }

    if(do_comms)
    gcparaxintV(p, f, gcv);
}

void ghostcell::startintVF(lexer *p, int *f, int gcv)
{
    FLOOP
    {
        if(p->flag7[FIm1JK]<0)
        f[FIm1JK] = f[FIJK];

        if(p->flag7[FIp1JK]<0)
        f[FIp1JK] = f[FIJK];

        if(p->flag7[FIJm1K]<0)
        f[FIJm1K] = f[FIJK];

        if(p->flag7[FIJp1K]<0)
        f[FIJp1K] = f[FIJK];

        if(p->flag7[FIJKm1]<0)
        f[FIJKm1] = f[FIJK];

        if(p->flag7[FIJKp1]<0)
        f[FIJKp1] = f[FIJK];
    }
    
    
    if(do_comms)
    {
        gcparax7int(p,f,7);
        
    }
}

void ghostcell::start7V(lexer *p, double *f, sliceint &bc, int gcv)
{
    if(do_comms)
    {
        gcparax7(p,f,7);
        gcparax7co(p,f,7);
    }

    
    if(gcv==250)
    fivec(p,f,bc);
        
    else if(gcv==150)
    fivec2D(p,f,bc);
        
    else if(gcv==210)
    fivec_vel(p,f,bc);
        
    else if(gcv==110)
    fivec2D_vel(p,f,bc);
    
}

void ghostcell::start7P(lexer *p, double *f, int gcv)
{
    GCBL_LOOP(7)
    {
    if(gcbl_h==1)
    {
        if(p->flag7[FIm1JK]<0)
        f[FIm1JK] = f[FIJK];

        if(p->flag7[FIp1JK]<0)
        f[FIp1JK] = f[FIJK];

        if(p->flag7[FIJm1K]<0)
        f[FIJm1K] = f[FIJK];

        if(p->flag7[FIJp1K]<0)
        f[FIJp1K] = f[FIJK];

        if(p->flag7[FIJKm1]<0)
        f[FIJKm1] = f[FIJK];

        if(p->flag7[FIJKp1]<0)
        f[FIJKp1] = 0.0;
    }
    else
    {
        if(p->flag7[FIJKm1]<0)
        f[FIJKm1] = f[FIJK];

        if(p->flag7[FIJKp1]<0)
        f[FIJKp1] = 0.0;
    
    }
    }

    if(do_comms)
    {
        gcparax7(p,f,7);
        gcparax7co(p,f,7);
    }
}

void ghostcell::start7S(lexer *p, double *f, int gcv)
{
    GCBL_LOOP(7)
    {
    if(gcbl_h==1)
    {
        if(p->flag7[FIm1JK]<0)
        f[FIm1JK] = f[FIJK];

        if(p->flag7[FIp1JK]<0)
        f[FIp1JK] = f[FIJK];

        if(p->flag7[FIJm1K]<0)
        f[FIJm1K] = f[FIJK];

        if(p->flag7[FIJp1K]<0)
        f[FIJp1K] = f[FIJK];

        if(p->flag7[FIJKm1]<0)
        f[FIJKm1] = f[FIJK];

        if(p->flag7[FIJKp1]<0)
        f[FIJKp1] = f[FIJK];
    }
    else
    {
        if(p->flag7[FIJKm1]<0)
        f[FIJKm1] = f[FIJK];

        if(p->flag7[FIJKp1]<0)
        f[FIJKp1] = f[FIJK];
    
    }
    }

    if(do_comms)
    {
        gcparax7(p,f,7);
        gcparax7co(p,f,7);
    }
}
