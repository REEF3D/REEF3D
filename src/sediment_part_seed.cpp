/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs & Alexander Hanke
--------------------------------------------------------------------*/

#include "sediment_part.h"
#include "partres.h"

#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include <math.h>

using std::cout;
using std::endl;

/// @brief Determines cells to be seeded with particles and maximum number of particles

/// Sets \p sediment_part::active_box and \p sediment_part::active_topo to `true` for cells to be seeded with particles.
/// \p sediment_part::active_box based on \p lexer::Q110 and \p sediment_part::active_topo based on \p fdm::topo.

/// \p lexer::Q24 is used as a guess for particles per cell, if it is not set, it is set to 0.
/// With \p sediment_part::ppcell and the number of cells active for seeding, \p sediment_part::maxparticle is calculated.
void sediment_part::seed_ini(lexer *p, fdm *a, ghostcell* pgc)
{
    // ini
    LOOP
    {
        active_box(i,j,k) = false;
        if(p->Q111>0)
        active_box_dummy(i,j,k) = false;
        active_topo(i,j,k) = false;
    }
    double minPPC=INT64_MAX;
    
    // Box
    size_t cellCountBox=0;
    for(int qn=0;qn<p->Q110;++qn)
    LOOP
    if(p->XN[IP]>=p->Q110_xs[qn] && p->XN[IP1]<=p->Q110_xe[qn]
    && p->YN[JP]>=p->Q110_ys[qn] && p->YN[JP1]<=p->Q110_ye[qn]
    && p->ZN[KP1]>=p->Q110_zs[qn] && p->ZN[KP1]<=p->Q110_ze[qn])
    {
        active_box(i,j,k) = true;
        ++cellCountBox;
    }
    size_t cellCountBoxDummy=0;
    for(int qn=0;qn<p->Q111;++qn)
    LOOP
    if(p->XN[IP]>=p->Q111_xs[qn] && p->XN[IP1]<=p->Q111_xe[qn]
    && p->YN[JP]>=p->Q111_ys[qn] && p->YN[JP1]<=p->Q111_ye[qn]
    && p->ZN[KP1]>=p->Q111_zs[qn] && p->ZN[KP1]<=p->Q111_ze[qn])
    {
        active_box_dummy(i,j,k) = true;
        ++cellCountBoxDummy;
    }

    // Topo
    size_t cellCountTopo = set_active_topo(p,a);

    // Make box to topo if inside topo?

    // guess particle demand
    if(p->Q24>0)
        ppcell = p->Q24;
    else
        ppcell = 0;
    // int gfminPPC=pgc->globalmin(floorf(minPPC));
    // if (gfminPPC<0) gfminPPC=0;
    // if(ppcell>gfminPPC)
    // {
    //     ppcell=gfminPPC;
    //     if(0==p->mpirank)
    //     cout<<"Reduced particles per cell to "<<ppcell<<" as min particles per cell is lower."<<endl;
    // }
    // double gQ41=pgc->globalmin(minPPC/ppcell);
    // if(ppcell!=0&&p->Q41*ppcell>gfminPPC)
    // {
    //     p->Q41=gQ41;
    //     if(0==p->mpirank)
    //     cout<<"Reduced packing factor to "<<p->Q41<<" to stay within max real particles per cell."<<endl;
    // }
    
    int partnum = (cellCountBox + cellCountTopo) * ppcell;
    maxparticle = ceil(p->Q25*double(pgc->globalisum(partnum)));
    PP2.reserve(ceil(p->Q25*double(cellCountBoxDummy*ppcell)));
}

/// @brief Calls seeding functions

/// Calls \p sediment_part::posseed_box and \p sediment_part::posseed_topo based on \p lexer::Q110 and \p lexer::Q101.
void sediment_part::seed(lexer *p, fdm *a)
{
    if(p->Q110>0)
        posseed_box(p,a);
	if(p->Q101>0)
        posseed_topo(p,a);
    if(p->Q111>0)
    {
        posseed_box_dummy(p,a);
        printDummyVTP(p,PP2);
        PP2.erase_all();
    }
}

/// @brief Seeds particle into boxes defined using `lexer::Q110`
void sediment_part::posseed_box(lexer *p, fdm *a)
{
    seed_srand(p);
	
    double x,y,z;
    int flag=1;
    size_t index;

    LOOP
        if(active_box(i,j,k)>0.0)
            for(int qn=0;qn<ppcell;++qn)
            {
                if(PP.size+1>0.9*PP.capacity)
                    PP.reserve();
                
                x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;

                index = PP.add(x,y,z,flag,a->u(i,j,k),a->v(i,j,k),a->w(i,j,k),p->Q41);
                pst->seeding(p, PP, index, ppcell);
            }
}

/// @brief Seeds particle into boxes defined using `lexer::Q111`
void sediment_part::posseed_box_dummy(lexer *p, fdm *a)
{
    seed_srand(p);
	
    double x,y,z;
    int flag=1;
    size_t index;

    LOOP
        if(active_box_dummy(i,j,k)>0.0)
            for(int qn=0;qn<ppcell;++qn)
            {
                if(PP2.size+1>0.9*PP2.capacity)
                    PP2.reserve();
                
                x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;

                index = PP2.add(x,y,z,flag,a->u(i,j,k),a->v(i,j,k),a->w(i,j,k),p->Q41);
            }
}

/// @brief Seeds particle in active topo cells
void sediment_part::posseed_topo(lexer *p, fdm *a)
{
    seed_srand(p);

    PLAINLOOP
    if(active_topo(i,j,k)>false)
    {
        seed_topo(p,a);
    }
}

/// @brief Seeds particle into suspension at inlet boundary
void sediment_part::posseed_suspended(lexer *p, fdm *a)
{
    seed_srand(p);
    
    double x,y,z;
    size_t index;
    for(int n=0;n<p->gcin_count;n++)
        if(p->gcin[n][3]>0)
        {
            i=p->gcin[n][0];
            j=p->gcin[n][1];
            k=p->gcin[n][2];
            if(a->topo(i,j,k)>=0.0)
            {
                if(PP.size+p->Q122>0.9*PP.capacity)
                    maxparticle=PP.reserve(PP.size+p->Q122);
                for(int qn=0;qn<p->Q122;++qn)
                {
                    x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                    y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                    z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
                    index=PP.add(x,y,z,1,a->u(i-1,j,k),a->v(i-1,j,k),a->w(i-1,j,k),p->Q41);
                    pst->seeding(p, PP, index, ppcell);
                }
            }
        }
}

/// @brief Seeds particle at points defined using `lexer::Q61`
void sediment_part::point_source(lexer *p, fdm *a)
{
    size_t index;
    for(size_t n=0;n<p->Q61;n++)
        if(p->count%p->Q61_i[n]==0)
        {
            index = PP.add(p->Q61_x[n],p->Q61_y[n],p->Q61_z[n],1,a->u(i,j,k),a->v(i,j,k),a->w(i,j,k),p->Q41);
            pst->seeding(p, PP, index, ppcell);
        }
}

/// @brief Seeds particles into active topo cells at inlet boundary
void sediment_part::topo_influx(lexer *p, fdm *a)
{
    seed_srand(p);
    for(int n=0;n<p->gcin_count;n++)
    {
        if(p->gcin[n][3]>0)
        {
            i=p->gcin[n][0];
            j=p->gcin[n][1];
            k=p->gcin[n][2];
            if(active_topo(i,j,k)>0.0)
            {
                seed_topo(p,a);
            }
        }
    }
}
/// @brief Seeds `rand()` function
void sediment_part::seed_srand(lexer *p)
{
    if(p->Q29>0)
        srand(p->Q29);

    if(p->Q29==0)
        srand((unsigned)time(0)*(p->mpirank+1));
}

/// @brief Seeds particle in relation to `fdm::topo`
void sediment_part::seed_topo(lexer *p, fdm *a)
{
    double tolerance = 5e-18;
    double x,y,z,ipolTopo,ipolSolid;
    int flag=0;
    size_t index;

    if(PP.size+ppcell>0.9*PP.capacity)
        PP.reserve();

    for(int qn=0;qn<ppcell*1000;++qn)
    {   
        x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
        y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
        k = 0;
        z = p->ZN[KP]+0.5*p->DZP[KP]-a->topo(i,j,k) - p->Q102*p->DZN[KP]*double(rand() % irand)/drand;
        k = p->posc_k(z);

        ipolTopo = p->ccipol4_b(a->topo,x,y,z);
        ipolSolid = p->ccipol4_b(a->solid,x,y,z);

        if (!(ipolTopo>tolerance||ipolTopo<-p->Q102*p->DZN[KP]||ipolSolid<0))
        {
            index=PP.add(x,y,z,flag,0,0,0,p->Q41);
            auto result = pst->seeding(p, PP, index, p->Q41*ppcell);
            if(result==seedReturn::STOP)
                break;
            if(result==seedReturn::REMOVE)
            {
                PP.erase(index);
                break;
            }
        }
    }
}

void sediment_part::solid_influx(lexer *p, fdm *a)
{
    seed_srand(p);
    PLAINLOOP
    if(active_topo(i,j,k)==10.0)
    {
        seed_topo(p,a);
    }
}

size_t sediment_part::set_active_topo(lexer *p, fdm *a)
{
    const double tolerance=5e-10;
    size_t cellCountTopo=0;
    BASELOOP
    {
        active_topo(i,j,k) = false;
        if( (a->topo(i,j,k)<0.5*p->DZN[KP]-tolerance) && (a->topo(i,j,k)>-p->DZN[KP]*ceil(p->Q102)-tolerance) && (a->solid(i,j,k)>=-p->DXM))
        {
            active_topo(i,j,k) = true;
            if(p->flag1[Im1JK]==SOLID_FLAG&&p->flag1[IJK]==WATER_FLAG)
            active_topo(i,j,k) = 10.0;
            cellCountTopo++;
        }
    }
    return cellCountTopo;
}

void sediment_part::seedDummy(lexer *p, fdm *a, particles_obj &PP)
{
    if(p->origin_i==0)
    {
        PP.reserve(PP.size+ppcell*p->knoy*p->knoz);
        i=0;
        JLOOP
        for(k=1; k<p->knoz; ++k)
            seedDummyCell(p,a,PP);
    }
    if(p->origin_i+p->knox == p->gknox)
    {
        PP.reserve(PP.size+ppcell*p->knoy*p->knoz);
        i = p->gknox;
        JLOOP
        for(k=1; k<p->knoz; ++k)
            seedDummyCell(p,a,PP);
    }
    if(p->origin_j==0)
    {
        PP.reserve(PP.size+ppcell*p->knox*p->knoz);
        j=0;
        ILOOP
        for(k=1; k<p->knoz; ++k)
            seedDummyCell(p,a,PP);
    }
    if(p->origin_j+p->knoy == p->gknoy)
    {
        PP.reserve(PP.size+ppcell*p->knox*p->knoz);
        j = p->gknoy;
        ILOOP
        for(k=1; k<p->knoz; ++k)
            seedDummyCell(p,a,PP);
    }
    if(p->origin_k==0)
    {
        PP.reserve(PP.size+ppcell*p->knoy*p->knox);
        k=0;
        ILOOP
        JLOOP
        seedDummyCell(p,a,PP);
    }
    cout<<"Dummy seeding done on "<<p->mpirank<<".\nSeeded "<<PP.size<<" particles."<<endl;
}

void sediment_part::seedDummyCell(lexer *p, fdm *a, particles_obj &PP)
{
    double tolerance = 5e-18;
    double x,y,z,ipolSolid;

    for(int qn=0;qn<ppcell;++qn)
    {   
        x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
        y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
        z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
        if(a->topo(i,j,k) +0.5*p->Q102*p->DZN[KP]>=0)
            break;
        ipolSolid = p->ccipol4_b(a->solid,x,y,z);
        if(!(ipolSolid<0))
            PP.add(x,y,z,0);
    }
}