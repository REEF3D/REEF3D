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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"spectral_store.h"
#include<algorithm>

spectral_store::spectral_store(int imin_, int jmin_, int ni_, int nj_, int nbin_, int tile_)
    : imin(imin_), jmin(jmin_), ni(ni_), nj(nj_), nbin(nbin_), tile(std::max(tile_,1)),
      ntile_alloc(0), ncell_active(0), ncell_alloc(0)
{
    ntx = (ni + tile - 1)/tile;
    nty = (nj + tile - 1)/tile;

    toff.assign(size_t(ntx)*nty,-1);
    tbase.assign(size_t(ntx)*nty,-1);
    act.assign(size_t(ni)*nj,0);
}

void spectral_store::build(const int *mask)
{
    std::fill(toff.begin(),toff.end(),-1);
    std::fill(tbase.begin(),tbase.end(),-1);
    ntile_alloc  = 0;
    ncell_active = 0;
    ncell_alloc  = 0;

    for(int ii=0; ii<ni; ++ii)
    for(int jj=0; jj<nj; ++jj)
    {
    const size_t n = size_t(ii)*nj + jj;
    act[n] = mask[n]>0 ? 1 : 0;

        if(act[n])
        {
        ++ncell_active;
        int &t = toff[size_t(ii/tile)*nty + jj/tile];

        if(t<0)
        t = ntile_alloc++;
        }
    }

    // first cell of every allocated tile; tiles are clipped at the end of the range
    for(int ti=0; ti<ntx; ++ti)
    for(int tj=0; tj<nty; ++tj)
    {
    const size_t t = size_t(ti)*nty + tj;

        if(toff[t]>=0)
        {
        tbase[t] = ncell_alloc;
        ncell_alloc += long(std::min(tile,ni-ti*tile))*std::min(tile,nj-tj*tile);
        }
    }

    // release the old pool before allocating the new one (no peak of two pools)
    std::vector<float>().swap(pool);
    pool.assign(size_t(ncell_alloc)*nbin,0.0f);
}

void spectral_store::fill(float val)
{
    for(int ii=0; ii<ni; ++ii)
    for(int jj=0; jj<nj; ++jj)
    {
    float *s = spec(ii+imin,jj+jmin);

    if(s!=nullptr)
    std::fill(s,s+nbin,act[size_t(ii)*nj+jj] ? val : 0.0f);
    }
}

bool spectral_store::copy_from(const spectral_store &s)
{
    if(s.ni!=ni || s.nj!=nj || s.nbin!=nbin || s.tile!=tile || s.tbase!=tbase || s.pool.size()!=pool.size())
    return false;

    std::copy(s.pool.begin(),s.pool.end(),pool.begin());
    return true;
}
