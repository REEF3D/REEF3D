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

#ifndef SPECTRAL_STORE_H_
#define SPECTRAL_STORE_H_

#include<cstddef>
#include<vector>

/*--------------------------------------------------------------------
REEF3D::Spectral - block-sparse storage of the wave action density
N(x,y,sigma,theta), single precision.

The 2D index range of a rank, ghost cells included (imin..imin+ni-1,
jmin..jmin+nj-1, the same range and cell order as the IJ macro), is cut
into square tiles of tile x tile cells. Storage is allocated only for
tiles that contain at least one active (sea) cell, so land and dry
areas cost no spectral memory. Inside a tile each cell holds one
contiguous spectrum of nbin floats (layout: see spectral_grid::bin).

  build(mask)   mask[(i-imin)*nj + (j-jmin)] > 0 marks an active cell
  spec(i,j)     pointer to the spectrum of cell (i,j), nullptr if the
                tile is not allocated
  active(i,j)   cell is active (allocated tiles also hold inactive
                cells, which are kept at zero)

No lexer or MPI dependency, so the class is unit-testable stand-alone.
--------------------------------------------------------------------*/

class spectral_store
{
public:
    spectral_store(int imin, int jmin, int ni, int nj, int nbin, int tile);

    void build(const int *mask);
    void fill(float val);

    float *spec(int i, int j)
    {
        const long o = offset(i,j);
        return o<0 ? nullptr : &pool[size_t(o)];
    }

    const float *spec(int i, int j) const
    {
        const long o = offset(i,j);
        return o<0 ? nullptr : &pool[size_t(o)];
    }

    bool active(int i, int j) const
    {
        if(i<imin || j<jmin || i>=imin+ni || j>=jmin+nj)
        return false;
        return act[size_t(i-imin)*nj + (j-jmin)]!=0;
    }

    int nbins() const {return nbin;}
    int tile_size() const {return tile;}
    int tiles_total() const {return ntx*nty;}
    int tiles_allocated() const {return ntile_alloc;}
    long cells_allocated() const {return long(ntile_alloc)*tile*tile;}
    long cells_active() const {return ncell_active;}
    size_t bytes() const {return pool.size()*sizeof(float) + toff.size()*sizeof(int) + act.size();}
    size_t bytes_dense() const {return size_t(ni)*nj*nbin*sizeof(float);}

private:
    long offset(int i, int j) const
    {
        if(i<imin || j<jmin || i>=imin+ni || j>=jmin+nj)
        return -1;

        const int ii = i-imin;
        const int jj = j-jmin;
        const int t  = toff[size_t(ii/tile)*nty + jj/tile];

        if(t<0)
        return -1;

        return (long(t)*tile*tile + long(ii%tile)*tile + (jj%tile))*nbin;
    }

    int imin, jmin, ni, nj, nbin, tile;
    int ntx, nty, ntile_alloc;
    long ncell_active;

    std::vector<int> toff;      // tile -> allocated tile number, -1: not allocated
    std::vector<char> act;      // active flag per cell
    std::vector<float> pool;    // allocated tiles x tile^2 cells x nbin
};

#endif
