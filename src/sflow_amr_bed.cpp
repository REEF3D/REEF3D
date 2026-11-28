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

#include"sflow_amr.h"
#include"lexer.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iostream>
#include<limits>
#include<string>

/*--------------------------------------------------------------------
Fine bed of the refined levels (G 60 1).

Without it a patch cell takes the bed of level 0, prolonged limited-linearly level by level: a
harbour, a pier or a street the level-0 cells do not resolve does not exist on the patches either.
With G 60 1 the file amr-bed.dat gives one or more tiles of bed samples on regular grids in the model
coordinates. A cell of level l >= 1 whose parent has all its children inside a tile (spacing not larger
than the children) takes

    bed = bed of the parent + (sample mean of the cell - mean of the sample means of the 2 x 2 children)

from the finest such tile: the detail below the parent comes from the tile, and the children average to
the parent as with the prolonged bed (the composite system stays well balanced: a lake at rest stays at
rest). Cells whose siblings are not all inside a tile keep the prolonged bed.

amr-bed.dat (text; '#' lines are comments):

    tile nx ny x0 y0 dx
    z(0,0) z(1,0) ... z(nx-1,0) z(0,1) ...        nx*ny values, x fastest, the same z as geo.dat

  sample (i,j) at x0 + (i+0.5) dx, y0 + (j+0.5) dy (dy = dx); nan: no data (the cell falls back to
  the next tile or to the prolongation). Several tiles in one file, in any order.

The mean over a parent therefore always comes from the coarser level, down to geo.dat on level 0: build
geo.dat from the same data (the cell means of the samples) where the tiles lie, so that level 0 holds the
same mean bed as the tile.

Cartesian grids only (not G 50, not multi-block); not with the sediment transport (S 10 1), whose bed
evolves.
--------------------------------------------------------------------*/

void sflow_amr::fbed_read(lexer *p)
{
    fbed.clear();

    if(p->geo_curv==1 || nblk>0 || p->S10!=0)
    {
        if(p->mpirank==0)
        cout<<"!!! SFLOW AMR fine bed (G 60 1): only on a Cartesian grid (not G 50, not multi-block) without S 10 -- fine bed off !!!"<<endl;
        return;
    }

    ifstream in("amr-bed.dat");

    if(!in.good())
    {
        if(p->mpirank==0)
        cout<<"!!! SFLOW AMR fine bed (G 60 1): amr-bed.dat cannot be opened !!!"<<endl;
        exit(1);
    }

    auto fail = [&](const string &what)
    {
        if(p->mpirank==0)
        cout<<"!!! SFLOW AMR fine bed (G 60 1): amr-bed.dat, tile "<<fbed.size()+1<<": "<<what<<" !!!"<<endl;
        exit(1);
    };

    string tok;
    while(in>>tok)
    {
        if(tok[0]=='#')
        {
            getline(in,tok);
            continue;
        }

        if(tok!="tile")
        fail("expected 'tile nx ny x0 y0 dx', found '"+tok+"'");

        fbed_tile t;
        if(!(in>>t.nx>>t.ny>>t.x0>>t.y0>>t.dx) || t.nx<1 || t.ny<1 || !(t.dx>0.0))
        fail("bad header");

        t.z.resize((size_t)t.nx*t.ny);

        for(size_t n=0; n<t.z.size(); ++n)
        {
            if(!(in>>tok))
            fail("fewer than nx*ny values");

            if(tok[0]=='#')
            {
                getline(in,tok);
                --n;
                continue;
            }

            char *e = nullptr;
            const double v = strtod(tok.c_str(),&e);

            if(e==tok.c_str())
            fail("not a number: '"+tok+"'");

            t.z[n] = float(v);
        }

        fbed.push_back(std::move(t));
    }

    // finest tile first
    stable_sort(fbed.begin(),fbed.end(),[](const fbed_tile &a, const fbed_tile &b){ return a.dx<b.dx; });

    if(p->mpirank==0)
    {
        cout<<"SFLOW AMR fine bed (G 60 1): "<<fbed.size()<<" tile(s) from amr-bed.dat"<<endl;

        for(const fbed_tile &t : fbed)
        cout<<"   "<<t.nx<<" x "<<t.ny<<" samples of "<<t.dx<<" m from ("<<t.x0<<", "<<t.y0<<")"<<endl;
    }
}

// mean of the samples of the finest tile that covers cell (I,J) of level l and resolves it (spacing <= cell);
// false where no tile does or a sample is missing (nan)
bool sflow_amr::fbed_mean(int l, int I, int J, double &z)
{
    const double xa = gnode(0,l,I), xb = gnode(0,l,I+1);
    double ya = gnode(1,l,J), yb = gnode(1,l,J+1);

    // one cell in y: the cell stands for every row
    if(p0->j_dir==0)
    {
        ya = -numeric_limits<double>::max();
        yb =  numeric_limits<double>::max();
    }

    const double w = MIN(xb-xa, p0->j_dir==0 ? xb-xa : yb-ya);

    for(const fbed_tile &t : fbed)
    {
        const double tol = 1.0e-6*t.dx;

        if(t.dx>w*(1.0+1.0e-9))
        continue;

        // the cell inside the tile
        if(xa<t.x0-tol || xb>t.x0+t.nx*t.dx+tol)
        continue;

        if(p0->j_dir==1 && (ya<t.y0-tol || yb>t.y0+t.ny*t.dx+tol))
        continue;

        // samples with their centre in [xa,xb) x [ya,yb)
        const int i0 = MAX(0, (int)ceil((xa-t.x0)/t.dx - 0.5 - 1.0e-9));
        const int i1 = MIN(t.nx, (int)ceil((xb-t.x0)/t.dx - 0.5 - 1.0e-9));
        const int j0 = p0->j_dir==0 ? 0 : MAX(0, (int)ceil((ya-t.y0)/t.dx - 0.5 - 1.0e-9));
        const int j1 = p0->j_dir==0 ? t.ny : MIN(t.ny, (int)ceil((yb-t.y0)/t.dx - 0.5 - 1.0e-9));

        if(i1<=i0 || j1<=j0)
        continue;

        double s = 0.0;
        bool ok = true;

        for(int j=j0; j<j1 && ok; ++j)
        for(int i=i0; i<i1; ++i)
        {
            const float v = t.z[(size_t)j*t.nx + i];

            if(!std::isfinite(v))
            {
                ok = false;
                break;
            }

            s += v;
        }

        if(!ok)
        continue;

        z = s/double((i1-i0)*(j1-j0));
        return true;
    }

    return false;
}
