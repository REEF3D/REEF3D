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

// Unit test of the open fraction of a cell face cut by the level set of a solid (D 1 2,
// src/solid_aperture.h): a plane level set cuts a line and a rectangle exactly.
//
//   g++ -std=c++17 -I ../../src solid_aperture_test.cpp -o solid_aperture_test && ./solid_aperture_test

#include"solid_aperture.h"
#include<cmath>
#include<cstdio>
#include<cstdlib>

static int nfail = 0;

static void check(const char *what, double val, double ref, double tol)
{
    if(!(std::fabs(val-ref)<=tol))
    {
    std::printf("FAIL %s: %.15g, expected %.15g\n",what,val,ref);
    ++nfail;
    }
}

// the part of the unit square where a + b x + c y > 0, by sampling
static double sampled(double a, double b, double c)
{
    const int n = 2000;
    long cnt = 0;
    for(int p=0; p<n; ++p)
    for(int q=0; q<n; ++q)
    if(a + b*(p+0.5)/n + c*(q+0.5)/n > 0.0)
    ++cnt;
    return double(cnt)/double(n)/double(n);
}

int main()
{
    // line
    check("line open",solid_cut_line(1.0,2.0),1.0,0.0);
    check("line closed",solid_cut_line(-1.0,-2.0),0.0,0.0);
    check("line on the surface",solid_cut_line(0.0,0.0),0.0,0.0);
    check("line 1/4",solid_cut_line(1.0,-3.0),0.25,1.0e-15);
    check("line 3/4",solid_cut_line(-1.0,3.0),0.75,1.0e-15);

    // rectangle
    check("quad open",solid_cut_quad(1.0,2.0,3.0,4.0),1.0,0.0);
    check("quad closed",solid_cut_quad(-1.0,-2.0,-3.0,0.0),0.0,0.0);
    check("quad half",solid_cut_quad(1.0,1.0,-1.0,-1.0),0.5,1.0e-15);
    check("quad 30 %",solid_cut_quad(0.3,0.3,-0.7,-0.7),0.3,1.0e-14);
    check("quad diagonal",solid_cut_quad(1.0,0.0,-1.0,0.0),0.5,1.0e-14);
    check("quad corner 1/8",solid_cut_quad(0.5,-0.5,-1.5,-0.5),0.125,1.0e-14);
    check("quad corner 7/8",solid_cut_quad(-0.5,0.5,1.5,0.5),0.875,1.0e-14);

    // planes at any angle against the sampled area (corners in order around the unit square)
    std::srand(7);
    for(int n=0; n<200; ++n)
    {
        const double a = 2.0*std::rand()/RAND_MAX - 1.0;
        const double b = 2.0*std::rand()/RAND_MAX - 1.0;
        const double c = 2.0*std::rand()/RAND_MAX - 1.0;
        check("quad plane",solid_cut_quad(a,a+b,a+b+c,a+c),sampled(a,b,c),2.0e-3);
        // symmetric under the sign of the level set
        check("quad complement",solid_cut_quad(a,a+b,a+b+c,a+c)+solid_cut_quad(-a,-a-b,-a-b-c,-a-c),1.0,1.0e-13);
    }

    if(nfail==0)
    std::printf("solid_aperture_test: all passed\n");

    return nfail==0 ? 0 : 1;
}
