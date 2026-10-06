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

#include"seastate_f.h"
#include"seastate_amr.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<fstream>
#include<iomanip>
#include<vector>
#include<sys/stat.h>
#include<sys/types.h>

/*--------------------------------------------------------------------
Handover of 2D spectra to the phase-resolved wave generation (A 760)

For every point A 760 x y the spectrum of the cell that contains the
point is written to

  REEF3D_SEASTATE_Spectra/spectrum-file-2d_P<n>.dat

in the format of REEF3D's 2D spectrum input (B 85 11, file
spectrum-file-2d.dat, wave_lib_spectrum_file.cpp):

  fq_dir  theta_1 ... theta_M              [rad]
  omega_1 S(omega_1,theta_1) ... S(omega_1,theta_M)
  ...                                      [rad/s, m^2 s/rad^2]

theta is the direction of propagation, counter-clockwise from +x (the
convention of SEASTATE and of the wave generation, cos(beta) x +
sin(beta) y). The directions are the sector of half-width A 761 around
the grid direction closest to the mean direction, written unwrapped and
increasing (they may be < 0 or > 2 pi). The wave generation builds one
component per (omega, theta) with amplitude sqrt(2 S domega dbeta),
domega the forward difference of the omegas (the last one repeated),
dbeta the mean direction step. S is therefore written as
E(sig,theta) dsig/domega, so that every component carries exactly the
energy of its SEASTATE bin: the Hs of the generated waves equals the
Hs of the sector.

A companion file seastate-wavegen_P<n>.txt lists the location, the
cell, the water depth, Hs of the full spectrum and of the sector, Tp,
the mean direction, and the ctrl keys for the wave generation.

The files are rewritten after every step (the latest spectrum).
--------------------------------------------------------------------*/

void seastate_f::handover(lexer *p, ghostcell *pgc)
{
    const seastate_grid &g = *e->grid;
    const double pi = 3.14159265358979323846;

    if(p->mpirank==0)
    mkdir("./REEF3D_SEASTATE_Spectra",0777);

    pgc->globalmax(0.0);    // the folder exists before any rank writes

    for(int n=0; n<p->A760; ++n)
    {
    const double xp = p->A760_x[n], yp = p->A760_y[n];

    // the rank that owns the cell containing the point writes the file
    int owner = 0, ci=-1, cj=-1;

        ILOOP
        JLOOP
        if(p->XN[IP]<=xp && xp<p->XN[IP1] && p->YN[JP]<=yp && yp<p->YN[JP1])
        {
        owner = 1;
        ci = i;
        cj = j;
        }

    const int found = int(pgc->globalisum(owner));

        if(found==0)
        {
            if(p->mpirank==0)
            cout<<"SEASTATE handover point "<<n+1<<" ("<<xp<<", "<<yp<<") is outside the domain"<<endl;

            continue;
        }

        if(owner==0)
        continue;

    // mesh refinement: the finest grid that holds the point
    lexer *q = p;
    fdm_seastate *ee = e;

    if(pamr!=nullptr)
    pamr->locate(xp,yp,q,ee,ci,cj);

    i = ci;
    j = cj;

    char name[256];
    sprintf(name,"./REEF3D_SEASTATE_Spectra/spectrum-file-2d_P%i.dat",n+1);
    ofstream out(name);

    char iname[256];
    sprintf(iname,"./REEF3D_SEASTATE_Spectra/seastate-wavegen_P%i.txt",n+1);
    ofstream info(iname);

        if(ee->wet(i,j)==0 || ee->N->spec(i,j)==nullptr)
        {
        info<<"REEF3D::SEASTATE handover point "<<n+1<<" ("<<xp<<", "<<yp<<"): the cell is dry or land, no spectrum"<<endl;
        continue;
        }

    const float *N = ee->N->spec(i,j);

    seastate_param sp;
    sp.compute(g,N);

    // direction sector around the grid direction closest to the mean direction
    const int m0 = int(std::lround(sp.dir*pi/180.0/g.dtheta)) % g.ndir;
    int nh = int(std::floor(p->A761*pi/180.0/g.dtheta + 1.0e-9));
    nh = std::max(nh,1);

    int ma, mb;
        if(2*nh+1>=g.ndir)
        {
        ma = m0 - g.ndir/2;
        mb = ma + g.ndir - 1;
        }
        else
        {
        ma = m0 - nh;
        mb = m0 + nh;
        }

    // domega of the wave generation: forward difference, the last one repeated
    std::vector<double> dw(g.nsig);
    for(int l=0; l<g.nsig-1; ++l)
    dw[l] = g.sig[l+1]-g.sig[l];
    dw[g.nsig-1] = dw[g.nsig-2];

    out<<"fq_dir";
    out<<setprecision(10);
    for(int m=ma; m<=mb; ++m)
    out<<" "<<double(m)*g.dtheta;
    out<<endl;

    double m0sec = 0.0;

        for(int l=0; l<g.nsig; ++l)
        {
        out<<g.sig[l];

            for(int m=ma; m<=mb; ++m)
            {
            const int mm = ((m%g.ndir)+g.ndir)%g.ndir;
            const double E = g.sig[l]*double(N[g.bin(l,mm)]);

            out<<" "<<E*g.dsig[l]/dw[l];
            m0sec += E*g.dsig[l]*g.dtheta;
            }

        out<<endl;
        }

    out.close();

    const double hs_sec = 4.0*std::sqrt(std::max(m0sec,0.0));

    info<<"REEF3D::SEASTATE handover point "<<n+1<<endl;
    info<<"location "<<xp<<" "<<yp<<", cell centre "<<q->XP[IP]<<" "<<q->YP[JP]<<", cell size "<<q->DXN[IP]<<" m, simtime "<<p->simtime<<endl;
    info<<"water depth "<<ee->depth(i,j)<<" m"<<endl;
    info<<"Hs "<<sp.Hs<<" m (full spectrum), "<<hs_sec<<" m (written sector "<<double(ma)*g.dtheta*180.0/pi<<" to "<<double(mb)*g.dtheta*180.0/pi<<" deg)"<<endl;
    info<<"Tp "<<sp.Tp<<" s, Tm01 "<<sp.Tm01<<" s, mean direction "<<sp.dir<<" deg, spread "<<sp.spread<<" deg"<<endl;
    info<<endl;
    info<<"wave generation (FNPF/NHFLOW): copy spectrum-file-2d_P"<<n+1<<".dat to spectrum-file-2d.dat and use"<<endl;
    info<<"B 85 11"<<endl;
    info<<"B 92 31          (linear; 32 or 33 for second order, with B 89 0)"<<endl;
    info<<"B 93 "<<sp.Hs<<" "<<sp.Tp<<endl;
    info<<"B 94 "<<ee->depth(i,j)<<endl;
    info<<"B 130 1          (needed with B 85 11; the spreading comes from the file)"<<endl;
    info<<"directions are absolute (propagation, ccw from +x): rotate the SEASTATE frame into the phase-resolved frame if needed"<<endl;
    info.close();
    }
}
