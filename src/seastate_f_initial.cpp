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
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"wave_lib_spectrum.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<string>
#include<vector>

/*--------------------------------------------------------------------
Parametric initial spectrum, A 710 1, from the wave keys of REEF3D's
irregular wave generation:

  B 85   frequency spectrum S(omega): 1 PM, 2 JONSWAP, 3 Torsethaugen
  B 93   Hs, Tp
  B 88   JONSWAP gamma
  B 130  directional spreading: 0 unidirectional, 1 cos^s (PNJ),
         2 Mitsuyasu cos^2s(theta/2)
  B 131  main direction [deg], direction of propagation, ccw from +x
  B 134  spreading shape parameter s (B 135 for Mitsuyasu with B 134 <= 0)

E(sig,theta) = S(sig) D(sig,theta), N = E/sig. D is normalised on the
discrete directional grid for every frequency, so the directional
spreading conserves the energy of S exactly.

  parametric_spectrum   the spectrum of one cell, also used as the
                        boundary spectrum (A 711 1)
  initial_parametric    the same spectrum in every active cell (A 710 1)
--------------------------------------------------------------------*/

void seastate_f::initial_parametric(lexer *p, ghostcell *pgc)
{
    std::vector<float> N;
    parametric_spectrum(p,pgc,N,"A 710 1");

    IMALOOP
    JMALOOP
    if(e->wet(i,j)==1)
    {
    float *s = e->N->spec(i,j);
    std::copy(N.begin(),N.end(),s);
    }
}

void seastate_f::parametric_spectrum(lexer *p, ghostcell *pgc, std::vector<float> &N, const char *key)
{
    const double pi = 3.14159265358979323846;
    std::string msg;

    if(!(p->B93_1>0.0) || !(p->B93_2>0.0))
    msg = std::string(key) + " needs B 93 (Hs, Tp)";
    else if(p->B85!=1 && p->B85!=2 && p->B85!=3)
    msg = std::string(key) + ": B 85 must be 1 (PM), 2 (JONSWAP) or 3 (Torsethaugen)";
    else if(p->B130!=0 && p->B130!=1 && p->B130!=2)
    msg = std::string(key) + ": B 130 must be 0, 1 or 2";

    if(!msg.empty())
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  "<<msg<<endl<<endl;

        pgc->final(true);
    }

    const seastate_grid &g = *e->grid;

    p->wHs = p->B93_1;
    p->wTp = p->B93_2;
    p->wwp = 2.0*pi/p->wTp;

    wave_lib_spectrum spec;

    const double main = p->B131*pi/180.0;

    // nearest direction bin of the main direction
    double mw = std::fmod(main,2.0*pi);
    if(mw<0.0)
    mw += 2.0*pi;
    const int mmain = int(std::lround(mw/g.dtheta))%g.ndir;

    N.assign(g.nbin,0.0f);
    std::vector<double> D(g.ndir);

    for(int l=0; l<g.nsig; ++l)
    {
    const double S = spec.wave_spectrum(p,g.sig[l]);
    double sum = 0.0;

        if(p->B130>0)
        {
        // spreading_function works with B 131 in radians
        const double b131 = p->B131;
        p->B131 = main;

            for(int m=0; m<g.ndir; ++m)
            {
            double beta = g.theta[m];

            while(beta-main>pi)
            beta -= 2.0*pi;

            while(beta-main<-pi)
            beta += 2.0*pi;

            D[m] = spec.spreading_function(p,beta,g.sig[l]);

            if(!(D[m]>0.0))
            D[m] = 0.0;

            sum += D[m]*g.dtheta;
            }

        p->B131 = b131;
        }

        // unidirectional, or spreading narrower than one bin
        if(!(sum>0.0))
        {
        std::fill(D.begin(),D.end(),0.0);
        D[mmain] = 1.0;
        sum = g.dtheta;
        }

        for(int m=0; m<g.ndir; ++m)
        N[g.bin(l,m)] = float(S*D[m]/sum/g.sig[l]);
    }

    seastate_param sp;
    sp.compute(g,N.data());

    if(p->mpirank==0)
    {
    cout<<"SEASTATE parametric spectrum ("<<key<<"): B 85 "<<p->B85<<", Hs "<<p->wHs<<" m, Tp "<<p->wTp<<" s, direction "<<p->B131<<" deg, B 130 "<<p->B130<<endl;
    cout<<"  discrete: Hs "<<setprecision(5)<<sp.Hs<<" m, Tp "<<sp.Tp<<" s, Tm01 "<<sp.Tm01<<" s, direction "<<sp.dir<<" deg, spread "<<sp.spread<<" deg"<<endl<<endl;
    cout<<setprecision(6);
    }
}
