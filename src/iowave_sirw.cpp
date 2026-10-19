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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"

// SIRW absorption set-up (B 98 8, B 99 8): gauges and band from B 534 / B 535, or from the reference period
// (B 528, else the peak period of the B 92 wave): flow = 0.4 fp, fhigh = 3 fp, dx = 0.4 L(fhigh) (k dx <= 0.8 pi),
// d = 0.7 L(flow) (Frigaard & Brorsen 1995). What the prediction does not explain is absorbed with the shallow-water
// velocity: at x+ the deviation of the boundary surface from the predicted outgoing wave; at x- (FNPF) the deviation
// from the target minus the predicted reflected wave after a first-order low pass with the corner frequency flow/2,
// because the full deviation from the linear target also contains the near field and bound waves of the generated
// waves (with it the reflection of the paddle was 6 % instead of 1.5 % in the regular-wave test). A higher-order low
// pass lags by more than 90 degrees below the band and drove the seiche modes of the tank.
// NHFLOW: the predicted wave (surface, velocity per level, non-hydrostatic pressure) is part of the outside state of
// an absorbing Riemann edge (x+, and x- with a wave theory, nhflow_open_edges); with a wavemaker at x- the active
// generation of B 98 4 gets the velocity correction and the pressure of the predicted reflected wave.
void iowave::sirw_ini(lexer *p)
{
    if(!sirw_in_req && !sirw_out_req)
    return;
    
    const double Tr = p->B528>0.0 ? p->B528 : (p->wTp>0.0 ? p->wTp : p->wT);
    const double h  = p->wd;
    const double g  = 9.81;
    
    auto Lof = [&](double f)
    {
        const double w = 2.0*PI*f;
        double k = std::max(w/sqrt(g*h), w*w/g);
        for(int q=0; q<60; ++q)
        {
        const double th = tanh(k*h);
        k -= (g*k*th - w*w)/(g*th + g*k*h*(1.0-th*th));
        }
        return 2.0*PI/k;
    };
    
    for(int side=-1; side<=1; side+=2)
    {
        if((side<0 && !sirw_in_req) || (side>0 && !sirw_out_req))
        continue;
        
        double d  = side<0 ? p->B534_d  : p->B535_d;
        double dx = side<0 ? p->B534_dx : p->B535_dx;
        double fl = side<0 ? p->B534_fl : p->B535_fl;
        double fh = side<0 ? p->B534_fh : p->B535_fh;
        
        if((fl<=0.0 || fh<=0.0) && Tr<=0.0)
        {
            if(p->mpirank==0)
            cout<<"iowave: B "<<(side<0?"98":"99")<<" 8 needs the band (B "<<(side<0?534:535)<<") or the reference period B 528; switched off"<<endl;
            continue;
        }
        
        if(fl<=0.0) fl = 0.4/Tr;
        if(fh<=0.0) fh = 3.0/Tr;
        if(dx<=0.0) dx = 0.4*Lof(fh);
        if(d<=0.0)  d  = 0.7*Lof(fl);
        
        const double xb = side<0 ? p->global_xmin : p->global_xmax;
        const double xg = xb - side*(d+dx);
        
        if(xg<=p->global_xmin || xg>=p->global_xmax)
        {
            if(p->mpirank==0)
            cout<<"iowave: B "<<(side<0?"98":"99")<<" 8: gauges at "<<xb-side*d<<" / "<<xg<<" m outside the domain; switched off"<<endl;
            continue;
        }
        
        if(side<0)
        {
        sirw_in.ini(p,-1,xb,h,d,dx,fl,fh);
        sw_tau = 1.0/(PI*fl);
        }
        else
        {
        sirw_out.ini(p,1,xb,h,d,dx,fl,fh);
        }
    }
}

// long band at x- (FNPF B 98 8): first-order low pass of the remaining deviation x of column n
double iowave::sw_lowpass(lexer *p, int n, double x)
{
    return ab2_lowpass(p,sw_l1,sw_c1,n,x,sw_tau);
}
