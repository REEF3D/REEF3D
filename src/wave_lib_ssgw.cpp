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
Author: Csaba Pakozdi
--------------------------------------------------------------------*/

#include"wave_lib_ssgw.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_ssgw::wave_lib_ssgw(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{
    N   = p->B170;  // default = 1024
    tol = 1e-14;
    
    if (p->B91 == 1)
    {
        // Wave length known
        allocated = false;
        allocated = resizing();
        setWave(wk, wdt, wH);
        getPhysicsParameters();
        surfaceCalculated = computeSurfaceVariables();
    } 
    else if (p->B93 == 1)
    {
        // Wave period known
        wk = 2.0*PI/(wT*sqrt(9.81*wdt));
        double delta_omega = 100000;
        int count = 0;
            
        allocated = false;
        allocated = resizing();
        setWave(wk, wdt, wH);
        getPhysicsParameters();
        surfaceCalculated = computeSurfaceVariables();
        
        while (fabs(delta_omega/ww) > 1e-7 && count < 30)
        {
            delta_omega = ww - ParameterValue.phaseVelocity*wk;
            wk += delta_omega/ParameterValue.groupVelocity;
         
            if (wk <= 0)
            {
                if (p->mpirank == 0)
                {
                    cout<<"SSGW has to use shallow water assumption to calculate wavenumber!"<<endl;
                }
                wk = ww/sqrt(9.81*wdt);
                break;
            }
            else
            {
                setWave(wk, wdt, wH);
                getPhysicsParameters();
                surfaceCalculated = computeSurfaceVariables();
            }

            count++;
        }

        wL = 2.0*PI/wk;
        p->wL = wL;
        p->wk = wk;
    }


    // Soliton
/*    xs.resize(8593);
    vector<double> xs_swapped = xs;
    std::reverse(xs_swapped.begin(), xs_swapped.end()); 
    xs_swapped.resize(8593-7792);
    std::reverse(xs_swapped.begin(), xs_swapped.end()); 
    xs = xs_swapped;
    ys.resize(8593);
    vector<double> ys_swapped = ys;
    std::reverse(ys_swapped.begin(), ys_swapped.end()); 
    ys_swapped.resize(8593-7792);
    std::reverse(ys_swapped.begin(), ys_swapped.end()); 
    ys = ys_swapped;
    phis.resize(8593);
    vector<double> phis_swapped = phis;
    std::reverse(phis_swapped.begin(), phis_swapped.end()); 
    phis_swapped.resize(8593-7792);
    std::reverse(phis_swapped.begin(), phis_swapped.end()); 
    phis = phis_swapped;
    wL = xs.back()-xs.front()+0.1;
*/

    if(p->mpirank==0)
    {
        cout<<"Wave_Lib: steady surface gravity waves; ";
        cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wdt: "<<wdt<<" kd: "<<wdt*wk<<endl;
        writeResult(".");
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_ssgw::~wave_lib_ssgw()
{
}

double wave_lib_ssgw::wave_eta(lexer *p, double x, double y)
{
    // Transform x location to current position xcurr at time instance p->simtime
    xcurr = fabs(modulo(x - ParameterValue.phaseVelocity*p->simtime, wL));
    xcurr = xcurr>0.5*wL ? xcurr-wL : xcurr;
   
    // Linear interpolation
    auto is = std::upper_bound(xs.begin(),xs.end(),xcurr);
    int index = std::distance(xs.begin(), is)-1;
    xl = xs[index]; 
    xr = xs[index + 1];
    yl = ys[index]; 
    yr = ys[index + 1];
    eta = yl + (xcurr - xl)/(xr - xl)*(yr - yl);

    return eta;
}

double wave_lib_ssgw::wave_fi(lexer *p, double x, double y, double z)
{
    // Transform x location to current position xcurr at time instance p->simtime
    xcurr = fabs(modulo(x - ParameterValue.phaseVelocity*p->simtime, wL));
    xcurr = xcurr>0.5*wL ? xcurr-wL : xcurr;

    // Linear interpolation
    auto is = std::upper_bound(xs.begin(),xs.end(),xcurr);
    int index = std::distance(xs.begin(), is)-1;
    xl = xs[index]; 
    xr = xs[index + 1];
    yl = phis[index]; 
    yr = phis[index + 1];
    fi = yl + (xcurr - xl)/(xr - xl)*(yr - yl);

    return fi;
}

double wave_lib_ssgw::wave_u(lexer *p, double x, double y, double z)
{
    double vel = 0.0;

    return cosgamma*vel;
}

double wave_lib_ssgw::wave_v(lexer *p, double x, double y, double z)
{
    double vel = 0.0;

    return singamma*vel;
}

double wave_lib_ssgw::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel = 0.0;

    return vel;
}

double wave_lib_ssgw::wave_w(lexer *p, double x, double y, double z)
{
    double vel = 0.0;

    return vel;
}

void wave_lib_ssgw::parameters(lexer *p, ghostcell *pgc){}

void wave_lib_ssgw::wave_prestep(lexer *p, ghostcell *pgc){}
