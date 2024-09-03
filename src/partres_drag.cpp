/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "partres.h"
#include "particles_obj.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"

    /// @brief Calculate drag force parameter
double partres::drag_model(lexer* p, double d, double du, double dv, double dw, double thetas) const
{
        double thetaf = 1.0-thetas;
        if(thetaf>1.0-theta_crit) // Saveguard
        thetaf=1.0-theta_crit;

        const double dU=sqrt(du*du+dv*dv+dw*dw);
        if(dU==0) // Saveguard
        return 0;

        const double Rep=dU*d*invKinVis;

        // const double Cd=24.0*(pow(thetaf,-2.65)+pow(Rep,2.0/3.0)*pow(thetaf,-1.78)/6.0)/Rep;
        const double Cd=24.0/Rep+4.0/sqrt(Rep)+0.4;
        const double Dp=Cd*3.0*drho*dU/d/4.0;

        return Dp;
}

double partres::drag_coefficient(double Re_p) const
{
    double Cd;
    if(Re_p<0.1)
    Cd = 24.0/Re_p;
    else if(Re_p<1.0)
    Cd = 22.73/Re_p+0.0903/pow(Re_p,2),+ 3.69;
    else if(Re_p<10.0)
    Cd = 29.1667/Re_p-3.8889/pow(Re_p,2)+ 1.222;
    else if(Re_p<100.0)
    Cd = 46.5/Re_p- 116.67/pow(Re_p,2)+0.6167;
    else if(Re_p<1000.0)
    Cd = 98.33/Re_p- 2778/pow(Re_p,2) +0.3644;
    else if(Re_p<5000.0)
    Cd = 148.62/Re_p-4.75e4/pow(Re_p,2)+0.357;
    else if(Re_p<10000.0)
    Cd = -490.546/Re_p +57.87e4/pow(Re_p,2) +0.46;
    // else if(Re_p<50000.0)
    else
    Cd = -1662.5/Re_p+5.4167e6/pow(Re_p,2)+0.5191;

    return Cd;
}