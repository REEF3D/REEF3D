/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double force::morison_cd(lexer* p, fdm *a, ghostcell *pgc, double Re, double Kc)
{
    cd=0.0;
    double cd6,cd8,cd10,cd15,cd20,cd40,cd60,cd100;

    if(Kc<=8.0)
    {
    cd6 = - 5.373e-51*pow(Re,9.0)
          + 2.486e-44*pow(Re,8.0)
          - 4.879e-38*pow(Re,7.0)
          + 5.300e-32*pow(Re,6.0)
          - 3.941e-26*pow(Re,5.0)
          + 1.437e-20*pow(Re,4.0)
          - 3.677e-15*pow(Re,3.0)
          + 5.605e-10*pow(Re,2.0)
          - 4.507e-06*Re
          + 1.913;
    cd6 = MIN(cd6,1.5);
    }

    if(Kc<=10.0 && Kc>6.0)
    {
    cd8 = - 3.111e-51*pow(Re,9.0)
          + 1.461e-44*pow(Re,8.0)
          - 2.947e-38*pow(Re,7.0)
          + 3.340e-32*pow(Re,6.0)
          - 2.336e-26*pow(Re,5.0)
          + 1.040e-20*pow(Re,4.0)
          - 2.927e-15*pow(Re,3.0)
          + 4.978e-10*pow(Re,2.0)
          - 4.491e-05*Re
          + 2.203;
    cd8 = MIN(cd8,1.765);
    }

    if(Kc<=15.0 && Kc>8.0)
    {
    cd10 = - 6.845e-51*pow(Re,9.0)
           + 3.125e-44*pow(Re,8.0)
           - 6.067e-38*pow(Re,7.0)
           + 6.540e-32*pow(Re,6.0)
           - 4.295e-26*pow(Re,5.0)
           + 1.775e-20*pow(Re,4.0)
           - 4.605e-15*pow(Re,3.0)
           + 7.212e-10*pow(Re,2.0)
           - 6.037e-05*Re
           + 2.69;
    cd10 = MIN(cd10,2.0);
    }

    if(Kc<=20.0 && Kc>10.0)
    {
    cd15 = - 2.991e-51*pow(Re,9.0)
           + 1.506e-44*pow(Re,8.0)
           - 3.222e-38*pow(Re,7.0)
           + 3.824e-32*pow(Re,6.0)
           - 2.755e-26*pow(Re,5.0)
           + 1.243e-20*pow(Re,4.0)
           - 3.505e-15*pow(Re,3.0)
           + 5.985e-10*pow(Re,2.0)
           - 5.585e-05*Re
           + 2.804;
    cd15 = MIN(cd15,2.0);
    }

    if(Kc<=40.0 && Kc>15.0)
    {
    cd20 = - 1.184e-51*pow(Re,9.0)
           + 6.604e-45*pow(Re,8.0)
           - 1.561e-38*pow(Re,7.0)
           + 2.044e-32*pow(Re,6.0)
           - 1.625e-26*pow(Re,5.0)
           + 8.103e-21*pow(Re,4.0)
           - 2.516e-15*pow(Re,3.0)
           + 4.660e-10*pow(Re,2.0)
           - 4.607e-05*Re
           + 2.481;
    cd20 = MIN(cd20,1.96);
    }

    if(Kc<=60.0 && Kc>20.0)
    cd40 = 2.480e-46*pow(Re,8.0)
         - 1.277e-39*pow(Re,7.0)
         + 2.762e-33*pow(Re,6.0)
         - 3.263e-27*pow(Re,5.0)
         + 2.290e-21*pow(Re,4.0)
         - 9.706e-16*pow(Re,3.0)
         + 2.406e-10*pow(Re,2.0)
         - 3.133e-05*Re
         + 2.185;

    if(Kc<=100.0 && Kc>40.0)
    {
    cd60 = 3.705e-46*pow(Re,8.0)
         - 1.644e-39*pow(Re,7.0)
         + 3.065e-33*pow(Re,6.0)
         - 3.133e-27*pow(Re,5.0)
         + 1.933e-21*pow(Re,4.0)
         - 7.503e-16*pow(Re,3.0)
         + 1.841e-10*pow(Re,2.0)
         - 2.644e-05*Re
         + 2.187;
    cd60 = MIN(cd60,1.8);
    }

    if(Kc>60.0)
    {
    cd100 = - 1.430e-51*pow(Re,9.0)
            + 6.290e-45*pow(Re,8.0)
            - 1.171e-38*pow(Re,7.0)
            + 1.204e-32*pow(Re,6.0)
            - 7.533e-27*pow(Re,5.0)
            + 3.004e-21*pow(Re,4.0)
            - 7.962e-16*pow(Re,3.0)
            + 1.507e-10*pow(Re,2.0)
            - 2.096e-05*Re
            + 2.097;
    cd100 = MIN(cd100,1.4);
    }




    if(Kc<=6.0)
    cd = cd6;

    if(Kc<=8.0 && Kc>6.0)
    cd = ((8.0-Kc)/2.0)*cd6 + ((Kc-6.0)/2.0)*cd8;

    if(Kc<=10.0 && Kc>8.0)
    cd = ((10.0-Kc)/2.0)*cd8 + ((Kc-8.0)/2.0)*cd10;

    if(Kc<=15.0 && Kc>10.0)
    cd = ((15.0-Kc)/5.0)*cd10 + ((Kc-10.0)/5.0)*cd15;

    if(Kc<=20.0 && Kc>15.0)
    cd = ((20.0-Kc)/5.0)*cd15 + ((Kc-15.0)/5.0)*cd20;

    if(Kc<=40.0 && Kc>20.0)
    cd = ((40.0-Kc)/20.0)*cd20 + ((Kc-20.0)/20.0)*cd40;

    if(Kc<=60.0 && Kc>40.0)
    cd = ((60.0-Kc)/20.0)*cd40 + ((Kc-40.0)/20.0)*cd60;

    if(Kc<=100.0 && Kc>60.0)
    cd = ((100.0-Kc)/40.0)*cd60 + ((Kc-60.0)/40.0)*cd100;

    if(Kc>100.0)
    cd = cd100;


    return cd;

}




