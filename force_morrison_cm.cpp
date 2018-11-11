/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

double force::morison_cm(lexer* p, fdm *a, ghostcell *pgc, double Re, double Kc)
{
    cm=0.0;
    double cm6,cm8,cm10,cm15,cm20,cm40,cm60,cm100;

    if(Kc<=8.0)
    cm6 = 5.795e-52*pow(Re,9.0)
        - 3.308e-45*pow(Re,8.0)
        + 7.959e-39*pow(Re,7.0)
        - 1.050e-32*pow(Re,6.0)
        + 8.253e-27*pow(Re,5.0)
        - 3.938e-21*pow(Re,4.0)
        + 1.103e-15*pow(Re,3.0)
        - 1.625e-10*pow(Re,2.0)
        + 8.968e-06*Re
        + 1.801;

    if(Kc<=10.0 && Kc>6.0)
    cm8 = - 1.097e-45*pow(Re,8.0)
          + 4.876e-39*pow(Re,7.0)
          - 9.013e-33*pow(Re,6.0)
          + 8.952e-27*pow(Re,5.0)
          - 5.144e-21*pow(Re,4.0)
          + 1.706e-15*pow(Re,3.0)
          - 3.039e-10*pow(Re,2.0)
          + 2.362e-05*Re
          + 1.226;

    if(Kc<=15.0 && Kc>8.0)
    cm10 = - 1.579e-45*pow(Re,8.0)
           + 7.073e-39*pow(Re,7.0)
           - 1.324e-32*pow(Re,6.0)
           + 1.340e-26*pow(Re,5.0)
           - 7.922e-21*pow(Re,4.0)
           + 2.745e-15*pow(Re,3.0)
           - 3.039e-10*pow(Re,2.0)
           + 2.362e-05*Re
           + 0.2974;

    if(Kc<=20.0 && Kc>10.0)
    cm15 = - 1.955e-46*pow(Re,8.0)
           + 1.275e-39*pow(Re,7.0)
           - 3.191e-33*pow(Re,6.0)
           + 4.123e-27*pow(Re,5.0)
           - 3.034e-21*pow(Re,4.0)
           + 1.296e-15*pow(Re,3.0)
           - 3.079e-10*pow(Re,2.0)
           + 3.492e-05*Re
           + 0.3195;

    if(Kc<=40.0 && Kc>15.0)
    cm20 = - 3.638e-34*pow(Re,6.0)
           + 9.752e-28*pow(Re,5.0)
           - 1.036e-21*pow(Re,4.0)
           + 5.546e-16*pow(Re,3.0)
           - 1.560e-10*pow(Re,2.0)
           + 2.149e-05*Re
           + 0.6446;

    if(Kc<=60.0 && Kc>20.0)
    cm40 = 1.754e-51*pow(Re,9.0)
         - 8.273e-45*pow(Re,8.0)
         + 1.636e-38*pow(Re,7.0)
         - 1.762e-32*pow(Re,6.0)
         + 1.129e-26*pow(Re,5.0)
         - 4.427e-21*pow(Re,4.0)
         + 1.076e-15*pow(Re,3.0)
         - 1.672e-10*pow(Re,2.0)
         + 1.753e-05*Re
         + 0.7768;

    if(Kc<=100.0 && Kc>40.0)
    cm60 = 1.368e-51*pow(Re,9.0)
         - 6.363e-45*pow(Re,8.0)
         + 1.266e-38*pow(Re,7.0)
         - 1.406e-32*pow(Re,6.0)
         + 9.520e-27*pow(Re,5.0)
         - 4.036e-21*pow(Re,4.0)
         + 1.066e-15*pow(Re,3.0)
         - 1.736e-10*pow(Re,2.0)
         + 1.802e-05*Re
         + 0.6902;

    if(Kc>60.0)
    {
    cm100 = 9.804e-52*pow(Re,9.0)
          - 4.683e-45*pow(Re,8.0)
          + 9.516e-39*pow(Re,7.0)
          - 1.075e-32*pow(Re,6.0)
          + 7.387e-27*pow(Re,5.0)
          - 3.193e-21*pow(Re,4.0)
          + 8.687e-16*pow(Re,3.0)
          - 1.470e-10*pow(Re,2.0)
          + 1.579e-05*Re
          + 0.6491;
    cm100 = MIN(cm100,1.9);
    }


    if(Kc<=6.0)
    cm = cm6;

    if(Kc<=8.0 && Kc>6.0)
    cm = ((8.0-Kc)/2.0)*cm6 + ((Kc-6.0)/2.0)*cm8;

    if(Kc<=10.0 && Kc>8.0)
    cm = ((10.0-Kc)/2.0)*cm8 + ((Kc-8.0)/2.0)*cm10;

    if(Kc<=15.0 && Kc>10.0)
    cm = ((15.0-Kc)/5.0)*cm10 + ((Kc-10.0)/5.0)*cm15;

    if(Kc<=20.0 && Kc>15.0)
    cm = ((20.0-Kc)/5.0)*cm15 + ((Kc-15.0)/5.0)*cm20;

    if(Kc<=40.0 && Kc>20.0)
    cm = ((40.0-Kc)/20.0)*cm20 + ((Kc-20.0)/20.0)*cm40;

    if(Kc<=60.0 && Kc>40.0)
    cm = ((60.0-Kc)/20.0)*cm40 + ((Kc-40.0)/20.0)*cm60;

    if(Kc<=100.0 && Kc>60.0)
    cm = ((100.0-Kc)/40.0)*cm60 + ((Kc-60.0)/40.0)*cm100;

    if(Kc>100.0)
    cm = cm100;

    return cm;

}



