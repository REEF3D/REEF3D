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
#include "sediment_fdm.h"

void partres::erosion(lexer *p, fdm &a, particles_obj &PP, sediment_fdm &s)
{
        double shear_eff;
        double shear_crit;
        int counter=0;
        double topoDist,u,v,w,du,dv,dw;
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==0)
            {

                    if(p->Q200==1)
                    {
                        // if(p->ccipol4_b(a.solid,PP.X[n],PP.Y[n],PP.Z[n])<0.6)
                        if(PP.X[n]>=p->S71 && PP.X[n]<=p->S72)
                        if(PP.Y[n]>=p->S77_xs && PP.Y[n]<=p->S77_xe)
                        {
                            shear_eff=p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]);
                            shear_crit=p->ccslipol4(s.tau_crit,PP.X[n],PP.Y[n]);
                            if(shear_eff>shear_crit)
                               if(p->ccipol4_b(a.topo,PP.X[n],PP.Y[n],PP.Z[n])+2.0*PP.d50>0)
                                {
                                    PP.Flag[n]=1;
                                    ++counter;

                                    PP.shear_eff[n]=shear_eff;
                                    PP.shear_crit[n]=shear_crit;

                                    i=p->posc_i(PP.X[n]);
                                    j=p->posc_j(PP.Y[n]);
                                    bedChange[IJ] -= PP.ParcelFactor[n];
                                }
                        }
                    }
                    
                    if(p->Q200==2)
                    {
                        relative_velocity(p,a,PP,n,du,dv,dw);
                        const double dU=sqrt(du*du+dv*dv+dw*dw);
                        const double Re_p = dU*PP.d50/(p->W2/p->W1);
                        const double Cd = drag_coefficient(Re_p);
                        const double Fd = Cd * PI/8.0 * pow(PP.d50,2) * p->W1 * pow(dU,2);
                        const double mu_s = tan(p->S81);
                        const double W = p->W1 * sqrt(p->W20*p->W20+p->W21*p->W21+p->W22*p->W22) * (p->S22/p->W1-1) * PI/6.0 *pow(PP.d50,3);
                        const double Fs = W * mu_s;
                        if(Fd > W * (mu_s*cos(s.teta[IJ])-sin(s.teta[IJ])))
                        {
                            PP.Flag[n]=1;
                            ++counter;
                            i=p->posc_i(PP.X[n]);
                            j=p->posc_j(PP.Y[n]);
                            bedChange[IJ] -= PP.ParcelFactor[n];
                        }
                    }
            }
            
            //for(size_t n=0;n<PP.loopindex;n++)
            //PP.Flag[n]=1;
        // if(counter>0)
        // cout<<"On rank "<<p->mpirank<<" were "<<counter<<" particles erosiond."<<endl;
}