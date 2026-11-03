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

#include"hydro_core.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

// the owned values of v (cells n0..n1-1) to every rank
void hydro_core::gather(ghostcell *pgc, vector<double> &v) const
{
    MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,v.data(),cnt.data(),displ.data(),MPI_DOUBLE,pgc->mpi_comm);
}

/*
kinematic wave on the channel of a cell (length L, A = alpha Q^0.6), implicit in time (nonlinear scheme of Chow, Maidment & Mays 1988, Applied Hydrology, ch. 9):
    dt/L Q + alpha Q^0.6 = dt/L (Qin + qgen) + alpha Qold^0.6
Qin: outflow of the upstream cells at the new time (they come first), qgen: lateral inflow of the cell [m3/s].
The channel volume L alpha Q^0.6 changes by dt (Qin + qgen - Q): the routing conserves the volume.
Solved with Newton, safeguarded by bisection in [0, min(R/c, (R/alpha)^(1/0.6))].
*/
double hydro_core::kinwave(double c, double a, double R) const
{
    if(R<=0.0)
    return 0.0;

    const double beta=0.6;
    double lo=0.0, hi=min(R/c,pow(R/a,1.0/beta));
    double Q=hi;

    for(int it=0; it<200; ++it)
    {
        const double f=c*Q+a*pow(Q,beta)-R;

        if(f>0.0)
        hi=Q;
        else
        lo=Q;

        const double df=c+a*beta*pow(max(Q,1.0e-300),beta-1.0);
        double Qn=Q-f/df;

        if(!(Qn>lo && Qn<hi))
        Qn=0.5*(lo+hi);

        if(fabs(Qn-Q)<=1.0e-14*Qn || hi-lo<=1.0e-15*hi)
        return Qn;

        Q=Qn;
    }

    return Q;
}

void hydro_core::route(lexer *p, ghostcell *pgc)
{
    gather(pgc,qgen);

    Qold=Q;
    fill(Qin.begin(),Qin.end(),0.0);

    // coupled (R 100 2): the outflow of SFLOW into the cells, the infiltration of SFLOW (vertical)
    double vfrom=0.0, vto=0.0;
    for(int n=0; n<N; ++n)
    {
        vfrom+=addv[n];

        if(qret[n]!=0.0)
        {
            Qin[n]+=qret[n];
            vfrom+=qret[n]*dt;
        }
    }

    double vout=0.0;
    for(int n=0; n<N; ++n)
    {
        // the share of the runoff under the SFLOW grid
        qsf[n]=share[n]*qgen[n];
        const double ql=qgen[n]-qsf[n];

        if(rescell[n]>=0)
        Q[n]=reservoir_step(p,resv[rescell[n]],Qin[n]+ql);
        else if(p->R106==0)
        Q[n]=Qin[n]+ql;
        else
        {
            const double c=dt/len[n];
            Q[n]=kinwave(c,alpha[n],c*(Qin[n]+ql)+alpha[n]*pow(Qold[n],0.6));
        }

        Qent[n]=0.0;
        if(enter[n])
        Qent[n]=Q[n];
        else if(down[n]>=0)
        Qin[down[n]]+=Q[n];
        else
        vout+=Q[n]*dt;

        vto+=(qsf[n]+Qent[n])*dt;
    }

    Vout+=vout;
    Vto+=vto;
    Vfrom+=vfrom;
}
