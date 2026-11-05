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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

/*--------------------------------------------------------------------
seepage flow in the bed (Q 69 1, S 10 1)

With S 10 1 the bed is a solid boundary for the fluid, the fluid solves no pressure in the
cells of its bed. Without seepage the pore pressure is hydrostatic. With Q 69 the pore
pressure follows Darcy's law in a rigid bed of uniform permeability: the piezometric pressure

    p* = p - rho_f g.x

is harmonic in the bed, div(grad p*) = 0, with p* of the fluid at the bed surface (Dirichlet at
the centres of the adjacent fluid cells), no flux through walls, solid bodies and the bottom,
or with Q 70 an excess pore pressure p*_s + Q 70 at the bottom (upflow tests; p*_s is the mean p*
of the fluid at the bed surface). The Laplace equation is solved by red-black SOR sweeps once per
time step, warm-started from the last step, until the change is below 1e-6 of the range of p*.

The seepage acts on the grains through the excess pore pressure. Per unit volume of the bed the
seepage force is -grad(p*) (Terzaghi); without a pore flow in the fluid (S 10 1) the parcels
take it as a whole, per unit mass of solid -grad(p*)/(rho_s theta). The contact network carries
the effective load: the overburden stress (stress_overburden) integrates

    (rho_s - rho_f) g theta + d(p*)/dz

downwards, so an upward gradient unloads the bed; where the excess gradient exceeds the
submerged weight (i = (s-1)(1-n) for a uniform bed), the effective stress vanishes, the bed
is fluidised. The jammed bed (Q 64) holds only where the contact network is loaded: the jam
weight is reduced with the ratio of the overburden with and without the seepage (Pov/Pnos,
full below 5 %, none above 30 %).
--------------------------------------------------------------------*/

// a cell of the bed for the seepage: inside the fluid's bed, not in a solid body
bool CPM::seep_bed(lexer *p, fdm *a, int ii, int jj, int kk)
{
    return a->topo(ii,jj,kk)<0.0 && !(p->solidread>0 && a->solid(ii,jj,kk)<0.0);
}

// piezometric pressure of the fluid in the cell (ii,jj,kk)
double CPM::pstar_fluid(lexer *p, fdm *a, int ii, int jj, int kk)
{
    return a->press(ii,jj,kk) - p->W1*(p->W20*p->XP[ii+marge] + p->W21*p->YP[jj+marge] + p->W22*p->ZP[kk+marge]);
}

void CPM::seepage_update(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->Q69<1 || p->S10==2)
    return;

    // reference: mean p* of the fluid cells on top of the bed; range for the tolerance
    double psum=0.0, pmin=1.0e20, pmax=-1.0e20;
    int pcount=0;

    BASELOOP
    if(seep_bed(p,a,i,j,k) && k+1<p->knoz && !seep_bed(p,a,i,j,k+1) && !(p->solidread>0 && a->solid(i,j,k+1)<0.0))
    {
        double ps = pstar_fluid(p,a,i,j,k+1);
        psum += ps;
        pmin = MIN(pmin,ps);
        pmax = MAX(pmax,ps);
        ++pcount;
    }

    psum = pgc->globalsum(psum);
    pcount = pgc->globalisum(pcount);
    pmin = pgc->globalmin(pmin);
    pmax = pgc->globalmax(pmax);

    if(pcount==0)
    return;

    const double pref = psum/double(pcount);
    const double pbot = pref + p->Q70;
    const double tol = 1.0e-6*MAX(pmax-pmin+fabs(p->Q70), 1.0);

    // start: p* of the fluid on top of the column (hydrostatic pore pressure)
    if(seep_ini==0)
    {
        ILOOP
        JLOOP
        {
            double ptop = pref;

            for(k=p->knoz-1; k>=0; --k)
            {
                if(p->flag4[IJK]<=OBJ_FLAG)
                continue;

                if(!seep_bed(p,a,i,j,k))
                {
                    if(!(p->solidread>0 && a->solid(i,j,k)<0.0))
                    ptop = pstar_fluid(p,a,i,j,k);

                    Pse(i,j,k) = ptop;
                }
                else
                Pse(i,j,k) = ptop;
            }
        }

        pgc->start4a(p,Pse,1);
        seep_ini=1;
    }

    // the fluid cells keep their p* (Dirichlet values for the bed)
    BASELOOP
    if(!seep_bed(p,a,i,j,k))
    Pse(i,j,k) = pstar_fluid(p,a,i,j,k);

    pgc->start4a(p,Pse,1);

    const double omega = 1.8;

    // at most 200 sweeps per time step (warm start: a few after the first steps)
    for(int it=0; it<200; ++it)
    {
        double change = 0.0;

        for(int color=0; color<2; ++color)
        {
            BASELOOP
            if(((p->origin_i+i + p->origin_j+j + p->origin_k+k)&1)==color && seep_bed(p,a,i,j,k))
            {
                double asum=0.0, bsum=0.0;

                auto nb = [&](int ii, int jj, int kk, double coef, bool bottom)
                {
                    // outside the domain: wall, or the bottom with the excess pressure Q 70
                    if(wallcell(p,ii,jj,kk))
                    {
                        if(bottom && fabs(p->Q70)>0.0)
                        {
                            asum += 2.0*coef;
                            bsum += 2.0*coef*pbot;
                        }
                        return;
                    }

                    // solid body: no flux
                    if(p->solidread>0 && a->solid(ii,jj,kk)<0.0)
                    return;

                    asum += coef;
                    bsum += coef*Pse(ii,jj,kk);
                };

                nb(i-1,j,k, 1.0/(p->DXN[IP]*p->DXP[IM1]), false);
                nb(i+1,j,k, 1.0/(p->DXN[IP]*p->DXP[IP]), false);

                if(p->j_dir==1)
                {
                nb(i,j-1,k, 1.0/(p->DYN[JP]*p->DYP[JM1]), false);
                nb(i,j+1,k, 1.0/(p->DYN[JP]*p->DYP[JP]), false);
                }

                nb(i,j,k-1, 1.0/(p->DZN[KP]*(k-1<0 ? p->DZN[KP] : p->DZP[KM1])), k-1<0);
                nb(i,j,k+1, 1.0/(p->DZN[KP]*p->DZP[KP]), false);

                if(asum>0.0)
                {
                    double pn = (1.0-omega)*Pse(i,j,k) + omega*bsum/asum;
                    change = MAX(change, fabs(pn-Pse(i,j,k)));
                    Pse(i,j,k) = pn;
                }
            }

            pgc->start4a(p,Pse,1);
        }

        change = pgc->globalmax(change);

        if(change<tol)
        break;
    }

    // vertical gradient of p* in the bed, for the effective stress of the contact network
    BASELOOP
    {
        Gsz(i,j,k) = 0.0;

        if(seep_bed(p,a,i,j,k))
        {
            bool okm = !wallcell(p,i,j,k-1) && !(p->solidread>0 && a->solid(i,j,k-1)<0.0);
            bool okp = !wallcell(p,i,j,k+1) && !(p->solidread>0 && a->solid(i,j,k+1)<0.0);

            double pm = okm ? Pse(i,j,k-1) : Pse(i,j,k);
            double pp = okp ? Pse(i,j,k+1) : Pse(i,j,k);

            // bottom with the excess pressure: the face value
            if(!okm && k-1<0 && p->nb5<0 && fabs(p->Q70)>0.0)
            pm = 2.0*pbot - Pse(i,j,k);

            double dz = (okm||(k-1<0 && fabs(p->Q70)>0.0) ? p->DZP[KM1] : 0.0) + (okp ? p->DZP[KP] : 0.0);

            Gsz(i,j,k) = dz>0.0 ? (pp-pm)/dz : 0.0;
        }
    }

    pgc->start4a(p,Gsz,1);
}
