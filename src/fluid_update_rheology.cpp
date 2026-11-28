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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"interface_width.h"
#include"fluid_update_rheology.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"rheology_f.h"
#include"rheology_cfd_core.h"
#include"heaviside.h"
#include"rheology_cfd_3phase.h"

fluid_update_rheology::fluid_update_rheology(lexer *p) : ro1(p->W1), ro2(p->W3), visc2(p->W4)
{
    iter=0;
    iocheck = true;
    pcore = nullptr;
    
    // rheology core (W 130), else the legacy model (W 90)
    if(p->W130>0)
    prheo = pcore = new rheology_cfd_core(p);
    
    else
    prheo = new rheology_f(p);

}

fluid_update_rheology::~fluid_update_rheology()
{
    delete prheo;
}

void fluid_update_rheology::start(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    const int gcval_ro = 1;
    const int gcval_visc = 1;

    double H_phi=0.0;
    p->volume1=0.0;
    p->volume2=0.0;
    
    if(p->count>iter)
    iocheck = 0;
    
    iter=p->count;
    
    // three-phase slide model (W 140)
    if(rheology_cfd_3phase::self!=nullptr && pcore!=nullptr)
    {
    start_3phase(p,a,pgc,u,v,w);
    return;
    }

    // W 131 3: the band treatment, for a material with a yield stress
    const bool band = (pcore!=nullptr && p->W131==3 && pcore->plastic(p,a));
    
    // density, viscosity & volumes
    LOOP
    {  
        const double eps = interface_width(p,a->phi,i,j,k);
        H_phi = heaviside(a->phi(i,j,k),eps);
        
        // W 131 3 (rheology core): the material's viscosity over the whole interface band (phi > -eps), not its
        // share H. The band is a heavy mixture with a small share of the material on its air side; with the
        // viscosity blended by H its strength was H mu_s p' against the weight of the mixture, and on an inclined
        // surface of an arrested deposit the band flowed downslope as a thin gravity current (0.15 m/s) that
        // carried the level set and let the deposit slump. With the full viscosity and the pressure of the band
        // (rheology_cfd_core: W 131 3) the strength is mu_s p' against the weight: the band holds below the
        // friction angle. At phi = -eps the material's viscosity is mu_min (p' = 0)
        // Only for a material with a yield stress (as in the three-phase model): a Newtonian or power-law material
        // has no deposit to hold, and its full viscosity in the light outer band cells made the explicit cross terms
        // of the diffusion unstable (case 12)
        double Hv = H_phi;
        
        if(band && a->phi(i,j,k)>-eps)
        Hv = 1.0;

        a->ro(i,j,k) = ro1*H_phi + ro2*(1.0-H_phi);

        visc1 = prheo->viscosity(p,a,pgc,u,v,w);
        a->visc(i,j,k) = visc1*Hv + visc2*(1.0-Hv);
        
        // dynamic form of the diffusion (D 24): rho nu is the blend of the dynamic viscosities
        if(p->D24>0)
        a->visc(i,j,k) = (ro1*visc1*Hv + ro2*visc2*(1.0-Hv))/a->ro(i,j,k);

        if(p->flagsf4[IJK]>0)
        {
            p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H_phi*a->porosity(i,j,k);
            p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H_phi)*a->porosity(i,j,k);
        }
    }

    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    p->volume1 = pgc->globalsum(p->volume1);
    p->volume2 = pgc->globalsum(p->volume2);

    
    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
        cout<<"Volume 1: "<<p->volume1<<endl;
        cout<<"Volume 2: "<<p->volume2<<endl;
    }
    
    ++iocheck;
}

// three-phase slide model (W 140): mixture of the slide (fraction alpha_s, bulk density W 141, the viscosity of
// the core) and the fluid of the level set (water W 1 W 2, air W 3 W 4), see rheology_cfd_3phase.h
void fluid_update_rheology::start_3phase(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    rheology_cfd_3phase *s3 = rheology_cfd_3phase::self;
    
    const int gcval_ro = 1;
    const int gcval_visc = 1;
    const double rhow = p->W1, rhoa = p->W3;
    const double muw = p->W1*p->W2, mua = p->W3*p->W4;
    const bool plastic = pcore->plastic(p,a);
    
    p->volume1=0.0;
    p->volume2=0.0;
    
    LOOP
    {
        const double H = heaviside(a->phi(i,j,k),interface_width(p,a->phi,i,j,k));
        const double al = s3->alphas(i,j,k);
        
        const double rof = rhow*H + rhoa*(1.0-H);
        const double muf = muw*H + mua*(1.0-H);
        
        a->ro(i,j,k) = al*s3->rho_slide(i,j,k) + (1.0-al)*rof;
        
        double mu = muf;
        
        if(al>1.0e-8)
        {
            const double mus = pcore->mu_eff(p,a,u,v,w);
            
            if(p->W149==2)
            mu = 1.0/(al/mus + (1.0-al)/muf);
            
            else
            mu = al*mus + (1.0-al)*muf;
            
            // W 131 3: the full viscosity of the slide over the smoothed band for a material with a yield stress
            // (rheology_cfd_3phase.h). Not for a Newtonian or power-law slide: it has no deposit to hold, and with
            // its full viscosity (75 Pa s) in the light outer band cells (rho 3 kg/m3) the explicit cross terms of
            // the diffusion blew up. The criterion must not be the state (stress below the yield stress): the
            // slumping band cells flow, so they stayed weak
            if(p->W131==3 && al>1.0e-3 && plastic)
            mu = MAX(mus,muf);
        }
        
        a->visc(i,j,k) = mu/a->ro(i,j,k);
        
        if(p->flagsf4[IJK]>0)
        {
            const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*a->porosity(i,j,k);
            p->volume1 += dV*(1.0-al)*H;
            p->volume2 += dV*(1.0-al)*(1.0-H);
        }
    }
    
    pgc->start4(p,a->ro,gcval_ro);
    pgc->start4(p,a->visc,gcval_visc);
    p->volume1 = pgc->globalsum(p->volume1);
    p->volume2 = pgc->globalsum(p->volume2);
    
    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
        cout<<"Volume 1: "<<p->volume1<<"  (water)"<<endl;
        cout<<"Volume 2: "<<p->volume2<<"  (air)"<<endl;
    }
    
    ++iocheck;
}
