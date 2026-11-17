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
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"patchBC_interface.h"
#include<vector>

void iowave::inflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    if(p->I230==0)
    {
    if(p->B98==0 || (p->W10>0.0 && p->B98==2))
    inflow_plain(p,a,pgc,u,v,w);

    
    // B 100 1: the relaxation method with the wave flux through the inlet face
    if(p->B98==2 && p->B100==1 && p->W10<=0.0)
    relax_inlet_flux(p,a,pgc,u);

	if(p->B98==3)
	dirichlet_wavegen(p,a,pgc,u,v,w);
	
	if(p->B98==4)
	active_wavegen(p,a,pgc,u,v,w);
	}
    
	if(p->B99==3||p->B99==4||p->B99==5)
	active_beach(p,a,pgc,u,v,w);
    
    if(p->I230>0)
    ff_inflow(p,a,pgc,u,v,w);
    
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void iowave::rkinflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
    u(i-1,j,k) = u(i-2,j,k) = u(i-3,j,k) = a->u(i-1,j,k);
    v(i-1,j,k) = v(i-2,j,k) = v(i-3,j,k) = a->v(i-1,j,k);
    w(i-1,j,k) = w(i-2,j,k) = w(i-3,j,k) = a->w(i-1,j,k);
    }
    
    if(p->B99==3||p->B99==4||p->B99==5)
	active_beach(p,a,pgc,u,v,w);
    
    pBC->patchBC_rkioflow(p,a,pgc,u,v,w);
}

void iowave::inflow_plain(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        u(i-1,j,k)=p->Ui;
        u(i-2,j,k)=p->Ui;
        u(i-3,j,k)=p->Ui;
		
		v(i-1,j,k)=0.0;
        v(i-2,j,k)=0.0;
        v(i-3,j,k)=0.0;
		
		w(i-1,j,k)=0.0;
        w(i-2,j,k)=0.0;
        w(i-3,j,k)=0.0;
        
        // Air inflow
        if(p->W50_air==1 && a->phi(i,j,k)<-0.6*p->DXM)
        {
        u(i-1,j,k)+=p->W50;
        u(i-2,j,k)+=p->W50;
        u(i-3,j,k)+=p->W50;
        }
    
    }
}

void iowave::flowfile(lexer *p, fdm* a, ghostcell* pgc, turbulence *pturb)
{
    if(p->gcin_count>0 && p->I230>0)
    flowfile_start(p,a,pgc,pturb);
}

// B 100 1, relaxation wave generation (B 98 2): the inlet face is no wall.
// With a wall there, the relaxed velocity of the first cells (the wave theory) has the divergence
// u_wave/dx in the first column of cells, in every step again. The projection removes it with a
// pressure of the order rho u_wave L/dt at the inlet (L about the water depth): this pressure
// grows as the time step falls and is not part of the wave.
// Here the inlet face in the water takes the velocity of the wave theory at the face (the same
// Heaviside function as u_relax). The face in the air takes the uniform velocity that balances the
// water flux, so that the net flux through the inlet is 0 in every step: the domain stays closed
// for the projection (no reference pressure drift, no compression of the air with W 32).
void iowave::relax_inlet_flux(lexer *p, fdm *a, ghostcell *pgc, field &u)
{
    const int nin = p->gcin_count;
    std::vector<double> uw(nin,0.0), ha(nin,0.0);
    double Qw=0.0, Aa=0.0;

    for(n=0;n<nin;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        phival = a->phi(i,j,k);

        if(phival>=-psi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(phival>=-epsi && phival<-psi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));

        const double zloc = p->pos_z();
        const double fsfloc = eta(i,j) + p->phimean;

        z = (zloc<=fsfloc) ? zloc-p->phimean : eta(i,j);

        double uwave = 0.0;
        if(zloc<=fsfloc+epsi)
        uwave = wave_u(p,pgc,xgen(p)-0.5*p->DXN[IP],ygen(p),z) + p->Ui;

        const double dA = p->DYN[JP]*p->DZN[KP];

        uw[n] = ramp(p)*uwave*H;
        ha[n] = 1.0-H;

        Qw += uw[n]*dA;
        Aa += ha[n]*dA;
    }

    Qw = pgc->globalsum(Qw);
    Aa = pgc->globalsum(Aa);

    const double ua = Aa>1.0e-12 ? -Qw/Aa : 0.0;

    for(n=0;n<nin;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        u(i-1,j,k) = u(i-2,j,k) = u(i-3,j,k) = uw[n] + ha[n]*ua;
    }
}
