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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ioflow_gravity.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"patchBC_interface.h"

ioflow_gravity::ioflow_gravity(lexer *p, ghostcell *pgc, patchBC_interface *ppBC) 
{
    pBC = ppBC;
    
	omega_x = 2.0*PI*p->B191_2;
	omega_y = 2.0*PI*p->B192_2;
	
	theta_x = p->B191_1*(PI/180.0);
	theta_y = p->B192_1*(PI/180.0);
}

ioflow_gravity::~ioflow_gravity()
{
}

void ioflow_gravity::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
    int count1,count2,n;

    count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
        p->gcin[count1][0]=p->gcb4[n][0];
        p->gcin[count1][1]=p->gcb4[n][1];
        p->gcin[count1][2]=p->gcb4[n][2];
        p->gcin[count1][3]=p->gcb4[n][3];
        p->gcin[count1][4]=p->gcb4[n][5];
        ++count1;
        }

        if(p->gcb4[n][4]==2)
        {
        p->gcout[count2][0]=p->gcb4[n][0];
        p->gcout[count2][1]=p->gcb4[n][1];
        p->gcout[count2][2]=p->gcb4[n][2];
        p->gcout[count2][3]=p->gcb4[n][3];
        p->gcout[count2][4]=p->gcb4[n][5];
        ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;
}

void ioflow_gravity::discharge(lexer *p, fdm* a, ghostcell* pgc)
{
    // patchBC
    pBC->patchBC_discharge(p,a,pgc);
}


void ioflow_gravity::inflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    double omega;
    a->gi = p->W20;
    a->gj = p->W21;
    a->gk = p->W22;
    
    
	// ------- 
    // translation 
    if(p->B181==1)
    a->gi += -p->B181_1*(2.0*PI*p->B181_2)*sin((2.0*PI*p->B181_2)*p->simtime + p->B181_3);
    
    if(p->B182==1)
    a->gj += -p->B182_1*(2.0*PI*p->B182_2)*sin((2.0*PI*p->B182_2)*p->simtime + p->B182_3);

    if(p->B183==1)
    a->gk += -p->B183_1*(2.0*PI*p->B183_2)*sin((2.0*PI*p->B183_2)*p->simtime + p->B183_3);
    
    
    // -------
    // rotation
    if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
    {
	a->gj += sin(theta_x*sin(omega_x*p->simtime))*p->W22;
    
    a->gk +=  cos(theta_x*sin(omega_x*p->simtime))*p->W22 - p->W22;
    }
    
    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
    {
	a->gi += sin(theta_y*sin(omega_y*p->simtime))*p->W22;
    
    a->gk +=  cos(theta_y*sin(omega_y*p->simtime))*p->W22  - p->W22;
    }
    
    // -------
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_gravity::rkinflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_gravity::fsfinflow(lexer *p, fdm *a, ghostcell *pgc)
{
    pBC->patchBC_waterlevel(p,a,pgc,a->phi);
}

void ioflow_gravity::fsfrkout(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
}

void ioflow_gravity::fsfrkin(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
    pBC->patchBC_waterlevel(p,a,pgc,f);
}

void ioflow_gravity::fsfrkoutV(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_gravity::fsfrkinV(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_gravity::fsfrkoutVa(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_gravity::fsfrkinVa(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
}

void ioflow_gravity::iogcb_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void  ioflow_gravity::isource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
	if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	{	
		n=0;
		ULOOP
		{
			dist_x = p->pos_x() - p->B192_3;
			dist_z = p->pos_z() - p->B192_4;
			a->rhsvec.V[n] += dist_z*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
						 + dist_x*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
						 //- a->w(i,j,k)*theta_y*omega_y*cos(omega_y*p->simtime);
		++n;
		}
		
	//a->gi = sin(theta_y*sin(omega_y*p->simtime))*p->W22;
	/*a->gi = p->W22*sin( theta_y*sin(omega_y*p->simtime) - 0.31*theta_y*theta_y*(1.0+cos(2.0*omega_y*p->simtime))
						+ pow(theta_y,3.0)*(0.16*cos(omega_y*p->simtime) - 0.16*cos(3.0*omega_y*p->simtime)
							+ 0.13*sin(omega_y*p->simtime) - 0.004*sin(3.0*omega_y*p->simtime)));*/

	}
}

void  ioflow_gravity::jsource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
	//if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	//a->gj = sin(theta_x*sin(omega_x*p->simtime))*p->W22;
}

void  ioflow_gravity::ksource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	a->rhsvec.V[n]=0.0;
	
	if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	{
		n=0;
		WLOOP
		{
			dist_y = p->pos_y() - p->B191_3;
			a->rhsvec.V[n] -= dist_y*theta_x*pow(omega_x,2.0)*sin(omega_x*p->simtime);
			
		++n;
		}
		
		//a->gk =  cos(theta_x*sin(omega_x*p->simtime))*p->W22;
	}
	
	
	
	if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	{
		n=0;
		WLOOP
		{
			dist_x = p->pos_x() - p->B192_3;
			dist_z = p->pos_z() - p->B192_4;
			a->rhsvec.V[n] +=  -dist_x*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
						 +  dist_z*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
						 //- a->u(i,j,k)*theta_y*omega_y*cos(omega_y*p->simtime);	
		++n;
		}
		
		//a->gk =  cos(theta_y*sin(omega_y*p->simtime))*p->W22;
		/*
		a->gk = p->W22*cos( theta_y*sin(omega_y*p->simtime) - 0.31*theta_y*theta_y*(1.0+cos(2.0*omega_y*p->simtime))
						+ pow(theta_y,3.0)*(0.16*cos(omega_y*p->simtime) - 0.16*cos(3.0*omega_y*p->simtime)
							+ 0.13*sin(omega_y*p->simtime) - 0.004*sin(3.0*omega_y*p->simtime)));*/

	}
}

void  ioflow_gravity::isource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	d->rhsvec.V[n]=0.0;
	
	if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	{	
		n=0;
		LOOP
		{
			dist_x = p->pos_x() - p->B192_3;
			dist_z = p->pos_z() - p->B192_4;
			d->rhsvec.V[n] += dist_z*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
						 + dist_x*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
		++n;
		}
		
	d->gi = sin(theta_y*sin(omega_y*p->simtime))*p->W22;
	}
}

void  ioflow_gravity::jsource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	d->rhsvec.V[n]=0.0;
	
	if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	d->gj = sin(theta_x*sin(omega_x*p->simtime))*p->W22;
}

void  ioflow_gravity::ksource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
	NLOOP4
	d->rhsvec.V[n]=0.0;
	
	if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	{
		n=0;
		LOOP
		{
			dist_y = p->pos_y() - p->B191_3;
			d->rhsvec.V[n] -= dist_y*theta_x*pow(omega_x,2.0)*sin(omega_x*p->simtime);
			
		++n;
		}
		
		//d->gk =  cos(theta_x*sin(omega_x*p->simtime))*p->W22;
	}
	
	
	
	if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)	
	{
		n=0;
		WLOOP
		{
			dist_x = p->pos_x() - p->B192_3;
			dist_z = p->pos_z() - p->B192_4;
			d->rhsvec.V[n] +=  -dist_x*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
						 +  dist_z*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
						 //- a->u(i,j,k)*theta_y*omega_y*cos(omega_y*p->simtime);	
		++n;
		}
		
		//d->gk =  cos(theta_y*sin(omega_y*p->simtime))*p->W22;
		/*
		a->gk = p->W22*cos( theta_y*sin(omega_y*p->simtime) - 0.31*theta_y*theta_y*(1.0+cos(2.0*omega_y*p->simtime))
						+ pow(theta_y,3.0)*(0.16*cos(omega_y*p->simtime) - 0.16*cos(3.0*omega_y*p->simtime)
							+ 0.13*sin(omega_y*p->simtime) - 0.004*sin(3.0*omega_y*p->simtime)));*/

	}
}


void ioflow_gravity::pressure_io(lexer *p, fdm *a, ghostcell *pgc)
{
    pBC->patchBC_pressure(p,a,pgc,a->press);
}

void ioflow_gravity::turbulence_io(lexer *p, fdm* a, ghostcell* pgc)
{
}

void ioflow_gravity::u_relax(lexer *p, fdm *a, ghostcell *pgc, field &uvel)
{
}

void ioflow_gravity::v_relax(lexer *p, fdm *a, ghostcell *pgc, field &vvel)
{
}

void ioflow_gravity::w_relax(lexer *p, fdm *a, ghostcell *pgc, field &wvel)
{
}

void ioflow_gravity::p_relax(lexer *p, fdm *a, ghostcell *pgc, field &press)
{
}

void ioflow_gravity::phi_relax(lexer *p, ghostcell *pgc, field &f)
{
}

void ioflow_gravity::vof_relax(lexer *p, ghostcell *pgc, field &f)
{
}

void ioflow_gravity::turb_relax(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
}

void ioflow_gravity::U_relax(lexer *p, ghostcell *pgc, double *U, double *UH)
{
}

void ioflow_gravity::V_relax(lexer *p, ghostcell *pgc, double *V, double *VH)
{
}

void ioflow_gravity::W_relax(lexer *p, ghostcell *pgc, double *W, double *WH)
{
}

void ioflow_gravity::P_relax(lexer *p, ghostcell *pgc, double *P)
{
}

void ioflow_gravity::WL_relax(lexer *p, ghostcell *pgc, slice &WL, slice &depth)
{
}

void ioflow_gravity::fi_relax(lexer *p, ghostcell *pgc, field &f, field &phi)
{
}

void ioflow_gravity::fivec_relax(lexer *p, ghostcell *pgc, double *f)
{
}
    
void ioflow_gravity::fifsf_relax(lexer *p, ghostcell *pgc, slice& f)
{
}

void ioflow_gravity::visc_relax(lexer *p, ghostcell *pgc, slice& f)
{
}


void ioflow_gravity::eta_relax(lexer *p, ghostcell *pgc, slice &f)
{
}

void ioflow_gravity::um_relax(lexer *p, ghostcell *pgc, slice &P, slice &bed, slice &eta)
{
}

void ioflow_gravity::vm_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
}

void ioflow_gravity::wm_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
}

void ioflow_gravity::ws_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
}

void ioflow_gravity::pm_relax(lexer *p, ghostcell *pgc, slice &f)
{
}

double ioflow_gravity::wave_fsf(lexer *p, ghostcell *pgc, double x)
{
    double val=0.0;

    return val;
}

int ioflow_gravity::iozonecheck(lexer *p, fdm*a)
{	
	int check = 1;
	
	return check;
}

void ioflow_gravity::inflow_walldist(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, reini *preini, ioflow *pflow)
{
}

void ioflow_gravity::discharge2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
}

void ioflow_gravity::Qin2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
}

void ioflow_gravity::Qout2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
}
void ioflow_gravity::inflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_gravity::rkinflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_gravity::isource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
	SLICELOOP1
	b->F(i,j)=0.0;
}

void ioflow_gravity::jsource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
	SLICELOOP2
	b->G(i,j)=0.0;
}

void ioflow_gravity::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    ALOOP
	a->porosity(i,j,k)=1.0;
	
	pgc->start4a(p,a->porosity,1);
}

void ioflow_gravity::full_initialize2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
	
}

void ioflow_gravity::flowfile(lexer *p, fdm* a, ghostcell* pgc, turbulence *pturb)
{
}

void ioflow_gravity::wavegen_precalc(lexer *p, ghostcell *pgc)
{
    
}

void ioflow_gravity::wavegen_precalc_ini(lexer *p, ghostcell *pgc)
{
    
}

void ioflow_gravity::wavegen_2D_precalc(lexer *p, fdm2D *b, ghostcell *pgc)
{
    
}

void ioflow_gravity::wavegen_2D_precalc_ini(lexer *p, ghostcell *pgc)
{
    
}

void ioflow_gravity::ini_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}

void ioflow_gravity::ini_ptf(lexer *p, fdm* a, ghostcell* pgc)
{
    
}

void ioflow_gravity::ini2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
    
}

void ioflow_gravity::veltimesave(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    
}

void ioflow_gravity::inflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, double *Fi, double *Uin,slice &Fifsf, slice &eta)
{

}

void ioflow_gravity::rkinflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, slice &frk, slice &f)
{
}

void ioflow_gravity::vrans_sed_update(lexer *p,fdm *a,ghostcell *pgc, vrans *pvrans)
{
    
}

void ioflow_gravity::ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{

}

void ioflow_gravity::discharge_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{

}

void ioflow_gravity::inflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{

}

void ioflow_gravity::rkinflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{

}

void ioflow_gravity::wavegen_precalc_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}

void ioflow_gravity::wavegen_precalc_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}

void ioflow_gravity::waterlevel2D(lexer *p, fdm2D *b, ghostcell* pgc, slice &eta)
{
    
}

void ioflow_gravity::fsfinflow_nhflow(lexer *p, fdm_nhf* d, ghostcell* pgc, slice &WL)
{

}