/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"fluid_update_fsf_heat.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heat.h"

fluid_update_fsf_heat::fluid_update_fsf_heat(lexer *p, fdm* a, ghostcell* pgc, heat *&ppheat) : dx(p->DXM)
{
    gcval_ro=1;
	gcval_visc=1;

	visc_2 = p->W4;
	visc_1 = p->W2;
	ro_2 = p->W3;
	ro_1 = p->W1;
	alpha_air = p->H2;
	alpha_water = p->H1;

	material(p,a,pgc);
	
	pheat = ppheat;
}

fluid_update_fsf_heat::~fluid_update_fsf_heat()
{
}

void fluid_update_fsf_heat::start(lexer *p, fdm* a, ghostcell* pgc)
{
	double H=0.0;
	double temp;
	p->volume1=0.0;
	p->volume2=0.0;

    if(p->count>iter)
    iocheck=0;
	iter=p->count;
    
    if(p->j_dir==0)        
    epsi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    epsi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);

   //
	LOOP
	{
        temp = pheat->val(i,j,k);
	    
        if(p->H9==1)
        {
	    ro_1 = material_ipol(water_density,water_density_num, temp);
	    ro_2 = material_ipol(air_density,air_density_num, temp);

	    visc_1 = material_ipol(water_viscosity,water_viscosity_num, temp);
	    visc_2 = material_ipol(air_viscosity,air_viscosity_num, temp);
        }
        
        if(p->H9==2)
        {
	    ro_1 = material_ipol(air_density,air_density_num, temp);
        ro_2 = material_ipol(water_density,water_density_num, temp);

	    visc_1 = material_ipol(air_viscosity,air_viscosity_num, temp);
        visc_2 = material_ipol(water_viscosity,water_viscosity_num, temp);
        }

		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

		a->ro(i,j,k)=      ro_1*H +   ro_2*(1.0-H);
		a->visc(i,j,k)= visc_1*H + visc_2*(1.0-H);

		p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
		p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
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

double fluid_update_fsf_heat::material_ipol(double **pm, int num, double temp)
{
    double val=0.0;
    int n;

    for(n=0; n<num-1; ++n)
    {
        if(temp>pm[n][0] && temp<= pm[n+1][0])
        val = ((pm[n+1][1]-pm[n][1])/ (pm[n+1][0]-pm[n][0]))*(temp-pm[n][0]) + pm[n][1];
    }

    if(temp<=pm[0][0])
    val=pm[0][1];

    if(temp>pm[num-1][0])
    val=pm[num-1][1];

    return val;
}

void fluid_update_fsf_heat::material(lexer *p, fdm* a, ghostcell* pgc)
{
	p->Darray(water_density,20,2);
	p->Darray(air_density,20,2);
	p->Darray(water_viscosity,20,2);
	p->Darray(air_viscosity,20,2);

    water_density[0][0] = -30.0;
    water_density[0][1] = 983.854;

    water_density[1][0] = -20.0;
    water_density[1][1] = 993.547;

    water_density[2][0] = -10.0;
    water_density[2][1] = 998.117;

    water_density[3][0] = 0.0;
    water_density[3][1] = 999.8395;

    water_density[4][0] = 4.0;
    water_density[4][1] = 999.972;

    water_density[5][0] = 10.0;
    water_density[5][1] = 999.7026;

    water_density[6][0] = 15.0;
    water_density[6][1] = 999.1026;

    water_density[7][0] = 20.0;
    water_density[7][1] = 998.2071;

    water_density[8][0] = 22.0;
    water_density[8][1] = 997.7735;

    water_density[9][0] = 25.0;
    water_density[9][1] = 997.0479;

    water_density[10][0] = 30.0;
    water_density[10][1] = 995.6502;

    water_density[11][0] = 40.0;
    water_density[11][1] = 992.2;

    water_density[12][0] = 60.0;
    water_density[12][1] = 983.2;

    water_density[13][0] = 80.0;
    water_density[13][1] = 971.8;

    water_density[14][0] = 100.0;
    water_density[14][1] = 958.4;

    water_density_num = 15;

    //----------------------

    air_density[0][0] = -150.0;
    air_density[0][1] = 2.793;

    air_density[1][0] = -100.0;
    air_density[1][1] = 1.98;

    air_density[2][0] = -50.0;
    air_density[2][1] = 1.534;

    air_density[3][0] = 0.0;
    air_density[3][1] = 1.293;

    air_density[4][0] = 20.0;
    air_density[4][1] = 1.205;

    air_density[5][0] = 40.0;
    air_density[5][1] = 1.127;

    air_density[6][0] = 60.0;
    air_density[6][1] = 1.067;

    air_density[7][0] = 80.0;
    air_density[7][1] = 1.0;

    air_density[8][0] = 100.0;
    air_density[8][1] = 0.946;

    air_density[9][0] = 120.0;
    air_density[9][1] = 0.898;

    air_density[10][0] = 140.0;
    air_density[10][1] = 0.854;

    air_density[11][0] = 160.0;
    air_density[11][1] = 0.815;

    air_density[12][0] = 180.0;
    air_density[12][1] = 0.779;

    air_density[13][0] = 200.0;
    air_density[13][1] = 0.746;

    air_density[14][0] = 250.0;
    air_density[14][1] = 0.675;

    air_density[15][0] = 300.0;
    air_density[15][1] = 0.616;

    air_density[16][0] = 350.0;
    air_density[16][1] = 0.566;

    air_density[17][0] = 400.0;
    air_density[17][1] = 0.524;

    air_density_num = 18;

    //----------------------

    water_viscosity[0][0] = 0.0;
    water_viscosity[0][1] = 1.787e-6;

    water_viscosity[1][0] = 5.0;
    water_viscosity[1][1] = 1.519e-6;

    water_viscosity[2][0] = 10.0;
    water_viscosity[2][1] = 1.307e-6;

    water_viscosity[3][0] = 20.0;
    water_viscosity[3][1] = 1.002e-6;

    water_viscosity[4][0] = 30.0;
    water_viscosity[4][1] = 0.798e-6;

    water_viscosity[5][0] = 40.0;
    water_viscosity[5][1] = 0.653e-6;

    water_viscosity[6][0] = 50.0;
    water_viscosity[6][1] = 0.547e-6;

    water_viscosity[7][0] = 60.0;
    water_viscosity[7][1] = 0.467e-6;

    water_viscosity[8][0] = 70.0;
    water_viscosity[8][1] = 0.404e-6;

    water_viscosity[9][0] = 80.0;
    water_viscosity[9][1] = 0.355e-6;

    water_viscosity[10][0] = 90.0;
    water_viscosity[10][1] = 0.315e-6;

    water_viscosity[11][0] = 100.0;
    water_viscosity[11][1] = 0.282e-6;

    water_viscosity_num = 12;

    //----------------------

    air_viscosity[0][0] = -150.0;
    air_viscosity[0][1] = 3.08e-6;

    air_viscosity[1][0] = -100.0;
    air_viscosity[1][1] = 5.95e-6;

    air_viscosity[2][0] = -50.0;
    air_viscosity[2][1] = 9.55e-6;

    air_viscosity[3][0] = 0.0;
    air_viscosity[3][1] = 13.3e-6;

    air_viscosity[4][0] = 20.0;
    air_viscosity[4][1] = 15.11e-6;

    air_viscosity[5][0] = 40.0;
    air_viscosity[5][1] = 16.97e-6;

    air_viscosity[6][0] = 60.0;
    air_viscosity[6][1] = 18.9e-6;

    air_viscosity[7][0] = 80.0;
    air_viscosity[7][1] = 20.94e-6;

    air_viscosity[8][0] = 100.0;
    air_viscosity[8][1] = 23.06e-6;

    air_viscosity[9][0] = 120.0;
    air_viscosity[9][1] = 25.23e-6;

    air_viscosity[10][0] = 140.0;
    air_viscosity[10][1] = 27.55e-6;

    air_viscosity[11][0] = 160.0;
    air_viscosity[11][1] = 29.85e-6;

    air_viscosity[12][0] = 180.0;
    air_viscosity[12][1] = 32.29e-6;

    air_viscosity[13][0] = 200.0;
    air_viscosity[13][1] = 34.63e-6;

    air_viscosity[14][0] = 250.0;
    air_viscosity[14][1] = 41.17e-6;

    air_viscosity[15][0] = 300.0;
    air_viscosity[15][1] = 47.85e-6;

    air_viscosity[16][0] = 350.0;
    air_viscosity[16][1] = 55.05e-6;

    air_viscosity[17][0] = 400.0;
    air_viscosity[17][1] = 62.53e-6;

    air_viscosity_num = 18;

}

int fluid_update_fsf_heat::iocheck;
int fluid_update_fsf_heat::iter;
