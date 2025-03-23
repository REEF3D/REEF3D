/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"ddweno_nug_sig.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"vec.h"
#include"cpt.h"

ddweno_nug_sig::ddweno_nug_sig(lexer* pp) : weno_nug_func(pp)
{
    p=pp;
}

ddweno_nug_sig::~ddweno_nug_sig()
{
}

double ddweno_nug_sig::ddwenox(double *F, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
    int check=0;
    
	grad=0.0;

	if(uw>0.0)
	{
	iqmin(F);

	is_min_x();
	weight_min_x();


	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	iqmax(F);
	is_max_x();
	weight_max_x();
    
    
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
    grad += 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((F[FIJKp1]-F[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));

	return grad;
}

double ddweno_nug_sig::ddwenoy(double *F, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
	grad=0.0;

	if(uw>0.0)
	{
	jqmin(F);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	jqmax(F);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}
    
    grad += 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((F[FIJKp1]-F[FIJKm1])/(p->DZN[KP]+p->DZN[KM1]));

	return grad;
}

double ddweno_nug_sig::ddwenoz(double *F, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    wf=0;

	grad=0.0;

	if(uw>0.0)
	{
	kqmin(F);
	is_min_z();
	weight_min_z();

	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}


	if(uw<0.0)
	{
	kqmax(F);
	is_max_z();
	weight_max_z();
    
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}
    
    grad *= p->sigz[IJ];
    
	return grad;
}
    
void ddweno_nug_sig::iqmin(double *F)
{
	q1 = (F[Im2JK] - F[Im3JK])/DX[IM3];
	q2 = (F[Im1JK] - F[Im2JK])/DX[IM2];
	q3 = (F[IJK]   - F[Im1JK])/DX[IM1];
	q4 = (F[Ip1JK] - F[IJK]  )/DX[IP];
	q5 = (F[Ip2JK] - F[Ip1JK])/DX[IP1];

}

void ddweno_nug_sig::jqmin(double *F)
{
	q1 = (F[IJm2K] - F[IJm3K])/DY[JM3];
	q2 = (F[IJm1K] - F[IJm2K])/DY[JM2];
	q3 = (F[IJK]   - F[IJm1K])/DY[JM1];
	q4 = (F[IJp1K] - F[IJK]  )/DY[JP];
	q5 = (F[IJp2K] - F[IJp1K])/DY[JP1];
}

void ddweno_nug_sig::kqmin(double *F)
{
	q1 = (F[IJKm2] - F[IJKm3])/DZ[KM3];
	q2 = (F[IJKm1] - F[IJKm2])/DZ[KM2];
	q3 = (F[IJK]   - F[IJKm1])/DZ[KM1];
	q4 = (F[IJKp1] - F[IJK]  )/DZ[KP];
	q5 = (F[IJKp2] - F[IJKp1])/DZ[KP1];
}

void ddweno_nug_sig::iqmax(double *F)
{
	q1 = (F[Im1JK] - F[Im2JK])/DX[IM2];
    q2 = (F[IJK]   - F[Im1JK])/DX[IM1];
    q3 = (F[Ip1JK] - F[IJK]  )/DX[IP];
    q4 = (F[Ip2JK] - F[Ip1JK])/DX[IP1];
    q5 = (F[Ip3JK] - F[Ip2JK])/DX[IP2];
}

void ddweno_nug_sig::jqmax(double *F)
{
	q1 = (F[IJm1K] - F[IJm2K])/DY[JM2];
    q2 = (F[IJK]   - F[IJm1K])/DY[JM1];
    q3 = (F[IJp1K] - F[IJK]  )/DY[JP];
    q4 = (F[IJp2K] - F[IJp1K])/DY[JP1];
    q5 = (F[IJp3K] - F[IJp2K])/DY[JP2];
}

void ddweno_nug_sig::kqmax(double *F)
{
	q1 = (F[IJKm1] - F[IJKm2])/DZ[KM2];
    q2 = (F[IJK]   - F[IJKm1])/DZ[KM1];
    q3 = (F[IJKp1] - F[IJK]  )/DZ[KP];
    q4 = (F[IJKp2] - F[IJKp1])/DZ[KP1];
    q5 = (F[IJKp3] - F[IJKp2])/DZ[KP2];
}
