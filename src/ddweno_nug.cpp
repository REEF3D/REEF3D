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

#include"ddweno_nug.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"vec.h"
#include"cpt.h"

ddweno_nug::ddweno_nug(lexer* pp):weno_nug_func(pp)
{
    p=pp;
}

ddweno_nug::~ddweno_nug()
{
}

double ddweno_nug::ddwenox(fdm* a, vec& b, double uw, int ipol, cpt &C)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
    int check=0;
    
	grad=0.0;

	if(uw>0.0)
	{
	iqmin(b,C);

	is_min_x();
	weight_min_x();


	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	iqmax(b,C);
	is_max_x();
	weight_max_x();
    
    
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}

	return grad;
}

double ddweno_nug::ddwenoy(fdm* a, vec& b, double uw, int ipol, cpt &C)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
	grad=0.0;

	if(uw>0.0)
	{
	jqmin(b,C);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	jqmax(b,C);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}

	return grad;
}

double ddweno_nug::ddwenoz(fdm* a, vec& b, double uw, int ipol, cpt &C)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    wf=0;

	grad=0.0;

	if(uw>0.0)
	{
	kqmin(b,C);
	is_min_z();
	weight_min_z();

	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}


	if(uw<0.0)
	{
	kqmax(b,C);
	is_max_z();
	weight_max_z();
    
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}
    
	return grad;
}
    
void ddweno_nug::iqmin(vec& f, cpt &C)
{
	q1 = (f.V[Im2_J_K] - f.V[Im3_J_K])/DX[IM3];
	q2 = (f.V[Im1_J_K] - f.V[Im2_J_K])/DX[IM2];
	q3 = (f.V[I_J_K]   - f.V[Im1_J_K])/DX[IM1];
	q4 = (f.V[Ip1_J_K] - f.V[I_J_K]  )/DX[IP];
	q5 = (f.V[Ip2_J_K] - f.V[Ip1_J_K])/DX[IP1];

}

void ddweno_nug::jqmin(vec& f, cpt &C)
{
	q1 = (f.V[I_Jm2_K] - f.V[I_Jm3_K])/DY[JM3];
	q2 = (f.V[I_Jm1_K] - f.V[I_Jm2_K])/DY[JM2];
	q3 = (f.V[I_J_K]   - f.V[I_Jm1_K])/DY[JM1];
	q4 = (f.V[I_Jp1_K] - f.V[I_J_K]  )/DY[JP];
	q5 = (f.V[I_Jp2_K] - f.V[I_Jp1_K])/DY[JP1];
}

void ddweno_nug::kqmin(vec& f, cpt &C)
{
	q1 = (f.V[I_J_Km2] - f.V[I_J_Km3])/DZ[KM3];
	q2 = (f.V[I_J_Km1] - f.V[I_J_Km2])/DZ[KM2];
	q3 = (f.V[I_J_K]   - f.V[I_J_Km1])/DZ[KM1];
	q4 = (f.V[I_J_Kp1] - f.V[I_J_K]  )/DZ[KP];
	q5 = (f.V[I_J_Kp2] - f.V[I_J_Kp1])/DZ[KP1];
}

void ddweno_nug::iqmax(vec& f, cpt &C)
{
	q1 = (f.V[Im1_J_K] - f.V[Im2_J_K])/DX[IM2];
    q2 = (f.V[I_J_K]   - f.V[Im1_J_K])/DX[IM1];
    q3 = (f.V[Ip1_J_K] - f.V[I_J_K]  )/DX[IP];
    q4 = (f.V[Ip2_J_K] - f.V[Ip1_J_K])/DX[IP1];
    q5 = (f.V[Ip3_J_K] - f.V[Ip2_J_K])/DX[IP2];
}

void ddweno_nug::jqmax(vec& f, cpt &C)
{
	q1 = (f.V[I_Jm1_K] - f.V[I_Jm2_K])/DY[JM2];
    q2 = (f.V[I_J_K]   - f.V[I_Jm1_K])/DY[JM1];
    q3 = (f.V[I_Jp1_K] - f.V[I_J_K]  )/DY[JP];
    q4 = (f.V[I_Jp2_K] - f.V[I_Jp1_K])/DY[JP1];
    q5 = (f.V[I_Jp3_K] - f.V[I_Jp2_K])/DY[JP2];
}

void ddweno_nug::kqmax(vec& f, cpt &C)
{
	q1 = (f.V[I_J_Km1] - f.V[I_J_Km2])/DZ[KM2];
    q2 = (f.V[I_J_K]   - f.V[I_J_Km1])/DZ[KM1];
    q3 = (f.V[I_J_Kp1] - f.V[I_J_K]  )/DZ[KP];
    q4 = (f.V[I_J_Kp2] - f.V[I_J_Kp1])/DZ[KP1];
    q5 = (f.V[I_J_Kp3] - f.V[I_J_Kp2])/DZ[KP2];
}
