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

#include"ddweno_f_nug.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"
#include"ghostcell.h"
#include"vec.h"
#include"cpt.h"

ddweno_f_nug::ddweno_f_nug(lexer* pp):weno_nug_func(pp)
{
    p=pp;
}

ddweno_f_nug::~ddweno_f_nug()
{
}

double ddweno_f_nug::ddwenox(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
	grad=0.0;

	if(uw>=0.0)
	{
	iqmin(p,f);
	is_min_x();
	weight_min_x();

	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	iqmax(p,f);
	is_max_x();
	weight_max_x();
    
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}

	return grad;
}

double ddweno_f_nug::ddwenoy(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
	grad=0.0;

	if(uw>=0.0)
	{
	jqmin(p,f);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	jqmax(p,f);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}

	return grad;
}

double ddweno_f_nug::ddwenoz(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    wf=0;
    
    
	grad=0.0;

	if(uw>=0.0)
	{
	kqmin(p,f);
	is_min_z();
	weight_min_z();

	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}


	if(uw<0.0)
	{
	kqmax(p,f);
	is_max_z();
	weight_max_z();
    
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}
    
	return grad;
}



// Slice WENO with all intermediates in registers. The previous version kept
// q1..q5, is1..is3 and w1..w3 in class members and re-derived the 4-level
// coefficient pointers for every term, so the compiler could not common up the
// weight denominators (12 instead of 4 divisions per call). Arithmetic and
// operand order are unchanged -> bit-identical results.
namespace
{
inline double weno_slice_kernel(double q1, double q2, double q3, double q4, double q5,
                                double *const *isf, const double *cf, double *const *qf,
                                int r, double eps, double psi, bool xsign)
{
    const double dq12 = q1 - q2;
    const double dq23 = q2 - q3;
    const double dq32 = q3 - q2;
    const double dq34 = q3 - q4;
    const double dq43 = q4 - q3;
    const double dq54 = q5 - q4;
    
    const double is1 = isf[r+0][0]*dq54*dq54 + isf[r+0][1]*(dq54)*(dq34) + isf[r+0][2]*dq34*dq34;
    const double is2 = isf[r+1][0]*dq23*dq23 + isf[r+1][1]*(dq43)*(dq23) + isf[r+1][2]*dq43*dq43;
    const double is3 = isf[r+2][0]*dq12*dq12 + isf[r+2][1]*(dq32)*(dq12) + isf[r+2][2]*dq32*dq32;
    
    const double a1 = (is1 + psi)*(is1 + psi);
    const double a2 = (is2 + psi)*(is2 + psi);
    const double a3 = (is3 + psi)*(is3 + psi);
    
    const double c1 = cf[r+0], c2 = cf[r+1], c3 = cf[r+2];
    const double sum = c1/a1 + c2/a2 + c3/a3;
    
    const double w1 = c1/(eps + a1*sum);
    const double w2 = c2/(eps + a2*sum);
    const double w3 = c3/(eps + a3*sum);
    
    if(r==0)
    return w1*(q4 + qf[0][0]*(q3-q4) - qf[0][1]*(q5-q4))
         + w2*(q3 + qf[1][0]*(q4-q3) - qf[1][1]*(q2-q3))
         + w3*(q2 + qf[2][0]*(q1-q2) + qf[2][1]*(q3-q2));
    
    // uw<0: x and y used different signs on the second coefficient of
    // stencils 1 and 3 in the original code; both variants are kept as-is.
    if(xsign)
    return w1*(q4 + qf[3][0]*(q3-q4) + qf[3][1]*(q5-q4))
         + w2*(q3 + qf[4][0]*(q2-q3) - qf[4][1]*(q4-q3))
         + w3*(q2 + qf[5][0]*(q3-q2) - qf[5][1]*(q1-q2));
    
    return w1*(q4 + qf[3][0]*(q3-q4) - qf[3][1]*(q5-q4))
         + w2*(q3 + qf[4][0]*(q2-q3) - qf[4][1]*(q4-q3))
         + w3*(q2 + qf[5][0]*(q3-q2) + qf[5][1]*(q1-q2));
}
}

double ddweno_f_nug::dswenox(slice& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
    const int ip = IP;
    const double *const dx = p->DXP;
    const double fm3=f(i-3,j), fm2=f(i-2,j), fm1=f(i-1,j), f0=f(i,j), fp1=f(i+1,j), fp2=f(i+2,j), fp3=f(i+3,j);
    
	grad=0.0;

	if(uw>=0.0)
	grad = weno_slice_kernel((fm2-fm3)/dx[ip-3], (fm1-fm2)/dx[ip-2], (f0-fm1)/dx[ip-1], (fp1-f0)/dx[ip], (fp2-fp1)/dx[ip+1],
                             isfx[ip][uf], cfx[ip][uf], qfx[ip][uf], 0, epsilon, psi, true);

	if(uw<0.0)
	grad = weno_slice_kernel((fm1-fm2)/dx[ip-2], (f0-fm1)/dx[ip-1], (fp1-f0)/dx[ip], (fp2-fp1)/dx[ip+1], (fp3-fp2)/dx[ip+2],
                             isfx[ip][uf], cfx[ip][uf], qfx[ip][uf], 3, epsilon, psi, true);

	return grad;
}

double ddweno_f_nug::dswenoy(slice& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
    const int jp = JP;
    const double *const dy = p->DYP;
    const double fm3=f(i,j-3), fm2=f(i,j-2), fm1=f(i,j-1), f0=f(i,j), fp1=f(i,j+1), fp2=f(i,j+2), fp3=f(i,j+3);
    
	grad=0.0;

	if(uw>=0.0)
	grad = weno_slice_kernel((fm2-fm3)/dy[jp-3], (fm1-fm2)/dy[jp-2], (f0-fm1)/dy[jp-1], (fp1-f0)/dy[jp], (fp2-fp1)/dy[jp+1],
                             isfy[jp][vf], cfy[jp][vf], qfy[jp][vf], 0, epsilon, psi, false);

	if(uw<0.0)
	grad = weno_slice_kernel((fm1-fm2)/dy[jp-2], (f0-fm1)/dy[jp-1], (fp1-f0)/dy[jp], (fp2-fp1)/dy[jp+1], (fp3-fp2)/dy[jp+2],
                             isfy[jp][vf], cfy[jp][vf], qfy[jp][vf], 3, epsilon, psi, false);

	return grad;
}



void ddweno_f_nug::iqmin(lexer *p,field& f)
{	
	q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];
	q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void ddweno_f_nug::jqmin(lexer *p,field& f)
{
	q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];
	q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void ddweno_f_nug::kqmin(lexer *p,field& f)
{
	q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];
	q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void ddweno_f_nug::iqmax(lexer *p,field& f)
{
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
	q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void ddweno_f_nug::jqmax(lexer *p,field& f)
{
	q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
	q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void ddweno_f_nug::kqmax(lexer *p,field& f)
{
	q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
	q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}



void ddweno_f_nug::isqmin(lexer *p,slice& f)
{	
	q1 = (f(i-2,j)-f(i-3,j))/DX[IM3];
	q2 = (f(i-1,j)-f(i-2,j))/DX[IM2];
	q3 = (f(i,j)-f(i-1,j))/DX[IM1];
	q4 = (f(i+1,j)-f(i,j))/DX[IP];
	q5 = (f(i+2,j)-f(i+1,j))/DX[IP1];
}

void ddweno_f_nug::jsqmin(lexer *p,slice& f)
{
	q1 = (f(i,j-2)-f(i,j-3))/DY[JM3];
	q2 = (f(i,j-1)-f(i,j-2))/DY[JM2];
	q3 = (f(i,j)-f(i,j-1))/DY[JM1];
	q4 = (f(i,j+1)-f(i,j))/DY[JP];
	q5 = (f(i,j+2)-f(i,j+1))/DY[JP1];
}

void ddweno_f_nug::isqmax(lexer *p,slice& f)
{
    q1 = (f(i-1,j)-f(i-2,j))/DX[IM2];
	q2 = (f(i,j)-f(i-1,j))/DX[IM1];
	q3 = (f(i+1,j)-f(i,j))/DX[IP];
	q4 = (f(i+2,j)-f(i+1,j))/DX[IP1];
	q5 = (f(i+3,j)-f(i+2,j))/DX[IP2];
}

void ddweno_f_nug::jsqmax(lexer *p,slice& f)
{
	q1 = (f(i,j-1)-f(i,j-2))/DY[JM2];
	q2 = (f(i,j)-f(i,j-1))/DY[JM1];
	q3 = (f(i,j+1)-f(i,j))/DY[JP];
	q4 = (f(i,j+2)-f(i,j+1))/DY[JP1];
	q5 = (f(i,j+3)-f(i,j+2))/DY[JP2];
}
