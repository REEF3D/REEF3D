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

#ifndef NHFLOW_STRAIN_H_
#define NHFLOW_STRAIN_H_

#include"nhflow_gradient.h"

class ghostcell;

using namespace std;

class nhflow_strain : public nhflow_gradient
{

public:
	nhflow_strain (lexer*,fdm_nhf*);
	virtual ~nhflow_strain();

	double sij(lexer*,fdm_nhf*,int,int);
	double qij(lexer*,fdm_nhf*,int,int);
	double pk(lexer*,fdm_nhf*);
    double Sij2(lexer*,fdm_nhf*);
    double Qij2(lexer*,fdm_nhf*);
	void Pk_update(lexer*,fdm_nhf*,ghostcell*);
    void Pk_b_update(lexer*,fdm_nhf*,ghostcell*);
	void wallf_update(lexer*,fdm_nhf*,ghostcell*,int*);
	double strainterm(lexer*,fdm_nhf*);
    double strainterm_v(lexer*,fdm_nhf*);   // vertical shear sqrt((du/dz)^2 + (dv/dz)^2) (A 561 1)
    double strainterm(lexer*,double*,double*,double*);
	double rotationterm(lexer*,fdm_nhf*);
    double rotationterm(lexer*,double*,double*,double*);
	double magSqrSd(lexer*,fdm_nhf*);
    double magSqrSd(lexer*,double*,double*,double*);
	double strainplain(lexer*,fdm_nhf*);
    
    // HLES split (A 561): horizontal eddy viscosity d->EVH = [nu_t] + (C_s Delta_h)^2 |S_h| + nu_b (A 562)
    void hles_update(lexer*,fdm_nhf*,ghostcell*);
    
    // QCR (A 563, Spalart 2000): extra stress dtau_ij = -C_cr1 (O_ik tau_jk + O_jk tau_ik) of the k-omega/k-eps
    // Boussinesq stress tau_ij = nu_t (du_i/dx_j + du_j/dx_i), O_ik = (du_i/dx_k - du_k/dx_i)/|grad u|;
    // qcr_update stores it at the cell centres, qcr_fx/fy/fz give D d(dtau_ij)/dx_j at (i,j,k) for the UH, VH, WH
    // momentum sources (isource/jsource/ksource)
    void qcr_update(lexer*,fdm_nhf*,ghostcell*);
    double qcr_fx(lexer*,fdm_nhf*);
    double qcr_fy(lexer*,fdm_nhf*);
    double qcr_fz(lexer*,fdm_nhf*);
    
	double *PK,*PK0,*PK_b;
    double *TQ11,*TQ22,*TQ33,*TQ12,*TQ13,*TQ23;

private:
    double s11,s22,s33,s12,s13,s23;
    double r11,r22,r33,r12,r13,r23;
    double ss11,ss22,ss33,ss12,ss13,ss23;
    double rr11,rr22,rr33,rr12,rr13,rr23;
    double q11,q22,q33,q12,q13,q23;
	double pkterm,s,q,val;
	const double epsi;

};

#endif
