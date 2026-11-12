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

#include"bcmom.h"
#include"lexer.h"
#include"bc_noflux.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

bcmom::bcmom(lexer* p):surftens(p),roughness(p),kappa(0.4)
{
	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

}

bcmom::~bcmom()
{
}

void bcmom::bcmom_start(fdm* a, lexer* p,ghostcell *pgc, turbulence *pturb,field& b,int gcval)
{
	int q;

	if(gcval==10 && p->B10!=0)
	{
	    QGC1LOOP
		if((p->gcb1[q][4]==5 || p->gcb1[q][4]==21 || p->gcb1[q][4]==22 || p->gcb1[q][4]==41 || p->gcb1[q][4]==42 || p->gcb1[q][4]==43) && p->gcb1[q][3]!=1 && p->gcb1[q][3]!=4 && !bc_periodic_face(p,p->gcb1[q][0],p->gcb1[q][1],p->gcb1[q][2],p->gcb1[q][3]))
		wall_law_u(a,p,pturb,b,p->gcb1[q][0], p->gcb1[q][1], p->gcb1[q][2], p->gcb1[q][3], p->gcb1[q][4], p->gcd1[q]);
        
        QGCDF1LOOP
        if(p->gcdf1[q][3]!=1 && p->gcdf1[q][3]!=4)   // tangential faces only, as in the gcb1 loop above
		wall_law_u(a,p,pturb,b,p->gcdf1[q][0], p->gcdf1[q][1], p->gcdf1[q][2], p->gcdf1[q][3], p->gcdf1[q][4],  0.5*p->DXM,
                   body_vel(a,p,b,0,p->gcdf1[q][0], p->gcdf1[q][1], p->gcdf1[q][2], p->gcdf1[q][3]));
	}

	if(gcval==11 && p->B10!=0 && p->j_dir==1)
	{
		QGC2LOOP
		if((p->gcb2[q][4]==5 || p->gcb2[q][4]==21 || p->gcb2[q][4]==22 || p->gcb2[q][4]==41 || p->gcb2[q][4]==42 || p->gcb2[q][4]==43) && p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3 && !bc_periodic_face(p,p->gcb2[q][0],p->gcb2[q][1],p->gcb2[q][2],p->gcb2[q][3]))
		wall_law_v(a,p,pturb,b,p->gcb2[q][0], p->gcb2[q][1], p->gcb2[q][2], p->gcb2[q][3], p->gcb2[q][4], p->gcd2[q]);
        
        QGCDF2LOOP
        if(p->gcdf2[q][3]!=2 && p->gcdf2[q][3]!=3)   // tangential faces only, as in the gcb2 loop above
		wall_law_v(a,p,pturb,b,p->gcdf2[q][0], p->gcdf2[q][1], p->gcdf2[q][2], p->gcdf2[q][3], p->gcdf2[q][4],  0.5*p->DXM,
                   body_vel(a,p,b,1,p->gcdf2[q][0], p->gcdf2[q][1], p->gcdf2[q][2], p->gcdf2[q][3]));
	}

	if(gcval==12 && p->B10!=0)
	{
		QGC3LOOP
		if((p->gcb3[q][4]==5 || p->gcb3[q][4]==21 || p->gcb3[q][4]==22 || p->gcb3[q][4]==41 || p->gcb3[q][4]==42 || p->gcb3[q][4]==43) && p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6 && !bc_periodic_face(p,p->gcb3[q][0],p->gcb3[q][1],p->gcb3[q][2],p->gcb3[q][3]))
		wall_law_w(a,p,pturb,b,p->gcb3[q][0], p->gcb3[q][1], p->gcb3[q][2], p->gcb3[q][3], p->gcb3[q][4], p->gcd3[q]);
        
        QGCDF3LOOP
        if(p->gcdf3[q][3]!=5 && p->gcdf3[q][3]!=6)   // tangential faces only, as in the gcb3 loop above
		wall_law_w(a,p,pturb,b,p->gcdf3[q][0], p->gcdf3[q][1], p->gcdf3[q][2], p->gcdf3[q][3], p->gcdf3[q][4],  0.5*p->DXM,
                   body_vel(a,p,b,2,p->gcdf3[q][0], p->gcdf3[q][1], p->gcdf3[q][2], p->gcdf3[q][3]));

	}
	surface_tension(a,p,a->phi,gcval);
}

// velocity of a 6DOF body (X 10) at the direct-forcing face cs of the velocity point (ii,jj,kk) of component c:
// the velocity of the point on the other side of the face if it lies inside a body (fb < 0 there; the direct
// forcing holds it at the body velocity); 0 at the faces of solids and topography. The wall law then acts on the
// velocity relative to the moving surface (it used the absolute velocity: a moving body braked the fluid next to
// it towards rest instead of dragging it along)
double bcmom::body_vel(fdm* a, lexer* p, field& b, int c, int ii, int jj, int kk, int cs)
{
    if(p->X10==0)
    return 0.0;

    if(cs==1) --ii;
    if(cs==4) ++ii;
    if(cs==3) --jj;
    if(cs==2) ++jj;
    if(cs==5) --kk;
    if(cs==6) ++kk;

    // level set of the bodies at the velocity point: the mean of its two cells
    const double phib = 0.5*(a->fb(ii,jj,kk) + a->fb(ii+(c==0),jj+(c==1),kk+(c==2)));

    if(phib<0.0)
    return b(ii,jj,kk);

    return 0.0;
}

void bcmom::wall_law_u(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist, double ub)
{
	i=ii;
	j=jj;
	k=kk;
    
    if(cs==2 || cs==3)
    deltaZ = p->DYN[JP];
    
    if(cs==5 || cs==6)
    deltaZ = p->DZN[KP];

    z0 = 0.5*deltaZ;
	
	ks=ks_val(p,a,ii,jj,kk,cs,bc);


		if(30.0*z0<ks)
		z0=ks/30.0;

		uplus = (1.0/kappa)*MAX(1.0,log(30.0*(z0/ks)));
    
	// b: velocity of the RK stage, relative to the wall (ub: a moving 6DOF body, 0 otherwise)
	const double ur = b(i,j,k) - ub;
	a->F(i,j,k) -= ((fabs(ur)*ur)/(uplus*uplus*deltaZ));
}

void bcmom::wall_law_v(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist, double ub)
{
	i=ii;
	j=jj;
	k=kk;
    
    if(cs==1 || cs==4)
    deltaZ = p->DXN[IP];
    
    if(cs==5 || cs==6)
    deltaZ = p->DZN[KP];
    
    z0 = 0.5*deltaZ;
    
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

		if(30.0*z0<ks)
		z0=ks/30.0;

		uplus = (1.0/kappa)*MAX(1.0,log(30.0*(z0/ks)));

	// b: velocity of the RK stage, relative to the wall (ub: a moving 6DOF body, 0 otherwise)
	const double ur = b(i,j,k) - ub;
	a->G(i,j,k) -= ((fabs(ur)*ur)/(uplus*uplus*deltaZ));
}

void bcmom::wall_law_w(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist, double ub)
{
	i=ii;
	j=jj;
	k=kk;
    
    if(cs==1 || cs==4)
    deltaZ = p->DXN[IP];
    
    if(cs==2 || cs==3)
    deltaZ = p->DYN[JP];
    
    z0 = 0.5*deltaZ;
	
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

		if(30.0*z0<ks)
		z0=ks/30.0;

		uplus = (1.0/kappa)*MAX(1.0,log(30.0*(z0/ks)));

	// b: velocity of the RK stage, relative to the wall (ub: a moving 6DOF body, 0 otherwise)
	const double ur = b(i,j,k) - ub;
	a->H(i,j,k) -= ((fabs(ur)*ur)/(uplus*uplus*deltaZ));
}




