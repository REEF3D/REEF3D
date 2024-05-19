/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"particle_corr.h"
#include"norm_vec.h"
#include"boundarycheck.h"
#include"field4.h"

using namespace std;

#ifndef PARTICLE_PLS_H_
#define PARTICLE_PLS_H_

class particle_pls : public particle_corr, public norm_vec, public boundarycheck
{
public:
	particle_pls(lexer*, fdm*, ghostcell*);
	virtual ~particle_pls();
	virtual void start(lexer*,fdm*,ghostcell*,ioflow*);
	void advect(lexer*,fdm*,ghostcell*,double**,int*,int);
	void seed(lexer*,fdm*,ghostcell*,double,double);
	void reseed(lexer*,fdm*,ghostcell*,double);
	void remove(lexer*,fdm*,ghostcell*);
	void random_delete(lexer*,fdm*,ghostcell*);
	void parcount(lexer*,fdm*,ghostcell*);
	void particlex(lexer*, fdm*, ghostcell*);
	void xupdate(lexer*,fdm*,ghostcell*);
	virtual void correct(lexer*, fdm*, ghostcell*,ioflow*);
	void parcorr(lexer*,fdm*,ioflow*,double,double*,int,int,int,int,int,int);
	virtual void picardmove(lexer*, fdm*, ghostcell*);

	virtual void setup(lexer*,fdm*,ghostcell*);
	void setradius(lexer*,fdm*);
	void posradius(lexer*,fdm*,int);
	void negradius(lexer*,fdm*,int);
	void inicorr(lexer*,fdm*,ghostcell*);
	void finalcorr(lexer*,fdm*,ghostcell*);
	int posseed(lexer*,fdm*,ghostcell*,double);
	int negseed(lexer*,fdm*,ghostcell*,double);
	

	double hside(fdm*);
	double phipol(lexer*,fdm*, double&,double&,double&);
	double upol(lexer*,fdm*, double&,double&,double&);
	double vpol(lexer*,fdm*, double&,double&,double&);
	double wpol(lexer*,fdm*, double&,double&,double&);
	double lint(field&,int&,int&,int&,double,double,double);
	double cint(double,double,double,double,double);
	double tricubic(lexer*,fdm*,field&,int&,int&,int&,double,double,double);
	void dgc_update(lexer*,fdm*,ghostcell*);
	void vel_setback(lexer*,fdm*,ghostcell*);
	void normal(fdm*, double&,double&,double&,double&);
	void normreg(fdm*, int,int,int);
	void allocate(lexer*,fdm*,ghostcell*);
	void print(lexer*,fdm*,ghostcell*);
	void print_ascii(lexer*,fdm*,ghostcell*);
	
	
	field4 phimax,phimin,phiold;
	field4 posnum, negnum;
	
	double **pos,**neg;
	double **posxs;
	double **negxs;
	double **posxr;
	double **negxr;
	int *pxs,*nxs;
	int *pxr,*nxr;
	int *posflag,*negflag;
	int *posmem,*negmem;
	int pcount, ncount;
	int n,nn,q,qq,count,check;
    double wa,wb,wc;
    double wx,wy,wz;
    double di,dj,dk,dnorm;
    double uvel,vvel,wvel;
    int posactive,negactive,maxparticle;
	int posactive_old, negactive_old, posbalance, negbalance;
	int gposactive_old, gnegactive_old, gposbalance, gnegbalance;
    int corrected,removed,xchange,reseeded;
    double val1,val2;

    double phix,phiy,phiz;
    double xs,ys,zs;
    double xp,yp,zp;
    double xc,yc,zc;
    double x1,x2,x3,x4,y1,y2;
    double xcell,ycell,zcell;
    double length,scalar,gamma;
    double coord1,coord2,coord3;
    double coord4,coord5,coord6;
    double u1,u2,v1,v2,w1,w2;
    double sp;
    int ii,jj,kk;

    double di0,di1,a0,a1,a2,a3,df;
    int i0,j0,k0,i3,j3,k3;
    int gnegactive,gposactive,gpcount,gncount,gcorrected,gremoved,greseeded,gxchange ;

	double H,Hval,nvec[3],phival,lambda,value,cosinus;
	const double zero,epsi,dx,rmin,rmax;
	const int pnum;
	const int ipolval;
	const int irand;
	const double drand;
	const double nu;
	
	double starttime;

	int maxpart, minpart;
	
	
	// PRINTVTU
	void print_vtu(lexer*,fdm*,ghostcell*,double**,int*,int,int);
	
	void pvtu_pos(fdm*,lexer*,ghostcell*);
	void pvtu_neg(fdm*,lexer*,ghostcell*);

    void header_pos(fdm*,lexer*,ghostcell*);
	void header_neg(fdm*,lexer*,ghostcell*);

    void piecename_pos(fdm*,lexer*,ghostcell*, int);
	void piecename_neg(fdm*,lexer*,ghostcell*, int);
	
	char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int gcval_phi;
    double printtime,printtime2;
    int printcount;
};

#endif

