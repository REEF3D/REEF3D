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

#ifndef WENO_NUG_FUNC_H_
#define WENO_NUG_FUNC_H_

#include"increment.h"

#include"lexer.h"
#include"field.h"
#include"slice.h"

using namespace std;

#include<cmath>

// Nonlinear WENO5 weights other than WENO-JS, shared by all weno_nug_func users and the
// NHFLOW register kernel. c1..c3 are the (non-uniform grid) ideal weights, is1..is3 the
// smoothness indicators in weno_nug_func ordering: is1 on {q3,q4,q5}, is3 on {q1,q2,q3},
// so tau5 = |is1 - is3| is the usual |beta_0 - beta_2|.
//   type 1: WENO-Z, p = 2 (Borges et al. 2008, JCP 227; Castro et al. 2011, JCP 230)
//   type 2: TENO5, C = 1, q = 6, cutoff ct (Fu, Hu & Adams 2016, JCP 305)
//   eps: regularisation of the smoothness indicators. It must not be tiny (1e-40 as in the
//   original papers): with a scale-free ratio tau/(is+eps), round-off level variations in an
//   otherwise uniform direction switch the weights at random, which can select downwind-biased
//   stencils and let transverse noise grow in 3D. With eps = psi (1e-6, as for WENO-JS)
//   near-uniform regions fall back to the linear upwind-biased scheme.
static inline __attribute__((always_inline)) void weno_weights_zt(const int type, const double ct, const double eps,
                                   const double c1, const double c2, const double c3,
                                   const double is1, const double is2, const double is3,
                                   double &w1, double &w2, double &w3)
{
    const double tau = std::fabs(is1 - is3);

    const double r1 = tau/(is1 + eps);
    const double r2 = tau/(is2 + eps);
    const double r3 = tau/(is3 + eps);

    if(type==1)
    {
        const double a1 = c1*(1.0 + r1*r1);
        const double a2 = c2*(1.0 + r2*r2);
        const double a3 = c3*(1.0 + r3*r3);
        const double inv = 1.0/(a1 + a2 + a3);

        w1 = a1*inv;
        w2 = a2*inv;
        w3 = a3*inv;
    }
    else
    {
        // (1+r)^6, capped so that the sum cannot overflow
        const double g1b = (r1 < 1.0e40) ? 1.0 + r1 : 1.0e40;
        const double g2b = (r2 < 1.0e40) ? 1.0 + r2 : 1.0e40;
        const double g3b = (r3 < 1.0e40) ? 1.0 + r3 : 1.0e40;
        const double g1 = (g1b*g1b)*(g1b*g1b)*(g1b*g1b);
        const double g2 = (g2b*g2b)*(g2b*g2b)*(g2b*g2b);
        const double g3 = (g3b*g3b)*(g3b*g3b)*(g3b*g3b);
        const double gct = ct*(g1 + g2 + g3);

        // sharp cutoff: a stencil is either kept with its ideal weight or discarded
        // (the largest chi is >= 1/3 > ct, so at least one stencil always survives)
        const double d1 = (g1 >= gct) ? c1 : 0.0;
        const double d2 = (g2 >= gct) ? c2 : 0.0;
        const double d3 = (g3 >= gct) ? c3 : 0.0;
        const double inv = 1.0/(d1 + d2 + d3);

        w1 = d1*inv;
        w2 = d2*inv;
        w3 = d3*inv;
    }
}

class weno_nug_func : public increment
{
public:
    weno_nug_func(lexer*);
    virtual ~weno_nug_func();

    void precalc_qf(lexer*);
    void precalc_cf(lexer*);
    void precalc_isf(lexer*);

    void ini(lexer*);

    // nonlinear weights: 0 WENO-JS (default, unchanged), 1 WENO-Z, 2 TENO5 with cutoff ct
    void set_weno_weights(int type, double ct)
    {
        wtype = (type==1 || type==2) ? type : 0;
        teno_ct = ct;
    }

    // Face divided differences dq(ii,j) = (f(ii+1,j)-f(ii,j))/DXP[ii] over every face the
    // WENO5 slice stencils of the interior cells touch, so that the q1..q5 of a cell are plain loads.
    void dsdiffx(slice&, slice&);
    void dsdiffy(slice&, slice&);

    // IS ----
    // x
    inline void is_min_x()
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1x = isfx[IP][uf][0][0]*dq54*dq54 + isfx[IP][uf][0][1]*(dq54)*(dq34) + isfx[IP][uf][0][2]*dq34*dq34;
        is2x = isfx[IP][uf][1][0]*dq23*dq23 + isfx[IP][uf][1][1]*(dq43)*(dq23) + isfx[IP][uf][1][2]*dq43*dq43;
        is3x = isfx[IP][uf][2][0]*dq12*dq12 + isfx[IP][uf][2][1]*(dq32)*(dq12) + isfx[IP][uf][2][2]*dq32*dq32;
    }
    inline void is_max_x()
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1x = isfx[IP][uf][3][0]*dq54*dq54 + isfx[IP][uf][3][1]*(dq54)*(dq34) + isfx[IP][uf][3][2]*dq34*dq34;
        is2x = isfx[IP][uf][4][0]*dq23*dq23 + isfx[IP][uf][4][1]*(dq43)*(dq23) + isfx[IP][uf][4][2]*dq43*dq43;
        is3x = isfx[IP][uf][5][0]*dq12*dq12 + isfx[IP][uf][5][1]*(dq32)*(dq12) + isfx[IP][uf][5][2]*dq32*dq32;
    }

    // y
    inline void is_min_y()
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1y = isfy[JP][vf][0][0]*dq54*dq54 + isfy[JP][vf][0][1]*(dq54)*(dq34) + isfy[JP][vf][0][2]*dq34*dq34;
        is2y = isfy[JP][vf][1][0]*dq23*dq23 + isfy[JP][vf][1][1]*(dq43)*(dq23) + isfy[JP][vf][1][2]*dq43*dq43;
        is3y = isfy[JP][vf][2][0]*dq12*dq12 + isfy[JP][vf][2][1]*(dq32)*(dq12) + isfy[JP][vf][2][2]*dq32*dq32;
    }
    inline void is_max_y()
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1y = isfy[JP][vf][3][0]*dq54*dq54 + isfy[JP][vf][3][1]*(dq54)*(dq34) + isfy[JP][vf][3][2]*dq34*dq34;
        is2y = isfy[JP][vf][4][0]*dq23*dq23 + isfy[JP][vf][4][1]*(dq43)*(dq23) + isfy[JP][vf][4][2]*dq43*dq43;
        is3y = isfy[JP][vf][5][0]*dq12*dq12 + isfy[JP][vf][5][1]*(dq32)*(dq12) + isfy[JP][vf][5][2]*dq32*dq32;
    }

    // z
    inline void is_min_z()
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1z = isfz[KP][wf][0][0]*dq54*dq54 + isfz[KP][wf][0][1]*(dq54)*(dq34) + isfz[KP][wf][0][2]*dq34*dq34;
        is2z = isfz[KP][wf][1][0]*dq23*dq23 + isfz[KP][wf][1][1]*(dq43)*(dq23) + isfz[KP][wf][1][2]*dq43*dq43;
        is3z = isfz[KP][wf][2][0]*dq12*dq12 + isfz[KP][wf][2][1]*(dq32)*(dq12) + isfz[KP][wf][2][2]*dq32*dq32;
    }
    inline void is_max_z()
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1z = isfz[KP][wf][3][0]*dq54*dq54 + isfz[KP][wf][3][1]*(dq54)*(dq34) + isfz[KP][wf][3][2]*dq34*dq34;
        is2z = isfz[KP][wf][4][0]*dq23*dq23 + isfz[KP][wf][4][1]*(dq43)*(dq23) + isfz[KP][wf][4][2]*dq43*dq43;
        is3z = isfz[KP][wf][5][0]*dq12*dq12 + isfz[KP][wf][5][1]*(dq32)*(dq12) + isfz[KP][wf][5][2]*dq32*dq32;
    }

    // Weights ----
    // x
    inline void weight_min_x()
    {
        // same arithmetic as before; the shared sum is formed once and the
        // coefficient row is loaded once (the member stores in between kept
        // the compiler from doing either)
        const double *const cf = cfx[IP][uf];
        const double c1 = cf[0], c2 = cf[1], c3 = cf[2];

        if(wtype!=0)
        {
            weno_weights_zt(wtype,teno_ct,psi,c1,c2,c3,is1x,is2x,is3x,w1x,w2x,w3x);
            return;
        }
        
        const double is1x_psi = is1x + psi;
        const double is2x_psi = is2x + psi;
        const double is3x_psi = is3x + psi;
        
        const double a1 = is1x_psi*is1x_psi;
        const double a2 = is2x_psi*is2x_psi;
        const double a3 = is3x_psi*is3x_psi;
        
        const double sum = c1/a1 + c2/a2 + c3/a3;

        w1x = c1/(epsilon + a1*sum);
        w2x = c2/(epsilon + a2*sum);
        w3x = c3/(epsilon + a3*sum);
    }
    inline void weight_max_x()
    {
        // same arithmetic as before; the shared sum is formed once and the
        // coefficient row is loaded once (the member stores in between kept
        // the compiler from doing either)
        const double *const cf = cfx[IP][uf];
        const double c1 = cf[3], c2 = cf[4], c3 = cf[5];

        if(wtype!=0)
        {
            weno_weights_zt(wtype,teno_ct,psi,c1,c2,c3,is1x,is2x,is3x,w1x,w2x,w3x);
            return;
        }
        
        const double is1x_psi = is1x + psi;
        const double is2x_psi = is2x + psi;
        const double is3x_psi = is3x + psi;
        
        const double a1 = is1x_psi*is1x_psi;
        const double a2 = is2x_psi*is2x_psi;
        const double a3 = is3x_psi*is3x_psi;
        
        const double sum = c1/a1 + c2/a2 + c3/a3;

        w1x = c1/(epsilon + a1*sum);
        w2x = c2/(epsilon + a2*sum);
        w3x = c3/(epsilon + a3*sum);
    }

    // y
    inline void weight_min_y()
    {
        // same arithmetic as before; the shared sum is formed once and the
        // coefficient row is loaded once (the member stores in between kept
        // the compiler from doing either)
        const double *const cf = cfy[JP][vf];
        const double c1 = cf[0], c2 = cf[1], c3 = cf[2];

        if(wtype!=0)
        {
            weno_weights_zt(wtype,teno_ct,psi,c1,c2,c3,is1y,is2y,is3y,w1y,w2y,w3y);
            return;
        }
        
        const double is1y_psi = is1y + psi;
        const double is2y_psi = is2y + psi;
        const double is3y_psi = is3y + psi;
        
        const double a1 = is1y_psi*is1y_psi;
        const double a2 = is2y_psi*is2y_psi;
        const double a3 = is3y_psi*is3y_psi;
        
        const double sum = c1/a1 + c2/a2 + c3/a3;

        w1y = c1/(epsilon + a1*sum);
        w2y = c2/(epsilon + a2*sum);
        w3y = c3/(epsilon + a3*sum);
    }
    inline void weight_max_y()
    {
        // same arithmetic as before; the shared sum is formed once and the
        // coefficient row is loaded once (the member stores in between kept
        // the compiler from doing either)
        const double *const cf = cfy[JP][vf];
        const double c1 = cf[3], c2 = cf[4], c3 = cf[5];

        if(wtype!=0)
        {
            weno_weights_zt(wtype,teno_ct,psi,c1,c2,c3,is1y,is2y,is3y,w1y,w2y,w3y);
            return;
        }
        
        const double is1y_psi = is1y + psi;
        const double is2y_psi = is2y + psi;
        const double is3y_psi = is3y + psi;
        
        const double a1 = is1y_psi*is1y_psi;
        const double a2 = is2y_psi*is2y_psi;
        const double a3 = is3y_psi*is3y_psi;
        
        const double sum = c1/a1 + c2/a2 + c3/a3;

        w1y = c1/(epsilon + a1*sum);
        w2y = c2/(epsilon + a2*sum);
        w3y = c3/(epsilon + a3*sum);
    }

    // z
    inline void weight_min_z()
    {
        // same arithmetic as before; the shared sum is formed once and the
        // coefficient row is loaded once (the member stores in between kept
        // the compiler from doing either)
        const double *const cf = cfz[KP][wf];
        const double c1 = cf[0], c2 = cf[1], c3 = cf[2];

        if(wtype!=0)
        {
            weno_weights_zt(wtype,teno_ct,psi,c1,c2,c3,is1z,is2z,is3z,w1z,w2z,w3z);
            return;
        }
        
        const double is1z_psi = is1z + psi;
        const double is2z_psi = is2z + psi;
        const double is3z_psi = is3z + psi;
        
        const double a1 = is1z_psi*is1z_psi;
        const double a2 = is2z_psi*is2z_psi;
        const double a3 = is3z_psi*is3z_psi;
        
        const double sum = c1/a1 + c2/a2 + c3/a3;

        w1z = c1/(epsilon + a1*sum);
        w2z = c2/(epsilon + a2*sum);
        w3z = c3/(epsilon + a3*sum);
    }
    inline void weight_max_z()
    {
        // same arithmetic as before; the shared sum is formed once and the
        // coefficient row is loaded once (the member stores in between kept
        // the compiler from doing either)
        const double *const cf = cfz[KP][wf];
        const double c1 = cf[3], c2 = cf[4], c3 = cf[5];

        if(wtype!=0)
        {
            weno_weights_zt(wtype,teno_ct,psi,c1,c2,c3,is1z,is2z,is3z,w1z,w2z,w3z);
            return;
        }
        
        const double is1z_psi = is1z + psi;
        const double is2z_psi = is2z + psi;
        const double is3z_psi = is3z + psi;
        
        const double a1 = is1z_psi*is1z_psi;
        const double a2 = is2z_psi*is2z_psi;
        const double a3 = is3z_psi*is3z_psi;
        
        const double sum = c1/a1 + c2/a2 + c3/a3;

        w1z = c1/(epsilon + a1*sum);
        w2z = c2/(epsilon + a2*sum);
        w3z = c3/(epsilon + a3*sum);
    }

    static double ****qfx,****qfy,****qfz;
    static double ***cfx,***cfy,***cfz;
    static double ****isfx,****isfy,****isfz;
    
	static int iniflag;
    
    
    
    double q1,q2,q3,q4,q5;

    const double epsilon,psi;
    int wtype;
    double teno_ct;
    double is1x,is2x,is3x;
    double is1y,is2y,is3y;
    double is1z,is2z,is3z;
    double w1x,w2x,w3x;
    double w1y,w2y,w3y;
    double w1z,w2z,w3z;

    int uf,vf,wf;
protected:
    inline void iqmin(field& f)
    {
        q1 = (f(i-2,j,k)-f(i-3,j,k))/p->DXP[IM3];
        q2 = (f(i-1,j,k)-f(i-2,j,k))/p->DXP[IM2];
        q3 = (f(i,j,k)-f(i-1,j,k))/p->DXP[IM1];
        q4 = (f(i+1,j,k)-f(i,j,k))/p->DXP[IP];
        q5 = (f(i+2,j,k)-f(i+1,j,k))/p->DXP[IP1];
    }
    inline void iqmax(field& f)
    {
        q1 = (f(i-1,j,k)-f(i-2,j,k))/p->DXP[IM2];
        q2 = (f(i,j,k)-f(i-1,j,k))/p->DXP[IM1];
        q3 = (f(i+1,j,k)-f(i,j,k))/p->DXP[IP];
        q4 = (f(i+2,j,k)-f(i+1,j,k))/p->DXP[IP1];
        q5 = (f(i+3,j,k)-f(i+2,j,k))/p->DXP[IP2];
    }

    inline void jqmin(field& f)
    {
        q1 = (f(i,j-2,k)-f(i,j-3,k))/p->DYP[JM3];
        q2 = (f(i,j-1,k)-f(i,j-2,k))/p->DYP[JM2];
        q3 = (f(i,j,k)  -f(i,j-1,k))/p->DYP[JM1];
        q4 = (f(i,j+1,k)-f(i,j,k)  )/p->DYP[JP];
        q5 = (f(i,j+2,k)-f(i,j+1,k))/p->DYP[JP1];
    }
    inline void jqmax(field& f)
    {
        q1 = (f(i,j-1,k)-f(i,j-2,k))/p->DYP[JM2];
        q2 = (f(i,j,k)-f(i,j-1,k))/p->DYP[JM1];
        q3 = (f(i,j+1,k)-f(i,j,k))/p->DYP[JP];
        q4 = (f(i,j+2,k)-f(i,j+1,k))/p->DYP[JP1];
        q5 = (f(i,j+3,k)-f(i,j+2,k))/p->DYP[JP2];
    }

    inline void kqmin(field& f)
    {
        q1 = (f(i,j,k-2)-f(i,j,k-3))/p->DZP[KM3];
        q2 = (f(i,j,k-1)-f(i,j,k-2))/p->DZP[KM2];
        q3 = (f(i,j,k)-f(i,j,k-1))/p->DZP[KM1];
        q4 = (f(i,j,k+1)-f(i,j,k))/p->DZP[KP];
        q5 = (f(i,j,k+2)-f(i,j,k+1))/p->DZP[KP1];
    }
    inline void kqmax(field& f)
    {
        q1 = (f(i,j,k-1)-f(i,j,k-2))/p->DZP[KM2];
        q2 = (f(i,j,k)-f(i,j,k-1))/p->DZP[KM1];
        q3 = (f(i,j,k+1)-f(i,j,k))/p->DZP[KP];
        q4 = (f(i,j,k+2)-f(i,j,k+1))/p->DZP[KP1];
        q5 = (f(i,j,k+3)-f(i,j,k+2))/p->DZP[KP2];
    }

    inline void isqmin(slice& f)
    {
        q1 = (f(i-2,j)-f(i-3,j))/p->DXP[IM3];
        q2 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];
        q3 = (f(i,j)-f(i-1,j))/p->DXP[IM1];
        q4 = (f(i+1,j)-f(i,j))/p->DXP[IP];
        q5 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];
    }
    inline void isqmax(slice& f)
    {
        q1 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];
        q2 = (f(i,j)-f(i-1,j))/p->DXP[IM1];
        q3 = (f(i+1,j)-f(i,j))/p->DXP[IP];
        q4 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];
        q5 = (f(i+3,j)-f(i+2,j))/p->DXP[IP2];
    }

    inline void jsqmin(slice& f)
    {
        q1 = (f(i,j-2)-f(i,j-3))/p->DYP[JM3];
        q2 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];
        q3 = (f(i,j)-f(i,j-1))/p->DYP[JM1];
        q4 = (f(i,j+1)-f(i,j))/p->DYP[JP];
        q5 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];
    }
    inline void jsqmax(slice& f)
    {
        q1 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];
        q2 = (f(i,j)-f(i,j-1))/p->DYP[JM1];
        q3 = (f(i,j+1)-f(i,j))/p->DYP[JP];
        q4 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];
        q5 = (f(i,j+3)-f(i,j+2))/p->DYP[JP2];
    }

    // WENO5 reconstruction from q1..q5
    inline double weno_min_x()
    {
        is_min_x();
        weight_min_x();

        return w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
             + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
    }
    inline double weno_max_x()
    {
        is_max_x();
        weight_max_x();

        return w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
             + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
    }

    inline double weno_min_y()
    {
        is_min_y();
        weight_min_y();

        return w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
             + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
    }
    inline double weno_max_y()
    {
        is_max_y();
        weight_max_y();

        return w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
             + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
    }

    inline double weno_min_z()
    {
        is_min_z();
        weight_min_z();

        return w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
             + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
             + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
    }
    inline double weno_max_z()
    {
        is_max_z();
        weight_max_z();

        return w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
             + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
             + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
    }

private:
    lexer *p;
};

#endif
