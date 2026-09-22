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

class weno_nug_func : public increment
{
public:
    weno_nug_func(lexer*);
    virtual ~weno_nug_func();

    void precalc_qf(lexer*);
    void precalc_cf(lexer*);
    void precalc_isf(lexer*);

    void ini(lexer*);

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
        const double is1x_psi = is1x + psi;
        const double is2x_psi = is2x + psi;
        const double is3x_psi = is3x + psi;

        w1x = cfx[IP][uf][0]/(epsilon + (is1x_psi*is1x_psi)*(cfx[IP][uf][0]/(is1x_psi*is1x_psi) + cfx[IP][uf][1]/(is2x_psi*is2x_psi) + cfx[IP][uf][2]/(is3x_psi*is3x_psi)));
        w2x = cfx[IP][uf][1]/(epsilon + (is2x_psi*is2x_psi)*(cfx[IP][uf][0]/(is1x_psi*is1x_psi) + cfx[IP][uf][1]/(is2x_psi*is2x_psi) + cfx[IP][uf][2]/(is3x_psi*is3x_psi)));
        w3x = cfx[IP][uf][2]/(epsilon + (is3x_psi*is3x_psi)*(cfx[IP][uf][0]/(is1x_psi*is1x_psi) + cfx[IP][uf][1]/(is2x_psi*is2x_psi) + cfx[IP][uf][2]/(is3x_psi*is3x_psi)));
    }
    inline void weight_max_x()
    {
        const double is1x_psi = is1x + psi;
        const double is2x_psi = is2x + psi;
        const double is3x_psi = is3x + psi;

        w1x = cfx[IP][uf][3]/(epsilon + (is1x_psi*is1x_psi)*(cfx[IP][uf][3]/(is1x_psi*is1x_psi) + cfx[IP][uf][4]/(is2x_psi*is2x_psi) + cfx[IP][uf][5]/(is3x_psi*is3x_psi)));
        w2x = cfx[IP][uf][4]/(epsilon + (is2x_psi*is2x_psi)*(cfx[IP][uf][3]/(is1x_psi*is1x_psi) + cfx[IP][uf][4]/(is2x_psi*is2x_psi) + cfx[IP][uf][5]/(is3x_psi*is3x_psi)));
        w3x = cfx[IP][uf][5]/(epsilon + (is3x_psi*is3x_psi)*(cfx[IP][uf][3]/(is1x_psi*is1x_psi) + cfx[IP][uf][4]/(is2x_psi*is2x_psi) + cfx[IP][uf][5]/(is3x_psi*is3x_psi)));
    }

    // y
    inline void weight_min_y()
    {
        const double is1y_psi = is1y + psi;
        const double is2y_psi = is2y + psi;
        const double is3y_psi = is3y + psi;

        w1y = cfy[JP][vf][0]/(epsilon + (is1y_psi*is1y_psi)*(cfy[JP][vf][0]/(is1y_psi*is1y_psi) + cfy[JP][vf][1]/(is2y_psi*is2y_psi) + cfy[JP][vf][2]/(is3y_psi*is3y_psi)));
        w2y = cfy[JP][vf][1]/(epsilon + (is2y_psi*is2y_psi)*(cfy[JP][vf][0]/(is1y_psi*is1y_psi) + cfy[JP][vf][1]/(is2y_psi*is2y_psi) + cfy[JP][vf][2]/(is3y_psi*is3y_psi)));
        w3y = cfy[JP][vf][2]/(epsilon + (is3y_psi*is3y_psi)*(cfy[JP][vf][0]/(is1y_psi*is1y_psi) + cfy[JP][vf][1]/(is2y_psi*is2y_psi) + cfy[JP][vf][2]/(is3y_psi*is3y_psi)));
    }
    inline void weight_max_y()
    {
        const double is1y_psi = is1y + psi;
        const double is2y_psi = is2y + psi;
        const double is3y_psi = is3y + psi;

        w1y = cfy[JP][vf][3]/(epsilon + (is1y_psi*is1y_psi)*(cfy[JP][vf][3]/(is1y_psi*is1y_psi) + cfy[JP][vf][4]/(is2y_psi*is2y_psi) + cfy[JP][vf][5]/(is3y_psi*is3y_psi)));
        w2y = cfy[JP][vf][4]/(epsilon + (is2y_psi*is2y_psi)*(cfy[JP][vf][3]/(is1y_psi*is1y_psi) + cfy[JP][vf][4]/(is2y_psi*is2y_psi) + cfy[JP][vf][5]/(is3y_psi*is3y_psi)));
        w3y = cfy[JP][vf][5]/(epsilon + (is3y_psi*is3y_psi)*(cfy[JP][vf][3]/(is1y_psi*is1y_psi) + cfy[JP][vf][4]/(is2y_psi*is2y_psi) + cfy[JP][vf][5]/(is3y_psi*is3y_psi)));
    }

    // z
    inline void weight_min_z()
    {
        const double is1z_psi = is1z + psi;
        const double is2z_psi = is2z + psi;
        const double is3z_psi = is3z + psi;

        w1z = cfz[KP][wf][0]/(epsilon + (is1z_psi*is1z_psi)*(cfz[KP][wf][0]/(is1z_psi*is1z_psi) + cfz[KP][wf][1]/(is2z_psi*is2z_psi) + cfz[KP][wf][2]/(is3z_psi*is3z_psi)));
        w2z = cfz[KP][wf][1]/(epsilon + (is2z_psi*is2z_psi)*(cfz[KP][wf][0]/(is1z_psi*is1z_psi) + cfz[KP][wf][1]/(is2z_psi*is2z_psi) + cfz[KP][wf][2]/(is3z_psi*is3z_psi)));
        w3z = cfz[KP][wf][2]/(epsilon + (is3z_psi*is3z_psi)*(cfz[KP][wf][0]/(is1z_psi*is1z_psi) + cfz[KP][wf][1]/(is2z_psi*is2z_psi) + cfz[KP][wf][2]/(is3z_psi*is3z_psi)));
    }
    inline void weight_max_z()
    {
        const double is1z_psi = is1z + psi;
        const double is2z_psi = is2z + psi;
        const double is3z_psi = is3z + psi;

        w1z = cfz[KP][wf][3]/(epsilon + (is1z_psi*is1z_psi)*(cfz[KP][wf][3]/(is1z_psi*is1z_psi) + cfz[KP][wf][4]/(is2z_psi*is2z_psi) + cfz[KP][wf][5]/(is3z_psi*is3z_psi)));
        w2z = cfz[KP][wf][4]/(epsilon + (is2z_psi*is2z_psi)*(cfz[KP][wf][3]/(is1z_psi*is1z_psi) + cfz[KP][wf][4]/(is2z_psi*is2z_psi) + cfz[KP][wf][5]/(is3z_psi*is3z_psi)));
        w3z = cfz[KP][wf][5]/(epsilon + (is3z_psi*is3z_psi)*(cfz[KP][wf][3]/(is1z_psi*is1z_psi) + cfz[KP][wf][4]/(is2z_psi*is2z_psi) + cfz[KP][wf][5]/(is3z_psi*is3z_psi)));
    }

    static double ****qfx,****qfy,****qfz;
    static double ***cfx,***cfy,***cfz;
    static double ****isfx,****isfy,****isfz;
    
	static int iniflag;
    
    
    
    double q1,q2,q3,q4,q5;

    const double epsilon,psi;
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

private:
    lexer *p;
};

#endif
