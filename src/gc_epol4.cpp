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

#include"ghostcell.h"
#include"lexer.h"

// ghost-cell label of the cell-centred fields (kernel chosen in gcdistro4); cs: side of the face
int ghostcell::gceval4(lexer *p, int gcv, int bc, int cs)
{
    const bool wall = (bc==21||bc==22||bc==5);
    const bool sixdof = (bc==41||bc==42||bc==43);
    const bool patch = (bc==111||bc==112||bc==121||bc==122);
    const bool para_w = (cs==1||cs==2||cs==3||cs==4);        // side faces

    switch(gcv)
    {
    // level set
    case 51:
    case 52:
    case 53:
    case 54:
        if(wall||sixdof||bc==7||bc==8||bc==9||bc==111||bc==121||bc==211||bc==221) return gclabel_lsm;
        if(bc==3)                                   return 4;
        if((bc==1||bc==6) && (gcv==52||gcv==54))    return 4;
        if(bc==2 && (gcv==51||gcv==54||(gcv==52 && p->B77==1))) return 4;
        if(bc==6)                                   return gclabel_lsm_in;
        return 0;

    case 50:
        return 4;

    // pressure
    case 40:
        if(wall||bc==3||bc==111||bc==112||bc==211||bc==212) return gclabel_press;
        if(bc==6 && pressin_lable==0)               return gclabel_press;
        if(bc==7 && awa_lable==0)                   return gclabel_press;
        if(bc==1 && pressin_lable==0)               return gclabel_press_in;
        if(bc==2 && pressout_lable==0)              return gclabel_press;
        if(bc==9)                                   return 21;
        return 0;

    // density
    case 1:
        return 4;

    case 2:
        return (cs!=5 && bc!=5 && bc!=21) ? 4 : 0;

    // turbulence: k, epsilon / omega, eddy viscosity
    case 20:
        if(wall||sixdof||bc==9)                     return gclabel_k;
        if(bc==2 || (bc==3 && cs!=6))               return 4;
        if(bc==3 || bc==6 || bc==7 || bc==8)        return 5;
        return 0;

    case 30:
        if(wall||sixdof||bc==6||bc==7||bc==8||bc==9) return gclabel_e;
        if(bc==1||bc==2||bc==3)                     return 4;
        return 0;

    case 24:
        return 4;

    // omega (sigma coordinate)
    case 12:
        if(wall || (bc==7 && awa_lable==0))
        {
            if(para_w)                              return 4;
            if(cs==6)                               return 5;
            return 0;
        }
        if(bc==6)                                   return 4;
        if(bc==2 && gclabel_outflow==1)             return para_w ? 4 : ((cs==5||cs==6) ? 5 : 0);
        if(patch)                                   return 4;
        if(bc==3 && para_w)                         return 4;
        if(bc==3 && (cs==5||cs==6) && p->A10==5)    return 4;
        return 0;

    case 3:
        return (bc==2 && gclabel_outflow==1 && para_w) ? 4 : 0;

    // VOF
    case 70:
    case 71:
    case 72:
    case 73:
    case 74:
    case 75:
        if(wall||sixdof||bc==3||bc==6||bc==7||bc==8||bc==9) return gcv==75 ? 3 : 4;
        if(bc==1 && (gcv==70||gcv==72||gcv==74))    return 4;
        if(bc==2 && (gcv==70||gcv==71||gcv==74))    return 4;
        return 0;

    // velocities for the production term
    case 101:
    case 102:
    case 103:
        if(wall||bc==41)                            return 5;
        if(bc==1||bc==2||bc==6||bc==7||bc==8)       return 4;
        if(bc==3)
        {
            const bool orth = (gcv==101 && (cs==1||cs==4)) || (gcv==102 && (cs==2||cs==3)) || (gcv==103 && (cs==5||cs==6));
            const bool para = (cs>=1 && cs<=6) && !orth;
            if(para)                                return 4;
            if(orth)                                return 5;
        }
        return 0;

    // suspended sediment, heat
    case 60:
        return 4;

    case 80:
        if((p->H61==1 && cs==1) || (p->H62==1 && cs==2) || (p->H63==1 && cs==3)
        || (p->H64==1 && cs==4) || (p->H65==1 && cs==5) || (p->H66==1 && cs==6))
                                                    return 61;
        return 4;

    case 81:
        return 4;

    // potential flow initialisation
    case 49:
        if(wall||sixdof||bc==9||bc==3)              return 4;
        if(bc==1||bc==2||bc==6||bc==7||bc==8)       return 7;
        return 0;

    // potential waves
    case 250:
        if(bc==1||bc==2||bc==6||bc==7)              return 4;
        if((wall||sixdof||bc==8||bc==9) && cs!=5)   return 4;
        if(bc==3 && cs!=6)                          return 4;
        return 0;

    // NHFLOW
    case 540:
        if(bc==3 && cs==6)                          return 11;
        if(cs!=6 && (wall||bc==3||bc==6||bc==111||bc==112||bc==211||bc==212
           || (bc==2 && pressout_lable==0) || (bc==7 && awa_lable==0))) return 4;
        return 0;
    }

    return 0;
}


void ghostcell::gcdistro4(lexer *p, field &f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval4(p,gcv,bc,cs);

    switch(bc_label)
    {
    case 3: extend(p,f,dist,gcv,bc,cs); break;
    case 4: neumann(f,gcv,bc,cs); break;
    case 5: noslip(f,dist,gcv,bc,cs); break;
    case 7: potentialbc(p,f,bc,cs); break;
    case 11: nhpress(p,f,dist,gcv,bc,cs); break;
    case 21: atmosphere(p,f,gcv,bc,cs); break;
    case 61: heatbc(p,f,gcv,bc,cs); break;
    }
}


