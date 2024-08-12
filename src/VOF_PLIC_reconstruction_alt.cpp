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
Author: Fabian Knooblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::reconstructPlane_alt_cube(fdm* a, lexer* p)
{
    double n_1, n_2, n_3, V0, V,d0, d;
    
    //get the normal vector
    calculateNormal_alt(a,p);
    
    //normalise normal vector (to be sure)
    double vecsum = sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
    nx(i,j,k)=nx(i,j,k)/vecsum;
    ny(i,j,k)=ny(i,j,k)/vecsum;
    nz(i,j,k)=nz(i,j,k)/vecsum;
    
    //sort vector components
    if(fabs(ny(i,j,k))>=fabs(nx(i,j,k)))
    {    
        n_1=fabs(nx(i,j,k));
        n_2=fabs(ny(i,j,k));
    }
    else
    {
        n_1=fabs(ny(i,j,k));
        n_2=fabs(nx(i,j,k));
    }
    
    if(fabs(nz(i,j,k))>=n_2)
    {
        n_3=fabs(nz(i,j,k));
    }
    else if(fabs(nz(i,j,k))<=n_1)
    {
        n_3=n_2;
        n_2=n_1;
        n_1=fabs(nz(i,j,k));
    }
    else
    {
        n_3=n_2;
        n_2=fabs(nz(i,j,k));
    }
    
    //reduced symmetry transformation
    V0=a->vof(i,j,k);
    V=0.5-fabs(V0-0.5);
    
    if(n_1+n_2<=2.0*V*n_3)   //case 5
    {
        d=V*n_3+(n_1+n_2)/2.0;
    }
    
    else if((3.0*n_2*(V*n_3+n_1-n_2)<=n_1*n_1) && (n_1*n_1<=6.0*V*n_2*n_3)) //case 2
    {
        d=n_1/2.0+sqrt(2*V*n_2*n_3-1.0/12.0*n_1*n_1);
    }
    
    else if(6.0*V*n_2*n_3<n_1*n_1) //case 1
    {
        d=cbrt(6.0*V*n_1*n_2*n_3);
    }
    
    else //cases 3 und 4
    {   
        double x3, y3, a3, b3, c3, f3 ,u3;
        x3=81.0*n_1*n_2*(n_1+n_2-2.0*V*n_3);
        
        if((23328.0*(n_1*n_2)*(n_1*n_2)*(n_1*n_2)-x3*x3)>=0.0)
        {
            y3=sqrt(23328.0*(n_1*n_2)*(n_1*n_2)*(n_1*n_2)-x3*x3);
            a3=cbrt(54.0)*n_1*n_2;
            b3=1.0/(cbrt(432.0));
            c3=n_1+n_2;
            u3=cbrt(x3*x3+y3*y3);
            f3= c3 -2.0*(a3+b3*cbrt(x3*x3+y3*y3))/(sqrt(u3)) *sin(PI/6.0-1.0/3.0*atan2(y3,x3));
        }
        else
            f3=-1000.0;
            
        if((n_2<=f3) && (f3<=min((n_1+n_2),n_3))) //case 3
            {
                d=f3;
            }
        else
        {
            double t4, x4, y4, a4, b4, c4, f4, u4;
            t4=9.0*(n_1+n_2+n_3)*(n_1+n_2+n_3)-18.0;
            x4=324.0*n_1*n_2*n_3*(1.0-2.0*V);
            if((4.0*t4*t4*t4-x4*x4)>=0.0)
            {
                y4=sqrt(4.0*t4*t4*t4-x4*x4);
                a4=1.0/(cbrt(864.0))*t4;
                b4=1.0/(cbrt(3456.0));
                c4=(n_1+n_2+n_3)/2.0;
                u4=cbrt(x4*x4+y4*y4);
                f4= c4 -2.0*(a4+b4*cbrt(x4*x4+y4*y4))/(sqrt(u4)) *sin(PI/6.0-1.0/3.0*atan2(y4,x4));
            }
            else
                f4=-1000.0;
            if(n_3 <= f4)   //case 4
            {
                d=f4;
            }
            else
            {
                d=0.5*p->DXN[IP];
                cout<<"Probleme beim Fall finden in der Inverse"<<endl;
            }
        }
    }
    
    //reverse reduced symmetry
    if((V0-0.5)>=0.0)
        d0=(n_1+n_2+n_3)/2.0-d;
    else
        d0=(-1,0)*((n_1+n_2+n_3)/2.0-d);
    
    alpha(i,j,k)=d0;
    
}
    

void VOF_PLIC::calculateNormal_alt(fdm* a, lexer* p)
{
    calcNormalLS(a, p);
}