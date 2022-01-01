/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Tobias Martin

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
--------------------------------------------------------------------*/

#include<sys/stat.h>
#include<iostream>
#include<fstream>
#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void fsi_strip::print_ini(lexer *p)
{
    // Ini print stl
	if(p->mpirank==0)
    {
        mkdir("./REEF3D_CFD_Beam_STL", 0777);
    }

    // Ini print parameter
    ofstream print;
    char str[1000];
    sprintf(str,"./REEF3D_CFD_Beam/REEF3D_Beam_position_%i.dat",nstrip);
    print.open(str);
	print<<"time \t x [m] \t y [m] \t z [m]"<<endl;
	print.close();
    
    sprintf(str,"./REEF3D_CFD_Beam/REEF3D_Beam_forces_%i.dat",nstrip);
    print.open(str);
	print<<"time \t Fx [N] \t Fy [N] \t Fz [N]"<<endl;
	print.close();

    // Ini triangulised strip
    tri_x = Matrix3Xd::Zero(3,2*Ne);
    tri_y = Matrix3Xd::Zero(3,2*Ne);
    tri_z = Matrix3Xd::Zero(3,2*Ne);

    // Ini print time
	printtime = 0.0;
    printcount_fsi = 0;
}


void fsi_strip::print_stl(lexer *p, fdm *a, ghostcell *pgc)
{
    int num = printcount_fsi;

    if(p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  || (p->simtime>printtime && p->P30>0.0)   || p->count==0))
    {
        printtime+=p->P30;
        
        char path[300];
        
        if(num<10)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-%i-00000%i.stl",nstrip,num);

        if(num<100&&num>9)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-%i-0000%i.stl",nstrip,num);

        if(num<1000&&num>99)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-%i-000%i.stl",nstrip,num);

        if(num<10000&&num>999)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-%i-00%i.stl",nstrip,num);

        if(num<100000&&num>9999)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-%i-0%i.stl",nstrip,num);

        if(num>99999)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-%i-%i.stl",nstrip,num);

        build_strip();

        ofstream result;
        result.open(path, ios::binary);
        
        result<<"solid"<<" "<<"ascii"<<endl;
        
        double zero=0.0;
        
        for(n=0; n<2*Ne; ++n)
        {
        result<<" facet normal "<<zero<<" "<<zero<<" "<<zero<<endl;
        result<<"  outer loop"<<endl;
        result<<"   vertex "<<tri_x(0,n)<<" "<<tri_y(0,n)<<" "<<tri_z(0,n)<<endl;
        result<<"   vertex "<<tri_x(1,n)<<" "<<tri_y(1,n)<<" "<<tri_z(1,n)<<endl;
        result<<"   vertex "<<tri_x(2,n)<<" "<<tri_y(2,n)<<" "<<tri_z(2,n)<<endl;
        result<<"  endloop"<<endl;
        result<<" endfacet"<<endl;
        }
        
        result<<"endsolid"<<endl;
        result.close();

        if(num<10)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-Lagrange-%i-00000%i.csv",nstrip,num);

        if(num<100&&num>9)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-Lagrange-%i-0000%i.csv",nstrip,num);

        if(num<1000&&num>99)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-Lagrange-%i-000%i.csv",nstrip,num);

        if(num<10000&&num>999)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-Lagrange-%i-00%i.csv",nstrip,num);

        if(num<100000&&num>9999)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-Lagrange-%i-0%i.csv",nstrip,num);

        if(num>99999)
        sprintf(path,"./REEF3D_CFD_Beam_STL/REEF3D-Beam-Lagrange-%i-%i.csv",nstrip,num);
        
        result.open(path, ios::binary);

	    result<<"X \t Y \t Z \t U \t V \t W"<<endl;
        for (int eI = 0; eI < Ne; eI++)
        {
            for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
            {
                const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);
                const Eigen::Vector3d& velI = lagrangeVel[eI].col(pI);
            
                result<<coordI(0)<<","<<coordI(1)<<","<<coordI(2)<<","<<velI(0)<<","<<velI(1)<<","<<velI(2)<<endl;
            }
        }
        result.close();

        printcount_fsi++;	
    }
}

void fsi_strip::print_parameter(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->mpirank == 0 && p->count%p->X19==0)
    {
        ofstream print;
        char str[1000];
        
        sprintf(str,"./REEF3D_CFD_Beam/REEF3D_Beam_position_%i.dat",nstrip);
        print.open(str, std::ofstream::out | std::ofstream::app);
        
        for (int eI = 0; eI < Ne+1; eI++)
        {
            print<<p->simtime<<" \t "<<x_el.col(eI).transpose()<<endl;
        }
        print.close();
        
        sprintf(str,"./REEF3D_CFD_Beam/REEF3D_Beam_forces_%i.dat",nstrip);
        print.open(str, std::ofstream::out | std::ofstream::app);
        Eigen::Vector3d sum = Eigen::Vector3d::Zero();
        for (int eI = 0; eI < Ne; eI++)
        {
            sum += F_el.col(eI);
        }
        print<<p->simtime<<" "<<sum.transpose()<<endl;
        print.close();
    }
}

void fsi_strip::build_strip()
{
    int ind = 0;

    for (int eI = 0; eI < Ne; eI++)
    {
        // Triangle 1
        tri_x.col(ind) << x_el(0,eI), x_el(0,eI+1), x_el(0,eI+1);
        tri_y.col(ind) << x_el(1,eI)-W_el/2, x_el(1,eI)-W_el/2, x_el(1,eI)+W_el/2;
        tri_z.col(ind) << x_el(2,eI), x_el(2,eI+1), x_el(2,eI+1);
       
        ind++;

        // Triangle 2
        tri_x.col(ind) << x_el(0,eI), x_el(0,eI+1), x_el(0,eI);
        tri_y.col(ind) << x_el(1,eI)+W_el/2, x_el(1,eI)+W_el/2, x_el(1,eI)-W_el/2;
        tri_z.col(ind) << x_el(2,eI), x_el(2,eI+1), x_el(2,eI);

        ind++;
    }
}
