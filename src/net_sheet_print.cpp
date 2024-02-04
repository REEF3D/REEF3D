/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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
#include"net_sheet.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void net_sheet::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof-1;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
 
    // Print forces
    if (p->mpirank==0)
    {
        char str[1000];

        sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_Net_Forces_%i.dat",nNet);
        ofstream header_out;
        header_out.open(str, std::ofstream::out | std::ofstream::app);
		header_out<<p->simtime<<" \t "<<Fx<<" "<<Fy<<" "<<Fz<<endl;
        header_out.close();
    }  
    
    // Print probe points
    if (p->mpirank==0 && p->X324 > 0)
    {
        for (int pp = 0; pp < p->X324; pp++)
        {
            char str[1000];
            sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_Net_%i_Point_Probe_%i.dat",nNet,pp+1);
            ofstream header_out;
            header_out.open(str, std::ofstream::out | std::ofstream::app);
            header_out<<p->simtime<<" \t "<<x_(probeKnot(pp),0)<<" \t "<<x_(probeKnot(pp),1)<<" \t "<<x_(probeKnot(pp),2)<<endl;
            header_out.close();
        }
    }  
    
	
    if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime += p->P30;
		
        if(p->P14==1)
		{
			if(num<10)
			sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%i-00000%i.stl",nNet,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%i-0000%i.stl",nNet,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%i-000%i.stl",nNet,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%i-00%i.stl",nNet,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%i-0%i.stl",nNet,num);

			if(num>99999)
			sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-%i-%i.stl",nNet,num);
		}	

        // Save net as .stl
        ofstream result;
        result.open(name, ios::binary);

	    result<<"solid"<<" "<<"ascii"<<endl;

        double x0, x1, x2, y0, y1, y2, z0, z1, z2, nx, ny, nz, mag;

        for(int n = 0; n < tend; ++n)
        {
            x0 = tri_x[n][0];
            x1 = tri_x[n][1];
            x2 = tri_x[n][2];
            
            y0 = tri_y[n][0];
            y1 = tri_y[n][1];
            y2 = tri_y[n][2];
            
            z0 = tri_z[n][0];
            z1 = tri_z[n][1];
            z2 = tri_z[n][2];          
            
            nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
            ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
            nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
                
            mag = sqrt(nx*nx + ny*ny + nz*nz);
            
            nx /= mag;
            ny /= mag;
            nz /= mag;
            
            result<<" facet normal "<<nx<<" "<<ny<<" "<<nz<<endl;
            result<<"  outer loop"<<endl;
            result<<"   vertex "<<tri_x[n][0]<<" "<<tri_y[n][0]<<" "<<tri_z[n][0]<<endl;
            result<<"   vertex "<<tri_x[n][1]<<" "<<tri_y[n][1]<<" "<<tri_z[n][1]<<endl;
            result<<"   vertex "<<tri_x[n][2]<<" "<<tri_y[n][2]<<" "<<tri_z[n][2]<<endl;
            result<<"  endloop"<<endl;
            result<<" endfacet"<<endl;
        }

        result<<"endsolid"<<endl;

        result.close(); 

        //- Print Lagrangian points
        if(num<10)
        sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-Lagrange-00000%i.csv",num);

        if(num<100&&num>9)
        sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-Lagrange-0000%i.csv",num);

        if(num<1000&&num>99)
        sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-Lagrange-000%i.csv",num);

        if(num<10000&&num>999)
        sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-Lagrange-00%i.csv",num);

        if(num<100000&&num>9999)
        sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-Lagrange-0%i.csv",num);

        if(num>99999)
        sprintf(name,"./REEF3D_CFD_6DOF_Net/REEF3D-Net-Lagrange-%i.csv",num);
        result.open(name, ios::binary);
        for (int ii = 0; ii < nK; ii++)
        {
            result<<x_(ii,0)<<","<<x_(ii,1)<<","<<x_(ii,2)<<","<<forces_knot(ii,0)<<","<<forces_knot(ii,1)<<","<<forces_knot(ii,2)<<endl;
        }
        result.close();
	}
}
