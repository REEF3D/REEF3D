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
#include"mooring_Catenary.h"
#include"lexer.h"
#include"ghostcell.h"


void mooring_Catenary::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
	if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime+=p->P30;
		

        sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%08i-%06i.vtk",line,num);
		
		// Reconstruct line
		buildLine(p);
		

		// Print results
		ofstream result;
		result.open(name, ios::binary);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Mooring line " << line << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << H << " float" <<endl;

		for (int n = 0; n < H; ++n)
		{
			result<<x[n]<<" "<<y[n]<<" "<<z[n]<<endl;
		}
		
		result << "\nCELLS " << H-1 << " " << (H-1)*3 <<endl;	
		
		for(int n = 0; n < (H-1); ++n)
		{
			result<<"2 "<< n << " " << n+1 << endl;
		}
		
		result << "\nCELL_TYPES " << H-1 << endl;	
		
		for(int n = 0; n < (H-1); ++n)
		{
			result<<"3"<<endl;
		}	

		result<<"\nPOINT_DATA " << H <<endl;
		result<<"SCALARS Tension float 1 \nLOOKUP_TABLE default"<<endl;
		
		for(int n = 0; n < H; ++n)
		{
			result<<T[n]<<endl;
		}
		
		result.close();


		eTout<<p->simtime<<" \t "<<T[H-1]<<endl;	
	}
}


void mooring_Catenary::buildLine(lexer *p)
{
	double d_xy,dl, segLen, alpha;
    
    lms = L - FV/w;
    segLen = L/(H-1);
    alpha = atan(dy/dx);

    for (int cnt = 0; cnt < H; cnt++)
    {
        if (segLen*cnt <= lms)
        {
            z[cnt] = p->X311_zs[line];
            T[cnt] = fabs(FH);
                    
            d_xy = segLen*cnt*(1.0 + FH/EA);
                    
            if (dx > 0)
            {
                x[cnt] = p->X311_xs[line] + d_xy*cos(alpha);
            }
            else
            {
                x[cnt] = p->X311_xs[line] - d_xy*cos(alpha);
            }
                
            if (dy > 0)
            {
                y[cnt] = p->X311_ys[line] + d_xy*fabs(sin(alpha));
            }
            else
            {
                y[cnt] = p->X311_ys[line] - d_xy*fabs(sin(alpha));
            }
        }
        else
        {
            z[cnt] = p->X311_zs[line] + FH/w*(sqrt(1+(w*(segLen*cnt-lms)/FH)*(w*(segLen*cnt-lms)/FH)) - 1) + w*(segLen*cnt-lms)*(segLen*cnt-lms)/(2*EA);		
            
            T[cnt] = sqrt(FH*FH + (w*(segLen*cnt - lms))*(w*(segLen*cnt - lms)));
                    
            d_xy = FH/w*log(w*(segLen*cnt-lms)/FH + sqrt(1+(w*(segLen*cnt-lms)/FH)*(w*(segLen*cnt-lms)/FH))) + FH*segLen*cnt/EA;
                
            if (dx > 0)
            {
                x[cnt] = p->X311_xs[line] + lms*cos(alpha) + d_xy*cos(alpha);
            }
            else
            {
                x[cnt] = p->X311_xs[line] - lms*cos(alpha) - d_xy*cos(alpha);
            }
                    
            if (dy > 0)
            {
                y[cnt] = p->X311_ys[line] + lms*fabs(sin(alpha)) + d_xy*fabs(sin(alpha));
            }
            else
            {
                y[cnt] = p->X311_ys[line] - lms*fabs(sin(alpha)) - d_xy*fabs(sin(alpha));
            }					
        }
    }
}
