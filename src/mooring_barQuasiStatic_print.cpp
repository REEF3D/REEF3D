/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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
#include"mooring_barQuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void mooring_barQuasiStatic::print(lexer *p,fdm *a, ghostcell *pgc)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
    // Check bottom and switch to catenary if necessary
    buildLine(p,a,pgc);
    checkBottom(p,a,pgc);
   

    // Print tension forces
    if (p->mpirank==0)
    {
        if(broken==true)
		eTout<<p->simtime<<" \t "<<0.0<<endl;
        
        if(broken==false)
		eTout<<p->simtime<<" \t "<<T[sigma+1]<<endl;
    }
    
    // Print line
	if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime+=p->P30;
		
		if(p->P14==1)
		{
			if(num<10)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-00000%i.vtk",line,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-0000%i.vtk",line,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-000%i.vtk",line,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-00%i.vtk",line,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-0%i.vtk",line,num);

			if(num>99999)
			sprintf(name,"./REEF3D_CFD_6DOF_Mooring/REEF3D-Mooring-%i-%i.vtk",line,num);
		}	


		ofstream result;
		result.open(name, ios::binary);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Mooring line " << line << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << sigma+2 << " float" <<endl;
		
		for(int n=0; n<sigma+2; ++n)
		{
			result<<x[n]<<" "<<y[n]<<" "<<z[n]<<endl;
		}
		
		result << "\nCELLS " << sigma+1 << " " << (sigma+1)*3 <<endl;	
		
		for(int n=0; n<sigma+1; ++n)
		{
			result<<"2 "<< n << " " << n+1 << endl;
		}
		
		result << "\nCELL_TYPES " << sigma+1 << endl;	
		
		for(int n=0; n<sigma+1; ++n)
		{
			result<<"3"<<endl;
		}	

		result<<"\nPOINT_DATA " << sigma+2 <<endl;
		result<<"SCALARS Tension float 1 \nLOOKUP_TABLE default"<<endl;
		
		for(int n=0; n<sigma+2; ++n)
		{
			result<<T[n]<<endl;
		}
		result.close();
	}
}


void mooring_barQuasiStatic::buildLine(lexer *p, fdm *a, ghostcell *pgc)
{
	x[0] = p->X311_xs[line];
	y[0] = p->X311_ys[line];
	z[0] = p->X311_zs[line];
	T[0] = fabs(A[0][0]);

	x[1] = x[0] + 0.5*l[1]*f[0][0];
	y[1] = y[0] + 0.5*l[1]*f[0][1];
	z[1] = z[0] + 0.5*l[1]*f[0][2];

	for (int cnt = 2; cnt < sigma+1; cnt++)
	{
		x[cnt] = x[cnt-1] + 0.5*(l[cnt] + l[cnt-1])*f[cnt-1][0];
		y[cnt] = y[cnt-1] + 0.5*(l[cnt] + l[cnt-1])*f[cnt-1][1];
		z[cnt] = z[cnt-1] + 0.5*(l[cnt] + l[cnt-1])*f[cnt-1][2];
			
		T[cnt-1] = fabs(A[cnt-1][cnt-1]);
	}		
	x[sigma+1] = x[sigma] + 0.5*(l[sigma])*f[sigma][0];
	y[sigma+1] = y[sigma] + 0.5*(l[sigma])*f[sigma][1];
	z[sigma+1] = z[sigma] + 0.5*(l[sigma])*f[sigma][2];

	T[sigma] = T[sigma+1] = fabs(A[sigma-1][sigma]);
}


void mooring_barQuasiStatic::checkBottom
(
    lexer *p, 
    fdm *a, 
    ghostcell *pgc
)
{
    for(int n=0; n<sigma+2; ++n)
    {
        if (z[n] < -0.01)
        {
            if (p->mpirank == 0 && broken == false) cout<<"Catenary solution"<<endl;
            pcatenary->getShape(p,a,pgc,x,y,z,T);
            break;
        }
    }
}
