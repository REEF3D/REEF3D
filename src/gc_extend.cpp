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
--------------------------------------------------------------------*/

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"
#include"cpt.h"

void ghostcell::extend(lexer *p,field& f,double dist,int gcv, int bc, int cs)
{
    weight=1.0;
    
    double dx;
    
    if(cs==1||cs==4)
    dx = p->DXP[IP];
    
    if(cs==2||cs==3)
    dx = p->DYP[JP];
    
    if(cs==5||cs==6)
    dx = p->DZP[KP];

//fill pos[]
    orderext=orderext2;

    if(bc_label==35)
    orderext=2;
	

	for(m=0;m<=orderext-3;m++)
	pos[m]=-dx*double(orderext-m-2);

	pos[orderext-2]=0.0;
	pos[orderext-1]=dx;

	for(m=0;m<margin;m++)
	x[m]=dx*double(m+2);
	
	

//fill y[]

	//Inflow
	if(cs==1 )
	for(m=0;m<orderext;m++)
	y[m]=f(i+orderext-m-1,j,k);

	//left
	if(cs==2)
	for(m=0;m<orderext;m++)
	y[m]=f(i,j-orderext+m+1,k);

	//right
	if(cs==3)
	for(m=0;m<orderext;m++)
	y[m]=f(i,j+orderext-m-1,k);

	//outlflow
	if(cs==4)
	for(m=0;m<orderext;m++)
	y[m]=f(i-orderext+m+1,j,k);

	//bottom
	if(cs==5)
	for(m=0;m<orderext;m++)
	y[m]=f(i,j,k+orderext-m-1);

	//top
	if(cs==6)
	for(m=0;m<orderext;m++)
	y[m]=f(i,j,k-orderext+m+1);
	
	
	for(q=0; q<margin; ++q)
	{
	    y[orderext+q]=0.0;

		for(m=0;m<orderext;m++)
		{
			weight=1.0;
			for(n=0;n<orderext;++n)
			{
			if(m!=n)
			weight*=(x[q]-pos[n])/(pos[m]-pos[n]+1.0e-20);
			}
		y[orderext+q]+=weight*y[m];
		}
	}

// write extrapolated ghost cell values into f()

	//inflow
	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=y[orderext+q];

	//left
	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=y[orderext+q];

	//right
	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=y[orderext+q];

	//outflow
	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=y[orderext+q];

	//bottom
	if(cs==5)
	{
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=y[orderext+q];
	
	//cout<<" LSM: "<<f(i,j,k+1)<<" "<<f(i,j,k)<<" "<<y[2]<<" "<<y[3]<<" "<<y[4]<<endl;
	}

	//top
	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=y[orderext+q];
}

void ghostcell::extendV(lexer *p,fdm *a,vec &f,double dist,int gcv, int bc, int cs)
{
    weight=1.0;
    
    double dx;
    
    if(cs==1||cs==4)
    dx = p->DXP[IP];
    
    if(cs==2||cs==3)
    dx = p->DYP[JP];
    
    if(cs==5||cs==6)
    dx = p->DZP[KP];

//fill pos[]
    orderext=2;

    if(bc_label==35)
    orderext=2;
	

	for(m=0;m<=orderext-3;m++)
	pos[m]=-dx*double(orderext-m-2);

	pos[orderext-2]=0.0;
	pos[orderext-1]=dx;

	for(m=0;m<margin;m++)
	x[m]=dx*double(m+2);
	
	
//fill y[]

	//Inflow
	if(cs==1 )
    {
    y[0]=f.V[Ip1_J_K_4]; //2-0-1
    y[1]=f.V[I_J_K_4]; //2-1-1
    }

	//left
	if(cs==2)
    {
    y[0]=f.V[I_Jm1_K_4]; //-2-0+1
    y[1]=f.V[I_J_K_4]; //-2-1+1
    }

	//right
	if(cs==3)
    {
	y[0]=f.V[I_Jp1_K_4]; //-2-0+1
    y[1]=f.V[I_J_K_4]; //-2-1+1
    }

	//outlflow
	if(cs==4)
	{
    y[0]=f.V[Im1_J_K_4]; //2-0-1
    y[1]=f.V[I_J_K_4]; //2-1-1
    }

	//bottom
	if(cs==5)
	{
    y[0]=f.V[I_J_Kp1_4]; //2-0-1
    y[1]=f.V[I_J_K_4]; //2-1-1
    }

	//top
	if(cs==6)
	{
    y[0]=f.V[I_J_Km1_4]; //2-0-1
    y[1]=f.V[I_J_K_4]; //2-1-1
    }
	
	
	for(q=0; q<margin; ++q)
	{
	    y[orderext+q]=0.0;

		for(m=0;m<orderext;m++)
		{
			weight=1.0;
			for(n=0;n<orderext;++n)
			{
			if(m!=n)
			weight*=(x[q]-pos[n])/(pos[m]-pos[n]+1.0e-20);
			}
		y[orderext+q]+=weight*y[m];
		}
	}

// write extrapolated ghost cell values into f()

	//inflow
	if(cs==1)
    {
    f.V[Im1_J_K_4]=y[orderext+0]; 
    f.V[Im2_J_K_4]=y[orderext+1];
    f.V[Im3_J_K_4]=y[orderext+2];
    }

	//left
	if(cs==2)
	{
	f.V[I_Jp1_K_4]=y[orderext+0]; 
    f.V[I_Jp2_K_4]=y[orderext+1];
    f.V[I_Jp3_K_4]=y[orderext+2];
    }

	//right
	if(cs==3)
	{
	f.V[I_Jm1_K_4]=y[orderext+0]; 
    f.V[I_Jm2_K_4]=y[orderext+1];
    f.V[I_Jm3_K_4]=y[orderext+2];
    }

	//outflow
	if(cs==4)
    {
	f.V[Ip1_J_K_4]=y[orderext+0]; 
    f.V[Ip2_J_K_4]=y[orderext+1];
    f.V[Ip3_J_K_4]=y[orderext+2];
    }

	//bottom
	if(cs==5)
	{
	f.V[I_J_Km1_4]=y[orderext+0]; 
    f.V[I_J_Km2_4]=y[orderext+1];
    f.V[I_J_Km3_4]=y[orderext+2];
    }

	//top
	if(cs==6)
	{
	f.V[I_J_Kp1_4]=y[orderext+0]; 
    f.V[I_J_Kp2_4]=y[orderext+1];
    f.V[I_J_Kp3_4]=y[orderext+2];
    }
}
