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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"wave_lib_spectrum.h"
#include"lexer.h"
#include"ghostcell.h"

void wave_lib_spectrum::recon_parameters(lexer *p, ghostcell *pgc)
{
    if(p->B94==0)
	wD=p->phimean;
	
	if(p->B94==1)
	wD=p->B94_wdt;
    
    double wL0;
    
    p->wN=wavenum;
    
	
	p->Darray(wi,p->wN);
	p->Darray(dw,p->wN);
	p->Darray(Ai,p->wN);
	p->Darray(Li,p->wN);
	p->Darray(ki,p->wN);
	p->Darray(Ti,p->wN);
	p->Darray(ei,p->wN);
    p->Darray(beta,p->wN);
    p->Darray(cosbeta,p->wN);
    p->Darray(sinbeta,p->wN);
    
    for(int n=0;n<p->wN;++n)
    {
    beta[n]=0.0;
    sinbeta[n]=0.0;
    cosbeta[n]=1.0;
    }
    
    
    // fillvalues for Ai, wi, Li, ki and ei
    for(int n=0;n<p->wN;++n)
	{
    // fill 
    Ai[n]=recon[n][0];
	wi[n]=recon[n][1];
    ei[n]=recon[n][2];
    
    
	// ki
	wL0 = (2.0*PI*9.81)/pow(wi[n],2.0);
	k0 = (2.0*PI)/wL0;
	S0 = sqrt(k0*wD) * (1.0 + (k0*wD)/6.0 + (k0*k0*wD*wD)/30.0); 
	Li[n] = wL0*tanh(S0);
        
    for(int qn=0; qn<100; ++qn)
    Li[n] = wL0*tanh(2.0*PI*wD/Li[n]);
    
	ki[n] = 2.0*PI/Li[n];
	}
    
}

void wave_lib_spectrum::recon_read(lexer *p, ghostcell* pgc)
{
	char name[100];
	double val,val0,val1,val2,val0n,val1n,val2n;
	int count;
	
	sprintf(name,"waverecon.dat");

// open file------------
	ifstream file(name, ios_base::in);
	
	if(!file)
	cout<<endl<<("no 'waverecon.dat' file found")<<endl<<endl;
	
    val0n=val1n=val2n=0.0;
	count=0;
	while(!file.eof())
	{
    val0n=val0;
    val1n=val1;
    val2n=val2;
    
	file>>val0>>val1>>val2;
	++count;
	}
    
    if(val0==val0n && val1==val1n && val2==val2n)
    --count;
	
	file.close();

    wavenum=count;
	
	p->Darray(recon,wavenum,3);
	
	file.open ("waverecon.dat", ios_base::in);
	
	count=0;
	for(int n=0; n<wavenum;++n)
	{
	file>>val0>>val1>>val2;

	recon[n][0] = val0;
	recon[n][1] = val1;
    recon[n][2] = val2;
	}	
    
    
	
	ts = recon[0][0];
	te = recon[wavenum-1][0];
}
