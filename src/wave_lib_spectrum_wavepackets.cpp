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
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void wave_lib_spectrum::wavepackets_parameters(lexer *p)
{
	double dd,w;
	double fac,Asum;
    double cmin,cmax;
    
    if(p->B94==0)
	wD=p->phimean;
    
	if(p->B94==1)
	wD=p->B94_wdt;
	
	p->wN=p->B86;
	
	p->Darray(Si,p->wN);
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
	
	ws=p->B87_1;
	we=p->B87_2;
	
	
	dd =(we-ws)/double(p->wN);	
	
	for(int n=0;n<p->wN;++n)
	dw[n] = dd;
	
	// Amplitudes Wavepackets
    if(p->B85==4)
	for(int n=0;n<p->wN;++n)
	{
	w = dw[n]*double(n) + ws;
	
	Ai[n] = (27.0*(w-ws)*pow(w-we,2.0)) / (4.0*1.5*pow(we-ws,3.0));
	}
    
    // Amplitudes Wavepackets Steep waves
    if(p->B85==5)
    {
        double maxA=0.0;
        double wP=0.0;
        for(int n=0;n<p->wN;++n)
        {
        w = dw[n]*double(n) + ws;
        
        
        Ai[n] = (27.0*(w-ws)*pow(w-we,2.0)) / (4.0*1.5*pow(we-ws,3.0));
        
            if(maxA>Ai[n])
            {
            maxA=Ai[n];
            wP=w;
            }

        }
        
        for(int n=0;n<p->wN;++n)
        {
        w = dw[n]*double(n) + ws;
        Ai[n] = (1.0 - pow((w-wP)/(we-wP),2.0))*(27.0*(w-ws)*pow(w-we,2.0)) / (4.0*1.5*pow(we-ws,3.0));
        }
    }
    
    // Amplitudes Gaussian Wavepackets
    if(p->B85==6)
	for(int n=0;n<p->wN;++n)
	{
	w = dw[n]*double(n) + ws;
	
	Ai[n] = (1.0/p->B83*sqrt(2.0*PI))*exp(-pow(w-p->wwp,2.0)/(2.0*p->B83*p->B83));
	}
	
	// Scale Amplitudes
	Asum=0.0;
	
	for(int n=0;n<p->wN;++n)
	Asum+=Ai[n];
	
    
	fac = p->wAs/Asum;

	for(int n=0;n<p->wN;++n)
	Ai[n]*=fac;	
    
	// Wave number
	w=ws;
	for(int n=0;n<p->wN;++n)
	{
	wL0 = (2.0*PI*9.81)/pow(w,2.0);
	k0 = (2.0*PI)/wL0;
	S0 = sqrt(k0*wD) * (1.0 + (k0*wD)/6.0 + (k0*k0*wD*wD)/30.0); 
	Li[n] = wL0*tanh(S0);
        
    for(int qn=0; qn<100; ++qn)
    Li[n] = wL0*tanh(2.0*PI*wD/Li[n]);

	ki[n] = 2.0*PI/Li[n];
	
	wi[n]=w;
	w+=dw[n];
	}
	
	// Focused Phases
    for(int n=0;n<p->wN;++n)
    {
        ei[n]  = ki[n]*p->B81_1 - wi[n]*p->B81_2;
    }

    
    // Group Velocities
    double duration,duration_tot;
    
    cmin = sqrt((9.81/ki[0])*tanh(ki[0]/wD));
    cmax = sqrt((9.81/ki[p->wN-1])*tanh(ki[p->wN-1]/wD));
    
    duration = 0.5*(p->B81_1/cmin + p->B81_1/cmax);
    
    if(p->B81_2<duration)
    duration_tot = p->B81_2; 
    
    if(p->B81_2>=duration)
    {
    duration_tot = p->B81_2; 
    
    p->wts = p->B81_2-duration;
    }
    
    p->wte = duration_tot;
    
    
    if(p->mpirank==0)
    cout<<"Wavepacket duration: "<<duration<<" wts: "<<p->wts<<" wte: "<<p->wte<<endl;
    
    print_amplitude_spectrum(p);
    
}


void wave_lib_spectrum::print_amplitude_spectrum(lexer *p)
{
	ofstream result;
	
	double xval=ws;

	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_Log-Wave",0777);
	
    if(p->mpirank==0)
    {
    // open file
	if(p->P14==0)
    result.open("REEF3D_amplitude-spectrum.dat");
	
	if(p->P14==1)
	result.open("./REEF3D_Log-Wave/REEF3D_amplitude-spectrum.dat");
	}	

	for(int n=0;n<p->wN;++n)
	{
	xval+=dw[n];
	result<<xval<<" "<<Ai[n]<<endl;
	}
	
	result.close();
}






