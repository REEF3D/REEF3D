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

#include<sstream>
#include<string>
#include<vector>
#include<fstream>
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
    
    // direction of the 4th column of waverecon.dat, else along x
    for(int n=0;n<p->wN;++n)
    {
    beta[n] = recon_dir ? recon[n][3]*(PI/180.0) : 0.0;
    sinbeta[n] = recon_dir ? sin(beta[n]) : 0.0;
    cosbeta[n] = recon_dir ? cos(beta[n]) : 1.0;
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
    // waverecon.dat: one component per line, "A omega phase" or "A omega phase direction",
    // direction [deg] from the x axis of the generation frame (all lines the same form)
    std::ifstream file("waverecon.dat", std::ios_base::in);
	
	if(!file)
	cout<<endl<<("no 'waverecon.dat' file found")<<endl<<endl;
    
    std::vector<std::vector<double>> rows;
    std::string line;
    int ncol = 0;
    
    while(std::getline(file,line))
    {
        std::istringstream ls(line);
        std::vector<double> v;
        double x;
        
        while(ls>>x)
        v.push_back(x);
        
        if(v.size()<3)
        continue;
        
        if(ncol==0)
        ncol = v.size()>=4 ? 4 : 3;
        
        v.resize(4,0.0);
        if(ncol==3)
        v[3] = 0.0;
        
        rows.push_back(v);
    }
	
	file.close();

    wavenum = int(rows.size());
    recon_dir = ncol==4;
	
	p->Darray(recon,wavenum,4);
	
	for(int n=0; n<wavenum;++n)
	for(int q=0; q<4; ++q)
	recon[n][q] = rows[n][q];
    
    if(wavenum>0)
    {
	ts = recon[0][0];
	te = recon[wavenum-1][0];
    }
}
