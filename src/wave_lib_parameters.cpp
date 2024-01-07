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

#include"wave_lib_parameters.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_parameters::wave_lib_parameters(lexer *p, ghostcell *pgc) : pshift(p->B120*(PI/180.0)), order(5)
{

  p->phiin=p->phimean;

	wtype=p->B92;

  p->wd = p->phimean;

  if(p->B94==0)
	wdt=p->phimean;

	if(p->B94==1)
	wdt=p->B94_wdt;

// Wave Length given
    if(p->B91==1 && (p->B92<30 || p->B92==70))
    {

     // define wave parameters
    wa = 0.5*p->B91_1;
    wH = p->B91_1;
    wL = p->B91_2;

    wk= (2.0*PI)/(wL>1.0e-20?wL:1.0e20);

    // define frequency
    if(wtype==1)
    ww= (2.0*PI*sqrt(9.81*wdt))/(wL);

    if(wtype==2)
    ww= sqrt(fabs(9.81*wk*tanh(wk*(wdt))));

    if(wtype==3)
    ww= sqrt(9.81*wk);

    if(wtype==4 || wtype>5)
    ww= sqrt(fabs(9.81*wk*tanh(wk*(wdt))));

  	wf = ww/(2.0*PI);
    wT = 1.0/wf;

        // 5th-order Stokes
        if(wtype==5)
        {
	    eps = 0.5*wk*wH;
	    S = 1.0/cosh(2*wk*wdt);
	    C = 1.0 - S;

	    c0 = sqrt(tanh(wk*wdt));

        c2 = (c0*(2.0 + 7.0*S*S)/(4.0*C*C));

        c4 = (c0*(4.0 + 32.0*S -116.0*S*S - 400.0*S*S*S - 71.0*pow(S,4.0) + 146.0*pow(S,5.0)))/(32.0*pow(C,5.0));

	    wT= (2.0*PI)/(sqrt(9.81*wk)*(c0 + eps*eps*c2 + eps*eps*eps*eps*c4));
	    wf = 1.0/wT;
	    ww = 2.0*PI*wf;

	    wC = ww/wk;
	    ubar = (c0 + eps*eps*c2 + eps*eps*eps*eps*c4)/sqrt(wk/9.81);
        
        //cout<<"C-Umean: "<<wC-ubar<<" C: "<<wC<<" Umean: "<<ubar<<endl;
	   }
       
    p->wT = wT;
    }

// Wave Period given
    if(p->B93==1 && (p->B92<30 || p->B92==70))
    {
		// define wave parameters
		wa = 0.5*p->B93_1;
		wH = p->B93_1;
		wT = p->B93_2;
        p->wT = wT;

		// define wave length
		if(wtype==1)
		wL= wT*sqrt(9.81*wdt);

		if(wtype==3)
		wL= (9.81/(2.0*PI))*wT*wT;


		if(wtype==2 || wtype==4 || wtype>5)
		{
		wL0 = (9.81/(2.0*PI))*wT*wT;
		k0 = (2.0*PI)/wL0;
		S0 = sqrt(k0*wdt) * (1.0 + (k0*wdt)/6.0 + (k0*k0*wdt*wdt)/30.0);

		wL = wL0*tanh(S0);

        for(int qn=0; qn<500; ++qn)
        wL = wL0*tanh(2.0*PI*wdt/wL);
		}


        // 5th-order Stokes
        if(wtype==5)
		{
		wL0 = (9.81/(2.0*PI))*wT*wT;
		k0 = (2.0*PI)/wL0;
		S0 = sqrt(k0*wdt) * (1.0 + (k0*wdt)/6.0 + (k0*k0*wdt*wdt)/30.0);

		wL = wL0*tanh(S0);

        for(int qn=0; qn<500; ++qn)
        wL = wL0*tanh(2.0*PI*wdt/wL);

        diff=10.0;
        int qn=0;
            do
            {
            wk_temp = (2.0*PI)/(wL>1.0e-20?wL:1.0e20); // wk: wavenumber k

            eps = 0.5*wk_temp*wH;
            S = 1.0/cosh(2*wk_temp*wdt);
            C = 1.0 - S;

            c0 = sqrt(tanh(wk_temp*wdt));
            c2 = (c0*(2.0 + 7.0*S*S)/(4.0*C*C));
            c4 = (c0*(4.0 + 32.0*S -116.0*S*S - 400.0*S*S*S - 71.0*pow(S,4.0) + 146.0*pow(S,5.0)))/(32.0*pow(C,5.0));

            wT_temp = (2.0*PI)/(sqrt(9.81*wk_temp)*(c0 + eps*eps*c2 + eps*eps*eps*eps*c4));

            diff=wT_temp-wT;

            //cout<<"wT_temp: "<<wT_temp<<" wT: "<<wT<<" wL: "<<wL<<" diff: "<<diff<<" qn: "<<qn<<endl;

            if(diff>0.0001)
            wL-=0.0005;

            if(diff<-0.0001)
            wL+=0.0005;
                ++qn;
            }while(fabs(diff)>0.0001 && qn<9000);

		}


		wf = 1.0/wT;
		ww = wf*2.0*PI;

		wk= (2.0*PI)/(wL>1.0e-20?wL:1.0e20);

        if(wtype==5)
        {
        wf = 1.0/wT;
        ww = 2.0*PI*wf;

        wC = ww/wk;
        ubar = (c0 + eps*eps*c2 + eps*eps*eps*eps*c4)/sqrt(wk/9.81);
        }
	   }



    p->wT = wT;
    p->wH = wH;
    p->wA = wa;
    p->wL = wL;
    p->wk = wk;
    p->ww = ww;

    if(p->B92>30 && p->B92!=70)
    {
       if(p->B91==1)
       {
         p->wHs = p->B91_1;
         p->wAs = 0.5*p->B91_1;
         p->wA = p->wAs;
         wL = p->B91_2;
         wk= (2.0*PI)/(wL>1.0e-20?wL:1.0e20);
         ww= sqrt(fabs(9.81*wk*tanh(wk*(wdt))));
         wf = ww/(2.0*PI);
         wT = 1.0/wf;
         p->wT = wT;
         p->wTp = wT;
         p->wwp = 2.0*PI/p->wTp;
       }

       if(p->B93==1)
       {
         p->wHs = p->B93_1;
         p->wAs = 0.5*p->B93_1;
         p->wA = p->wAs;
         p->wTp = p->B93_2;
         wT = p->B93_2;
         p->wT = wT;
         p->wwp = 2.0*PI/p->wTp;
       }
	   }

    p->Darray(factcos,order);

    // trig functions
    for(n=0; n<order; ++n)
    {
        factcos[n]=1.0;
        s=1.0;
        for(q=1;q<=(2*n);++q)
        {
        factcos[n] *= s;
        s+=1.0;
        }
    }

}

wave_lib_parameters::~wave_lib_parameters()
{
}

double wave_lib_parameters::sinfunc(double x)
{

    f = 0.0;


    for(n=0; n<order; ++n)
    {
        factorial=1;
        for(q=1;q<=(2*n+1);++q)
        factorial *= q;

    f += pow(-1.0,n)*pow(x,2*n+1)/(factorial);

    }

    return f;
}


double wave_lib_parameters::cosfunc(double x)
{

    f = 0.0;


    for(n=0; n<order; ++n)
    {
        r=1.0;
        for(int qr=1; qr<=n;++qr)
        r*=-1.0;


        s=1.0;
        for(int qr=1; qr<=2*n;++qr)
        s*=x;


    f += r*s/(factcos[n]);
    }

    return f;
}

double wave_lib_parameters::coshfunc(double x)
{
    f = 0.0;

    for(n=0; n<order; ++n)
    {

        s=1.0;
        for(int qr=1; qr<=2*n;++qr)
        s*=x;


    f += s/(factcos[n]);
    }

    return f;
}
