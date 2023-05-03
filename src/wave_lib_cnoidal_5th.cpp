/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"wave_lib_cnoidal_5th.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_cnoidal_5th::wave_lib_cnoidal_5th(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: 5th-order cnoidal waves; ";
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wdt: "<<wdt<<" kd: "<<wdt*wk<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_cnoidal_5th::~wave_lib_cnoidal_5th()
{
}

double wave_lib_cnoidal_5th::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);
	
    return cosgamma*vel;
}

double wave_lib_cnoidal_5th::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);
	
    return singamma*vel;
}

double wave_lib_cnoidal_5th::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel;
	double yh;
	double sn,cn,dn;
	
	yh = (z+wdt)/wht; 
		
	teta = 2.0*Km*(x/wL - p->simtime/wT) + pshift;
	
	elliptic(p,teta,sn,cn,dn);
	
	vel =  wC + sqrt(9.81*wht)*(-1.0 + delta*(-0.5 + cn*cn) 
	    + pow(delta,2.0)*((-(19.0/40.0) + (3.0/2.0)*cn*cn - pow(cn,4.0)) + yh*yh*(-(3.0/2.0)*cn*cn + (9.0/4.0)*pow(cn,4.0)))
	
		+ pow(delta,3.0)*(-(55.0/112.0) + (71.0/40.0)*cn*cn - (27.0/10.0)*pow(cn,4.0) + (6.0/5.0)*pow(cn,6.0)
		
				+ pow(yh,2.0)*(-(9.0/4.0)*cn*cn + (75.0/8.0)*pow(cn,4.0) - (15.0/2.0)*pow(cn,6.0))
				
				+ pow(yh,4.0)*((3.0/8.0)*cn*cn - (45.0/16.0)*pow(cn,4.0) + (45.0/16.0)*pow(cn,6.0)))
			
		+ pow(delta,4.0)*(-(11813.0/22400.0) + (53327.0/42000.0)*cn*cn - (13109.0/3000.0)*pow(cn,4.0)
						+ (1763.0/375.0)*pow(cn,6.0) - (197.0/125.0)*pow(cn,8.0)
						
				+ pow(yh,2.0)*(-(213.0/80.0)*pow(cn,2.0) + (3231.0/160.0)*pow(cn,4.0) - (729.0/20.0)*pow(cn,6.0) + (189.0/10.0)*pow(cn,8.0))
				
				+ pow(yh,4.0)*((9.0/16.0)*cn*cn - (327.0/32.0)*pow(cn,4.0) + (915.0/32.0)*pow(cn,6.0) - (315.0/16.0)*pow(cn,8.0))
				
				+ pow(yh,6.0)*(-(3.0/80.0)*cn*cn + (189.0/160.0)*pow(cn,4.0) - (63.0/16.0)*pow(cn,6.0) + (189.0/64.0)*pow(cn,8.0)))
		
		+ pow(delta,5.0)*(-(57159.0/98560.0) - (144821.0/156800.0)*cn*cn - (1131733.0/294000.0)*pow(cn,4.0) + (757991.0/73500.0)*pow(cn,6.0)
						- (298481.0/36750.0)*pow(cn,8.0) + (13438.0/6125.0)*pow(cn,10.0)
						
				+ pow(yh,2.0)*(-(53327.0/28000.0)*cn*cn + (1628189.0/56000.0)*pow(cn,4.0) - (192481.0/2000.0)*pow(cn,6.0)
						+ (11187.0/100.0)*pow(cn,8.0) - (5319.0/125.0)*pow(cn,10.0))
						
				+ pow(yh,4.0)*((213.0/320.0)*cn*cn - (13563.0/640.0)*pow(cn,4.0) + (68643.0/640.0)*pow(cn,6.0)
						- (5481.0/32.0)*pow(cn,8.0) + (1701.0/20.0)*pow(cn,10.0))
						
				+ pow(yh,6.0)*(-(9.0/160.0)*cn*cn + (267.0/64.0)*pow(cn,4.0) - (987.0/32.0)*pow(cn,6.0) 
						+ (7875.0/128.0)*pow(cn,8.0) - (567.0/16.0)*pow(cn,10.0))
						
				+ pow(yh,8.0)*((9.0/4480.0)*cn*cn - (459.0/1792.0)*pow(cn,4.0) + (567.0/256.0)*pow(cn,6.0)
						- (1215.0/256.0)*pow(cn,8.0) + (729.0/256.0)*pow(cn,10.0))));		
	
    return vel;
}

double wave_lib_cnoidal_5th::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
	double yh;
	double sn,cn,dn;
	
	yh = (z+wdt)/wht;
	
	teta = 2.0*Km*(x/wL - p->simtime/wT) + pshift;

	elliptic(p,teta,sn,cn,dn);
	
	vel =  2.0*acn*dn*sn*cn*sqrt(9.81*wht)*( delta*yh 
	
	    + pow(delta,2.0)*(yh*((3.0/2.0) - 2.0*cn*cn) + pow(yh,3.0)*(-0.5 + 1.5*cn*cn))
		
		+ pow(delta,3.0)*(yh*((71.0/40.0) - (27.0/5.0)*cn*cn + (18.0/5.0)*pow(cn,4.0))
		
						+ pow(yh,3.0)*(-(3.0/4.0) + (25.0/4.0)*cn*cn - 7.5*pow(cn,4.0))
						
						+ pow(yh,5.0)*((3.0/40.0) - (9.0/8.0)*cn*cn + (27.0/16.0)*pow(cn,4.0)))
					
		+ pow(delta,4.0)*(yh*((53327.0/42000.0) - (26218.0/30000.0)*cn*cn + (1763.0/125.0)*pow(cn,4.0) - (788.0/125.0)*pow(cn,6.0))
		
						+ pow(yh,3.0)*(-(71.0/80.0) + (1077.0/80.0)*cn*cn - (729.0/20.0)*pow(cn,4.0) + (126.0/5.0)*pow(cn,6.0))
						
						+ pow(yh,5.0)*((9.0/80.0) - (327.0/80.0)*cn*cn + (549.0/32.0)*pow(cn,4.0) - (63.0/4.0)*pow(cn,6.0))
						
						+ pow(yh,7.0)*(-(30.0/560.0) + (27.0/80.0)*cn*cn - (27.0/16.0)*pow(cn,4.0) + (27.0/16.0)*pow(cn,6.0)))
						
		+ pow(delta,5.0)*(yh*(-(144821.0/156800.0) - (2263466.0/294000.0)*cn*cn + (2273973.0/735000)*pow(cn,4.0)
							- (1193927.0/36750.0)*pow(cn,6.0) + (67190.0/1225.0)*pow(cn,8.0))
							
						+ pow(yh,3.0)*(-(53327.0/84000.0) + (1628189.0/84000.0)*cn*cn - (192481.0/2000)*pow(cn,4.0)
							+ (3729.0/25.0)*pow(cn,6.0) - (1773.0/25.0)*pow(cn,8.0))
							
						+ pow(yh,5.0)*((213.0/1600.0) - (27126.0/1600.0)*cn*cn + (205929/3200.0)*pow(cn,4.0)
							- (5481.0/40.0)*pow(cn,6.0) - (1773.0/25.0)*pow(cn,8.0))
							
						+ pow(yh,7.0)*(-(9.0/1120.0) + (267.0/224.0)*cn*cn -(423.0/32.0)*pow(cn,4.0) 
						    + (1125.0/32.0)*pow(cn,6.0) + (1701.0/20.0)*pow(cn,8.0))
							
						+ pow(yh,9.0)*((1.0/4480.0) - (51.0/896.0)*cn*cn + (189.0/256.0)*pow(cn,4.0)
							- (135.0/64.0)*pow(cn,6.0) + (405.0/25.0)*pow(cn,8.0))));
							
    return vel;
}

double wave_lib_cnoidal_5th::wave_eta(lexer *p, double x, double y)
{
    double eta;
	double sn,cn,dn;
	
	teta = 2.0*Km*(x/wL - p->simtime/wT) + pshift;
	
	elliptic(p,teta,sn,cn,dn);
	
	eta =    -wdt + wht + wht*(epsilon*cn*cn + pow(epsilon,2.0)*(-0.75*cn*cn + 0.75*pow(cn,4.0)) 
	
			+ pow(epsilon,3.0)*((5.0/8.0)*cn*cn - (151.0/80.0)*pow(cn,4.0) + (101.0/80.0)*pow(cn,6.0))
			
			+ pow(epsilon,4.0)*(-(8209.0/6000.0)*cn*cn + (11641.0/3000.0)*pow(cn,4.0) - (112393.0/24000.0)*pow(cn,6.0)  
								+ (17367.0/8000.0)*pow(cn,8.0))
			
			+ pow(epsilon,5.0)*((364671.0/196000.0)*cn*cn - (2920931.0/392000.0)*pow(cn,4.0) + (2001361.0/156800.0)*pow(cn,6.0)  
								- (17906339.0/1568000.0)*pow(cn,8.0) + (1331817.0/313600.0)*pow(cn,10.0)));
	
	return eta;	
}

double wave_lib_cnoidal_5th::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_cnoidal_5th::parameters(lexer *p, ghostcell *pgc)
{
    double diff=1.0;
	int qq,maxiter;
	double modulus_old=0.5;
	double Ur;
	modulus = 0.9;
	maxiter =5000;
	
	Ur = (wH*wL*wL)/pow(wdt,3.0);
	if(p->mpirank==0)
	cout<<"Ursell number: "<<Ur<<endl;
	
	qq=1;
	

	while(fabs(diff)>1.0e-12)
	{
		++qq;
		modulus_old = modulus;
		
		
		modulus = (3.0/16.0)* ((wH*wL*wL)/(wdt*wdt*wdt*modulus_old*pow(K_elliptic_5(modulus_old),2.0)));
		
		diff = modulus - modulus_old;
		
		modulus = modulus_old + 0.1*diff;
		
		modulus=MIN(0.999, modulus);
	
		if(qq>maxiter)
		break;
	}
	
	if(p->mpirank==0)
	cout<<"MODULUS: "<<modulus<<"    qq: "<<qq<<endl;	
	
	Km = K_elliptic_5(modulus);
	Em = E_elliptic_5(modulus);
	
	ell = Em/Km;
	
	// coefficients
	
	double md = wdt;
	
	wht =  wdt*(1.0 - (wH/md)*ell + pow(wH/md,2.0)*(ell/4.0) + pow(wH/md,3.0)*((1.0/25.0)*ell + 0.25*ell*ell)
	    + pow(wH/md,4.0)*((573.0/2000.0)*ell - (57.0/400.0)*pow(ell,2.0) + 0.25*pow(ell,3.0))
		+ pow(wH/md,5.0)*(-(302159.0/1470000.0)*ell + (1779.0/2000.0)*ell*ell - (123.0/400.0)*pow(ell,3.0) + 0.25*pow(ell,4.0))); 

	if(p->mpirank==0)	
	cout<<"WAVE TROUGH: "<<wht<<endl;
	
	epsilon = wH/wht;
	
	acn = sqrt((3.0*epsilon)/4.0) 
	      * (1.0 - (5.0*epsilon)/8.0 + (71.0*pow(epsilon,2.0))/128.0 - (100627.0*pow(epsilon,3.0))/179200.0 
		    + (16259737.0*pow(epsilon,4.0))/28672000.0);
	
	delta = (4.0/3.0)*acn*acn;
	
	wC = sqrt(9.81*wht)*(1.0 + (wH/wht)*(0.5-ell) + pow(wH/wht,2.0)*(-3.0/20.0 + (5.0/12.0)*ell) + pow(wH/wht,3.0)*(3.0/56.0 - (19.0/600.0)*ell)
					    + pow(wH/wht,4.0)*(-309.0/5600.0  + (3719.0/21000.0)*ell) + pow(wH/wht,5.0)*(12237.0/616000.0  - (997699.0/8820000.0)*ell));
						
						
	wR = 9.81*wht*(1.5 + 0.5*(wH/wht) - pow(wH/wht,2.0)*(1.0/40.0) - pow(wH/wht,3.0)*(3.0/140.0) - pow(wH/wht,4.0)*(3.0/175.0)
						- pow(wH/wht,5.0)*(2427.0/154000.0));
	if(p->mpirank==0)
	cout<<"wC: "<<wC<<" wR: "<<wR<<endl;
}

void wave_lib_cnoidal_5th::wave_prestep(lexer *p, ghostcell *pgc)
{
}
