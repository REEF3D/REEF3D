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

class lexer;
class fdm;
class ghostcell;

#ifndef WAVE_LIB_SPECTRUM_H_
#define WAVE_LIB_SPECTRUM_H_

using namespace std;

class wave_lib_spectrum
{
public:
    wave_lib_spectrum();
	virtual ~wave_lib_spectrum();

    void irregular_parameters(lexer*);
    double wave_spectrum(lexer*, double);
	  void wavepackets_parameters(lexer*);

	  double PM(lexer*, double);
	  double JONSWAP(lexer*, double);
    double Goda_JONSWAP(lexer*, double);
    double TMA(lexer*, double);
	  double Torsethaugen(lexer*, double);
	  double spectrum_file(lexer*, double);

	  void spectrum_file_read(lexer*);

    void amplitudes_irregular(lexer*);
    void amplitudes_focused(lexer*);

    void phases_irregular(lexer*);
    void phases_focused(lexer*);

    void directional_spreading(lexer*);
    double spreading_function(lexer*,double,double);

    void print_spectrum(lexer*);
    void print_amplitude_spectrum(lexer*);
    void print_components(lexer*);
    void print_spreading(lexer*);

    double *Si,*Ai,*Li,*ki,*Ti,*wi,*ei,*ww,*cdf,*wee,*dee,*Sn;
    int NN, ND;
	  double wp,ws,we,*dw;

    // directional spreading
    double *beta,*cosbeta,*sinbeta,*Dn,*Dcdf,*Ddee,*betat,*beta_n,*Di,*Di_n;
    double dbeta;

    // recon
    void recon_read(lexer*, ghostcell*);
    void recon_parameters(lexer*, ghostcell*);
    int wavenum;
    double ts,te;
    double **recon;

private:
    double S,Sval,sigma,wD;
    double wL0,k0,S0;
    double s_f,G_0;
    double beta_J;

	int ptnum;
    int numcomp;


	double **spectrum;
};

#endif
