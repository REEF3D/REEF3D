/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include<string>
#include<sstream>
using namespace std;

double wave_lib_spectrum::spectrum_file(lexer *p, double w)
{

	Sval = 0.0;

	for(int n=0;n<ptnum-1;++n)
	if(w>spectrum[n][0] && w<spectrum[n+1][0])
	Sval = ((spectrum[n+1][1]-spectrum[n][1])/(spectrum[n+1][0]-spectrum[n][0]))*(w-spectrum[n][0]) + spectrum[n][1];

	if(w<spectrum[0][0] || w>spectrum[ptnum-1][0])
	Sval = 0.0;

	return Sval;
}

void wave_lib_spectrum::spectrum_file_read(lexer *p)
{
  char name[100];
	double val,val0,val1,Sval;
	int count;


	double ts,te;
	int timecount;

	sprintf(name,"spectrum-file.dat");

// open file------------
	ifstream file(name, ios_base::in);

	if(!file)
	{
		cout<<endl<<("no 'spectrum-file.dat' file found")<<endl<<endl;

	}

	count=0;
	while(!file.eof())
	{
	file>>val0>>val1;
	++count;
	}

	file.close();


    ptnum=count;

	p->Darray(spectrum,ptnum,2);

	file.open ("spectrum-file.dat", ios_base::in);

	count=0;
	while(!file.eof())
	{

	file>>val0>>val1;


	spectrum[count][0] = val0;
	spectrum[count][1] = val1;
	++count;
	}

	ts = spectrum[0][0];
	te = spectrum[ptnum-1][0];

}

double wave_lib_spectrum::spectrum_file_2d(lexer *p, double w, double beta_in)

{

	Sval = 0.0;
	int n, m;
	double S_w1_b1, S_w1_b2, S_w2_b1, S_w2_b2;
	double w1, w2, b1, b2;
	int w_idx = -1, b_idx = -1;

	// Find frequency indices for interpolation
	for(n = 0; n < ptnum_freq_2d - 1; ++n)
	{
		if(w >= freq_2d[n] && w <= freq_2d[n+1])
		{
			w_idx = n;
			w1 = freq_2d[n];
			w2 = freq_2d[n+1];
			break;
		}
	}

	// Find direction indices for interpolation
	for(m = 0; m < ptnum_dir_2d - 1; ++m)
	{
		if(beta_in >= dir_2d[m] && beta_in <= dir_2d[m+1])
		{
			b_idx = m;
			b1 = dir_2d[m];
			b2 = dir_2d[m+1];
			break;
		}
	}

	// If out of bounds, return 0
	if(w_idx == -1 || b_idx == -1)
		return 0.0;

	// Get the 4 corner values for bilinear interpolation
	S_w1_b1 = spectrum_2d[w_idx][b_idx];
	S_w1_b2 = spectrum_2d[w_idx][b_idx+1];
	S_w2_b1 = spectrum_2d[w_idx+1][b_idx];
	S_w2_b2 = spectrum_2d[w_idx+1][b_idx+1];

	// Bilinear interpolation
	double dw = w2 - w1;
	double db = b2 - b1;

	if(dw > 0.0 && db > 0.0)
	{
		double tw = (w - w1) / dw;
		double tb = (beta_in - b1) / db;

		Sval = (1.0 - tw) * (1.0 - tb) * S_w1_b1 +
		       (1.0 - tw) * tb * S_w1_b2 +
		       tw * (1.0 - tb) * S_w2_b1 +
		       tw * tb * S_w2_b2;
	}
	else if(dw > 0.0)
	{
		double tw = (w - w1) / dw;
		Sval = (1.0 - tw) * S_w1_b1 + tw * S_w2_b1;
	}
	else if(db > 0.0)
	{
		double tb = (beta_in - b1) / db;
		Sval = (1.0 - tb) * S_w1_b1 + tb * S_w1_b2;
	}
	else
	{
		Sval = S_w1_b1;
	}

	return Sval;
}

 

void wave_lib_spectrum::spectrum_file_2d_read(lexer *p)

{
	char name[100];
	int n, m;
	double val;
	string header;
    
	sprintf(name,"spectrum-file-2d.dat");

	// open file
	ifstream file(name, ios_base::in);

	if(!file)
	{
		cout<<endl<<("no 'spectrum-file-2d.dat' file found")<<endl<<endl;
		return;
	}

	// Read header line to determine number of directions
	getline(file, header);
	istringstream iss(header);
	string temp;
	iss >> temp; // Skip "fq_dir"
	ptnum_dir_2d = 0;

	while(iss >> val)
		++ptnum_dir_2d;

	// Allocate direction array
	p->Darray(dir_2d, ptnum_dir_2d);

	// Re-read header to get direction values
	file.clear();
	file.seekg(0, ios::beg);
	getline(file, header);
	istringstream iss2(header);
	iss2 >> temp; // Skip "fq_dir"
 
	for(m = 0; m < ptnum_dir_2d; ++m)
		iss2 >> dir_2d[m];
        
	// Count number of frequency lines
	ptnum_freq_2d = 0;
	while(!file.eof())
	{
		getline(file, header);
		if(header.length() > 0)
			++ptnum_freq_2d;
	}
	file.close();

	// Allocate arrays
	p->Darray(freq_2d, ptnum_freq_2d);
	p->Darray(spectrum_2d, ptnum_freq_2d, ptnum_dir_2d);

	// Re-open and read data
	file.open(name, ios_base::in);
	getline(file, header); // Skip header line

	n = 0;
	while(n < ptnum_freq_2d && !file.eof())
	{
		file >> freq_2d[n];
		for(m = 0; m < ptnum_dir_2d; ++m)
			file >> spectrum_2d[n][m];
		++n;
	}

	file.close();
 
	if(p->mpirank == 0)
	{
		cout << "2D Spectrum file read successfully" << endl;
		cout << "Frequencies: " << ptnum_freq_2d << " Directions: " << ptnum_dir_2d << endl;
		cout << "Freq range: " << freq_2d[0] << " to " << freq_2d[ptnum_freq_2d-1] << endl;
		cout << "Dir range: " << dir_2d[0] << " to " << dir_2d[ptnum_dir_2d-1] << endl;
	}

}