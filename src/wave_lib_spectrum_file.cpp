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

	// Debug first call
	static int call_count = 0;
	if(p->mpirank == 0 && call_count < 5)
	{
		cout<<"DEBUG spectrum_file_2d call "<<call_count<<": w="<<w<<" beta="<<beta_in<<endl;
		cout<<"  w_idx="<<w_idx<<" b_idx="<<b_idx<<endl;
		if(w_idx >= 0)
			cout<<"  freq range: ["<<freq_2d[w_idx]<<", "<<freq_2d[w_idx+1]<<"]"<<endl;
		if(b_idx >= 0)
			cout<<"  dir range: ["<<dir_2d[b_idx]<<", "<<dir_2d[b_idx+1]<<"]"<<endl;
		++call_count;
	}

	// If out of bounds, return 0
	if(w_idx == -1 || b_idx == -1)
	{
		if(p->mpirank == 0 && call_count <= 5)
			cout<<"  OUT OF BOUNDS - returning 0"<<endl;
		return 0.0;
	}

	// Get the 4 corner values for bilinear interpolation (using 1D indexing)
	S_w1_b1 = spectrum_2d_1d[w_idx * ptnum_dir_2d + b_idx];
	S_w1_b2 = spectrum_2d_1d[w_idx * ptnum_dir_2d + (b_idx+1)];
	S_w2_b1 = spectrum_2d_1d[(w_idx+1) * ptnum_dir_2d + b_idx];
	S_w2_b2 = spectrum_2d_1d[(w_idx+1) * ptnum_dir_2d + (b_idx+1)];

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

	ifstream file(name, ios_base::in);

	if(!file)
	{
		if(p->mpirank==0)
			cout<<endl<<("no 'spectrum-file-2d.dat' file found")<<endl<<endl;
	}

	// Read header line to determine number of directions
	if(file.is_open())
		getline(file, header);
	else
		header = "";
	istringstream iss(header);
	string temp;
	iss >> temp; // Skip "fq_dir"

	ptnum_dir_2d = 0;
	if(file.is_open())
	{
		while(iss >> val)
			++ptnum_dir_2d;
	}

	// Ensure minimum size even if file failed
	if(ptnum_dir_2d == 0)
		ptnum_dir_2d = 1;

	// Allocate direction array (use standard new[] to avoid Darray alignment issues)
	dir_2d = new double[ptnum_dir_2d];

	// Re-read header to get direction values
	if(file.is_open())
	{
		file.clear();
		file.seekg(0, ios::beg);
		getline(file, header);
		istringstream iss2(header);
		iss2 >> temp; // Skip "fq_dir"

		for(m = 0; m < ptnum_dir_2d; ++m)
			iss2 >> dir_2d[m];

	}
	else
	{
		// Initialize with dummy values
		for(m = 0; m < ptnum_dir_2d; ++m)
			dir_2d[m] = 0.0;
	}

	// Count number of frequency lines
	ptnum_freq_2d = 0;
	if(file.is_open())
	{
		while(!file.eof())
		{
			getline(file, header);
			if(header.length() > 0)
				++ptnum_freq_2d;
		}
		file.close();
	}

	// Ensure minimum size even if file failed
	if(ptnum_freq_2d == 0)
		ptnum_freq_2d = 1;

	// MPI: Broadcast dimensions from rank 0 to all other ranks
	int temp_vals[2];
	if(p->mpirank == 0)
	{
		temp_vals[0] = ptnum_freq_2d;
		temp_vals[1] = ptnum_dir_2d;
	}
	MPI_Bcast(temp_vals, 2, MPI_INT, 0, MPI_COMM_WORLD);
	ptnum_freq_2d = temp_vals[0];
	ptnum_dir_2d = temp_vals[1];

	// Allocate arrays (use standard new[] to avoid Darray alignment issues)
	freq_2d = new double[ptnum_freq_2d];
	spectrum_2d_1d = new double[ptnum_freq_2d * ptnum_dir_2d];

	// Allocate direction array on all ranks (was only allocated on rank 0 before)
	if(p->mpirank != 0)
		dir_2d = new double[ptnum_dir_2d];

	// Re-open and read data
	file.open(name, ios_base::in);

	if(file.is_open())
	{
		getline(file, header); // Skip header line

		n = 0;
		while(n < ptnum_freq_2d && !file.eof())
		{
			file >> freq_2d[n];
			for(m = 0; m < ptnum_dir_2d; ++m)
				file >> spectrum_2d_1d[n * ptnum_dir_2d + m];  // 1D indexing
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
	else
	{
		// Initialize with dummy values if file couldn't be read
		for(n = 0; n < ptnum_freq_2d; ++n)
		{
			freq_2d[n] = 0.0;
			for(m = 0; m < ptnum_dir_2d; ++m)
				spectrum_2d_1d[n * ptnum_dir_2d + m] = 0.0;  // 1D indexing
		}
	}

	// MPI: Broadcast all data arrays from rank 0 to all other ranks
	MPI_Bcast(freq_2d, ptnum_freq_2d, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(dir_2d, ptnum_dir_2d, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(spectrum_2d_1d, ptnum_freq_2d * ptnum_dir_2d, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
