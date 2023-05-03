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

#include"wave_lib_spectrum.h"
#include"lexer.h"
#include"ghostcell.h"


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
