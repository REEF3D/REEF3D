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

#include<iostream>
#include<string>

using namespace std;

#ifndef RESIZE_CLASS_H_
#define RESIZE_CLASS_H_

class resize_class
{
public:

	resize_class();
	virtual ~resize_class();
	
	void rank(int);

	void Darray(double*&, int);
	void Darray(double**&, int, int);
	void Darray(double***&, int, int, int);
	void Darray(double****&, int, int, int, int);
	
	void Darray(double**&, int, int*);

	void Iarray(int*&, int);
	void Iarray(int**&, int, int);
	void Iarray(int***&, int, int, int);
	void Iarray(int****&, int, int, int, int);
	void Iarray(int*****&, int, int, int, int, int);
	
	void Iarray(int**&, int, int*);
	void Iarray(int***&, int, int*, int);

	void cvec(char*&, int);
	void cvec(char**&, int, int);
	void cvec(char***&, int, int, int);


	void del_Darray(double*&,int);
	void del_Darray(double**&,int, int);
	void del_Darray(double***&,int, int, int);
	void del_Darray(double****&,int,int, int, int);

	void del_Iarray(int*&,int);
	void del_Iarray(int**&,int, int);
	void del_Iarray(int***&,int, int, int);
	void del_Iarray(int****&,int, int, int, int);
	void del_Iarray(int*****&,int,int, int, int, int);
	
	void del_Iarray(int***&,int, int*, int);

	void del_cvec(char*&,int);
	void del_cvec(char**&,int, int);
	void del_cvec(char***&,int, int, int);
	
	void Dresize(double*&,int,int);
	void Dresize(double**&,int,int,int,int);
	void Dresize(double***&,int,int,int,int,int,int);
	void Dresize(double****&,int,int,int,int,int,int,int,int);
	
	void Iresize(int*&,int,int);
	void Iresize(int**&,int,int,int,int);
	void Iresize(int***&,int,int,int,int,int,int);
	void Iresize(int****&,int,int,int,int,int,int,int,int);
	
	void Iresize(int***&,int,int,int*,int*,int,int);
	
	int pararank;

};

#endif
