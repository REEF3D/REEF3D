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

#include"driver.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4.h"
#include"vec.h"
#include"vector"
#include <Eigen/Dense>

void driver::vec_test(lexer *p, fdm *a, ghostcell *pgc, field &f)
{	
	int qn,n;
	double t1,t2,t3,t4,t5,t6,t7;
	double val;
	
	vec vec(p);
	
	double *d;
	p->Darray(d,p->cellnum);
	
	vector<double> stdvector(p->cellnum,0.0);
	cout<<p->cellnum<<endl;
	
    Eigen::VectorXd eigenVector(p->cellnum);

    //std::array<double, 36000> stdarray;


    //--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	LOOP
	f(i,j,k)=0.0;
	
	for(qn=0; qn<1000; ++qn)
	LOOP
	val=f(i,j,k);
	
	t1 = pgc->timer() - starttime;


	//--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	f.V[IJK]=0.0;
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	val=f.V[IJK];
	
	t2 = pgc->timer() - starttime;
	

	//--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	vec.V[n]=0.0;
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	val=vec.V[n];
	
	t3 = pgc->timer() - starttime;
	

	//--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	d[n]=0.0;
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	val=d[n];
	
	t4 = pgc->timer() - starttime;
	

	//--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	stdvector[n]=0.0;
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	val=stdvector[n];
	
	t5 = pgc->timer() - starttime;	


	//--
    starttime = pgc->timer();
	/*
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	stdarray[n]=0.0;
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	val=stdarray[n];
	*/
	t6 = pgc->timer() - starttime;		
	
	
    //--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	eigenVector(n)=0.0;
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	val=eigenVector(n);
	
	t7 = pgc->timer() - starttime;	
	

	t1=pgc->globalmax(t1);
	t2=pgc->globalmax(t2);
	t3=pgc->globalmax(t3);
	t4=pgc->globalmax(t4);
	t5=pgc->globalmax(t5);
	t6=pgc->globalmax(t6);
	t7=pgc->globalmax(t7);
	
	if(p->mpirank==0)
	cout<<"t_field: "<<setprecision(15)<<t1<<"  t_field.V: "<<setprecision(9)<<t2<<"\n  t_vec: "<<setprecision(9)
    <<t3<<"  t_double: "<<setprecision(9)<<t4<<"    "<<"  t_stdvector: "<<setprecision(9)<<t5<<"  t_EigenVector: "<<setprecision(9)<<t7<<endl;


    //--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	LOOP
	f(i,j,k)=2.0*f(i,j,k);
	
	for(qn=0; qn<1000; ++qn)
	LOOP
	val=2.0*f(i,j,k);
	
	t1 = pgc->timer() - starttime;

    //--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	eigenVector=2.0*eigenVector;
	
	for(qn=0; qn<1000; ++qn)
	for(n=0; n<p->cellnum; ++n)
	val=2.0*eigenVector(n);
	
	t7 = pgc->timer() - starttime;	

	t1=pgc->globalmax(t1);
	t7=pgc->globalmax(t7);
	
    if(p->mpirank==0)
	cout<<"sum_field: "<<setprecision(9)<<t1<<"  sum_EigenVector: "<<setprecision(9)<<t7<<endl;
}

void driver::func_test(lexer *p, fdm *a, ghostcell *pgc, field &f)
{	
	int qn,n;
	double t1,t2;
	
	double val;
	nom=9.0;
	

	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	LOOP
	f(i,j,k)=(9.0 + nom*5.0)/nom;
	
	
	t1 = pgc->timer() - starttime;
	
	//--
	starttime = pgc->timer();
	
	for(qn=0; qn<1000; ++qn)
	LOOP
	f(i,j,k)=calc();
	
	t2 = pgc->timer() - starttime;
	
		
	t1=pgc->globalmax(t1);
	t2=pgc->globalmax(t2);

	
	if(p->mpirank==0)
	cout<<"t_inline: "<<setprecision(9)<<t1<<"  t_func: "<<setprecision(9)<<t2<<endl;
	
}

double driver::calc()
{
	
	val = (9.0 + nom*5.0)/nom;
	
	return val;
}
