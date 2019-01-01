/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"print_runup.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void print_runup::print(lexer *p, fdm* a, ghostcell *pgc)
{
	//if(p->P15==1)
    num = p->printcount;

    //if(p->P15==2)
    num = p->count;	
	
	
	if(p->mpirank==0)
	{
			// open file
			if(p->P14==0)
			{
			if(num<10)
			sprintf(name,"REEF3D-runup-fsf-00000%d.dat",num);

			if(num<100&&num>9)
			sprintf(name,"REEF3D-runup-fsf-0000%d.dat",num);

			if(num<1000&&num>99)
			sprintf(name,"REEF3D-runup-fsf-000%d.dat",num);

			if(num<10000&&num>999)
			sprintf(name,"REEF3D-runup-fsf-00%d.dat",num);

			if(num<100000&&num>9999)
			sprintf(name,"REEF3D-runup-fsf-0%d.dat",num);

			if(num>99999)
			sprintf(name,"REEF3D-runup-fsf-%d.dat",num);
			}
			
			if(p->P14==1)
			{
			if(num<10)
			sprintf(name,"./REEF3D_Runup-FSF/REEF3D-runup-fsf-00000%d.dat",num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_Runup-FSF/REEF3D-runup-fsf-0000%d.dat",num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_Runup-FSF/REEF3D-runup-fsf-000%d.dat",num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_Runup-FSF/REEF3D-runup-fsf-00%d.dat",num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_Runup-FSF/REEF3D-runup-fsf-0%d.dat",num);

			if(num>99999)
			sprintf(name,"./REEF3D_Runup-FSF/REEF3D-runup-fsf-%d.dat",num);
			}
			
			fsfout.open(name);
			
			// -----------------
			
			// open file
			if(p->P14==0)
			{
			if(num<10)
			sprintf(name,"REEF3D-runup-max-00000%d.dat",num);

			if(num<100&&num>9)
			sprintf(name,"REEF3D-runup-max-0000%d.dat",num);

			if(num<1000&&num>99)
			sprintf(name,"REEF3D-runup-max-000%d.dat",num);

			if(num<10000&&num>999)
			sprintf(name,"REEF3D-runup-max-00%d.dat",num);

			if(num<100000&&num>9999)
			sprintf(name,"REEF3D-runup-max-0%d.dat",num);

			if(num>99999)
			sprintf(name,"REEF3D-runup-max-%d.dat",num);
			}
			
			if(p->P14==1)
			{
			if(num<10)
			sprintf(name,"./REEF3D_Runup-MAX/REEF3D-runup-max-00000%d.dat",num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_Runup-MAX/REEF3D-runup-max-0000%d.dat",num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_Runup-MAX/REEF3D-runup-max-000%d.dat",num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_Runup-MAX/REEF3D-runup-max-00%d.dat",num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_Runup-MAX/REEF3D-runup-max-0%d.dat",num);

			if(num>99999)
			sprintf(name,"./REEF3D_Runup-MAX/REEF3D-runup-max-%d.dat",num);
			}
			
			maxout.open(name);
	}
	
	// instantaneous runup
	for(n=0;n<line_num;++n)
	if(cut_active[n]>0)
	fsfout<<cut_x[n]<<" \t "<<cut_y[n]<<" \t "<<cut_z[n]<<endl;
	
	fsfout.close();
	
	
	// max runup
	for(n=0;n<line_num;++n)
	if(runup_active[n]>0)
	maxout<<runup_x[n]<<" \t "<<runup_y[n]<<" \t "<<runup_z[n]<<endl;
	
	maxout.close();
}














