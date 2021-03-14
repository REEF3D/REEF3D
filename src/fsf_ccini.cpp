/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"fsf_vtp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fsf_vtp::ccini(lexer *p,fdm* a, ghostcell *pgc)
{
	double **ptcoor;
	int **ptijk;
	
	p->Darray(ptcoor,p->pointnum+p->ccptnum,3);
	p->Iarray(ptijk,p->pointnum,3);
	p->Darray(ccell,p->ccellnum,8,3);
	p->Iarray(ccijk,p->ccellnum,8,3);
	p->Iarray(ccnode,p->ccellnum);
	p->Iarray(ccid,p->ccellnum,8);
	p->Iarray(ccflag,p->ccellnum);
	p->Darray(lscc,p->ccellnum,8);
	p->Iarray(vertice_cc,p->ccellnum,8);
	
	for(n=0;n<p->ccellnum;++n)
	for(q=0;q<8;++q)
	ccid[n][q]=0;
	
	count=0;
    TPLOOP
	{
	ptcoor[count][0] = p->XN[IP1];
	ptcoor[count][1] = p->YN[JP1];
	ptcoor[count][2] = p->ZN[KP1];
	
	ptijk[count][0] = i;
	ptijk[count][1] = j;
	ptijk[count][2] = k;
	++count;
	}

	for(n=0;n<p->ccptnum;++n)
	{
    ptcoor[count][0] = p->ccpoint[n][0];
	ptcoor[count][1] = p->ccpoint[n][1];
	ptcoor[count][2] = p->ccpoint[n][2];
	++count;
	}
	
	// ccell
	for(n=0;n<p->ccellnum;++n)
	{
	ccnode[n]=a->ccedge[n];
	
    ccell[n][0][0]=ptcoor[abs(a->pvccnode[n][0])][0];
	ccell[n][0][1]=ptcoor[abs(a->pvccnode[n][0])][1];
	ccell[n][0][2]=ptcoor[abs(a->pvccnode[n][0])][2];
	
	if(abs(a->pvccnode[n][0])<p->pointnum)
	{
	ccijk[n][0][0]=ptijk[abs(a->pvccnode[n][0])][0];
	ccijk[n][0][1]=ptijk[abs(a->pvccnode[n][0])][1];
	ccijk[n][0][2]=ptijk[abs(a->pvccnode[n][0])][2];
	
	ccid[n][0]=1;
	}
	
	
	ccell[n][1][0]=ptcoor[abs(a->pvccnode[n][1])][0];
	ccell[n][1][1]=ptcoor[abs(a->pvccnode[n][1])][1];
	ccell[n][1][2]=ptcoor[abs(a->pvccnode[n][1])][2];
	
	if(abs(a->pvccnode[n][1])<p->pointnum)
	{
	ccijk[n][1][0]=ptijk[abs(a->pvccnode[n][1])][0];
	ccijk[n][1][1]=ptijk[abs(a->pvccnode[n][1])][1];
	ccijk[n][1][2]=ptijk[abs(a->pvccnode[n][1])][2];
	
	ccid[n][1]=1;
	}
	
	
	ccell[n][2][0]=ptcoor[abs(a->pvccnode[n][2])][0];
	ccell[n][2][1]=ptcoor[abs(a->pvccnode[n][2])][1];
	ccell[n][2][2]=ptcoor[abs(a->pvccnode[n][2])][2];
	
	if(abs(a->pvccnode[n][2])<p->pointnum)
	{
	ccijk[n][2][0]=ptijk[abs(a->pvccnode[n][2])][0];
	ccijk[n][2][1]=ptijk[abs(a->pvccnode[n][2])][1];
	ccijk[n][2][2]=ptijk[abs(a->pvccnode[n][2])][2];
	
	ccid[n][2]=1;
	}
	
	
	ccell[n][3][0]=ptcoor[abs(a->pvccnode[n][3])][0];
	ccell[n][3][1]=ptcoor[abs(a->pvccnode[n][3])][1];
	ccell[n][3][2]=ptcoor[abs(a->pvccnode[n][3])][2];
	
	if(abs(a->pvccnode[n][3])<p->pointnum)
	{
	ccijk[n][3][0]=ptijk[abs(a->pvccnode[n][3])][0];
	ccijk[n][3][1]=ptijk[abs(a->pvccnode[n][3])][1];
	ccijk[n][3][2]=ptijk[abs(a->pvccnode[n][3])][2];
	
	ccid[n][3]=1;
	}
	

    if(a->ccedge[n]>4)
    {
	ccell[n][4][0]=ptcoor[abs(a->pvccnode[n][4])][0];
	ccell[n][4][1]=ptcoor[abs(a->pvccnode[n][4])][1];
	ccell[n][4][2]=ptcoor[abs(a->pvccnode[n][4])][2];
	
	if(abs(a->pvccnode[n][4])<p->pointnum)
	{
	ccijk[n][4][0]=ptijk[abs(a->pvccnode[n][4])][0];
	ccijk[n][4][1]=ptijk[abs(a->pvccnode[n][4])][1];
	ccijk[n][4][2]=ptijk[abs(a->pvccnode[n][4])][2];
	
	ccid[n][4]=1;
	}
	
        if(a->ccedge[n]>5)
        {
        ccell[n][5][0]=ptcoor[abs(a->pvccnode[n][5])][0];
		ccell[n][5][1]=ptcoor[abs(a->pvccnode[n][5])][1];
		ccell[n][5][2]=ptcoor[abs(a->pvccnode[n][5])][2];
		
		if(abs(a->pvccnode[n][5])<p->pointnum)
		{
		ccijk[n][5][0]=ptijk[abs(a->pvccnode[n][5])][0];
		ccijk[n][5][1]=ptijk[abs(a->pvccnode[n][5])][1];
		ccijk[n][5][2]=ptijk[abs(a->pvccnode[n][5])][2];
		
		ccid[n][5]=1;
		}
		
            if(a->ccedge[n]>6)
            {
            ccell[n][6][0]=ptcoor[abs(a->pvccnode[n][6])][0];
			 ccell[n][6][1]=ptcoor[abs(a->pvccnode[n][6])][1];
			 ccell[n][6][2]=ptcoor[abs(a->pvccnode[n][6])][2];
			 
			 if(abs(a->pvccnode[n][6])<p->pointnum)
			 {
			 ccijk[n][6][0]=ptijk[abs(a->pvccnode[n][6])][0];
			 ccijk[n][6][1]=ptijk[abs(a->pvccnode[n][6])][1];
			 ccijk[n][6][2]=ptijk[abs(a->pvccnode[n][6])][2];
			
			 ccid[n][6]=1;
			 }

                if(a->ccedge[n]>7)
                {
                ccell[n][7][0]=ptcoor[abs(a->pvccnode[n][7])][0];
				 ccell[n][7][1]=ptcoor[abs(a->pvccnode[n][7])][1];
				 ccell[n][7][2]=ptcoor[abs(a->pvccnode[n][7])][2];
				 
				 if(abs(a->pvccnode[n][7])<p->pointnum)
				 {
				 ccijk[n][7][0]=ptijk[abs(a->pvccnode[n][7])][0];
				 ccijk[n][7][1]=ptijk[abs(a->pvccnode[n][7])][1];
				 ccijk[n][7][2]=ptijk[abs(a->pvccnode[n][7])][2];
				
				 ccid[n][7]=1;
				 }
                }
            }
        }
    }

	}
	
	p->del_Darray(ptcoor,p->pointnum+p->ccptnum,3);
	p->del_Iarray(ptijk,p->pointnum,3);
}
