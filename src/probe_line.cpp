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

#include"probe_line.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<sys/stat.h>
#include<sys/types.h>

probe_line::probe_line(lexer *p, fdm* a, ghostcell *pgc) : probenum(p->P62), eps(1.0e-10*p->DXM)
{
	linecount=0;
	
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_ProbeLine",0777);
	
	if(p->mpirank==0 && p->P62>0)
	cout<<"probeline_num: "<<probenum<<endl;

	lineout = new ofstream[p->P62];
	
	p->Darray(length,p->P62);
	p->Darray(ds,p->P62);
	p->Darray(norm,p->P62);
	p->Iarray(elnum,p->P62);
	p->Iarray(totelnum,p->P62);
	p->Iarray(totelnum_all,p->P62);
	p->Iarray(elnum_all,p->P62,p->M10);
	
	for(n=0;n<p->P62;++n)
	length[n] = sqrt(pow(p->P62_xe[n]-p->P62_xs[n],2.0) + pow(p->P62_ye[n]-p->P62_ys[n],2.0) + pow(p->P62_ze[n]-p->P62_zs[n],2.0));
	
	for(n=0;n<p->P62;++n)
	totelnum[n] = conv(length[n]/p->DXM) + 1;
    
    // -------
	
	for(n=0;n<p->P62;++n)
	ds[n] = length[n]/double(totelnum[n]-1);

	maxelnum=0;
	for(n=0;n<p->P62;++n)
	maxelnum = MAX(maxelnum,totelnum[n]);
    
    maxelnum += p->M10+6;
 
    // -------
	
	p->Iarray(active,p->P62,maxelnum);
    p->Iarray(active_all,p->P62,maxelnum);
    p->Darray(xs,p->P62,maxelnum);
    p->Darray(xs_all,p->P62,maxelnum);
    p->Iarray(iloc,p->P62,maxelnum);
    p->Iarray(iloc_all,p->P62,maxelnum);
    p->Iarray(kloc,p->P62,maxelnum);
    p->Iarray(kloc_all,p->P62,maxelnum);
	p->Iarray(elid,p->P62,maxelnum);
	p->Iarray(elid_all,p->P62,maxelnum);
	
	p->Iarray(displ,p->P62,p->M10+1);
	
	domain_xs = p->originx;
	domain_ys = p->originy;
	domain_zs = p->originz;
	
	domain_xe = p->endx;
	domain_ye = p->endy;
	domain_ze = p->endz;
	

	// xdir
	if(p->nb1==-2)
	eps_xs = -1.0e-10*p->DXM;
	
	if(p->nb1>=0)
	eps_xs = 1.0e-10*p->DXM;
	
	if(p->nb4==-2)
	eps_xe = 1.0e-10*p->DXM;
	
	if(p->nb4>=0)
	eps_xe = 1.0e-10*p->DXM;
	
	// ydir
	if(p->nb3==-2)
	eps_ys = -1.0e-10*p->DXM;
	
	if(p->nb3>=0)
	eps_ys = 1.0e-10*p->DXM;
	
	if(p->nb2==-2)
	eps_ye = 1.0e-10*p->DXM;
	
	if(p->nb2>=0)
	eps_ye = 1.0e-10*p->DXM;
	
	// zdir
	if(p->nb5==-2)
	eps_zs = -1.0e-10*p->DXM;
	
	if(p->nb5>=0)
	eps_zs = 1.0e-10*p->DXM;
	
	if(p->nb6==-2)
	eps_ze = 1.0e-10*p->DXM;
	
	if(p->nb6>=0)
	eps_ze = 1.0e-10*p->DXM;
	

	//------------------------------
	ini_global_location(p,a,pgc);
	//------------------------------
	
	maxlocalelnum=0;
	
	for(n=0;n<p->P62;++n)
	maxlocalelnum = MAX(maxlocalelnum,elnum[n]);
	
    p->Darray(U,p->P62,maxlocalelnum);
	p->Darray(V,p->P62,maxlocalelnum);
	p->Darray(W,p->P62,maxlocalelnum);
	p->Darray(P,p->P62,maxlocalelnum);
	p->Darray(K,p->P62,maxlocalelnum);
	p->Darray(E,p->P62,maxlocalelnum);
	p->Darray(VT,p->P62,maxlocalelnum);
	p->Darray(LS,p->P62,maxlocalelnum);
	
    p->Darray(U_all,p->P62,maxelnum);
	p->Darray(V_all,p->P62,maxelnum);
	p->Darray(W_all,p->P62,maxelnum);
	p->Darray(P_all,p->P62,maxelnum);
	p->Darray(K_all,p->P62,maxelnum);
	p->Darray(E_all,p->P62,maxelnum);
	p->Darray(VT_all,p->P62,maxelnum);
	p->Darray(LS_all,p->P62,maxelnum);
	p->Darray(VAL,p->P62,maxelnum);
	
	p->Iarray(flag,p->P62,maxlocalelnum);
    p->Iarray(flag_all,p->P62,maxelnum);
	
	//------------------------------
	ini_location(p,a,pgc);
	//------------------------------
}

probe_line::~probe_line()
{
	for(n=0;n<probenum;++n)
    lineout[n].close();
}

void probe_line::start(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb)
{
	char name[250];
	int num;
	
	if(p->P15==1)
    num = linecount;

    if(p->P15==2)
    num = p->count;

	for(n=0;n<p->P62;++n)
	{
		
		if(p->mpirank==0)
		{
			// open file
			if(p->P14==0)
			{
			if(num<10)
			sprintf(name,"REEF3D-CFD-probeline-%i-00000%i.dat",n+1,num);

			if(num<100&&num>9)
			sprintf(name,"REEF3D-CFD-probeline-%i-0000%i.dat",n+1,num);

			if(num<1000&&num>99)
			sprintf(name,"REEF3D-CFD-probeline-%i-000%i.dat",n+1,num);

			if(num<10000&&num>999)
			sprintf(name,"REEF3D-CFD-probeline-%i-00%i.dat",n+1,num);

			if(num<100000&&num>9999)
			sprintf(name,"REEF3D-CFD-probeline-%i-0%i.dat",n+1,num);

			if(num>99999)
			sprintf(name,"REEF3D-CFD-probeline-%i-%i.dat",n+1,num);
			}
			
			if(p->P14==1)
			{
			if(num<10)
			sprintf(name,"./REEF3D_CFD_ProbeLine/REEF3D-CFD-probeline-%i-00000%i.dat",n+1,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_CFD_ProbeLine/REEF3D-CFD-probeline-%i-0000%i.dat",n+1,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_CFD_ProbeLine/REEF3D-CFD-probeline-%i-000%i.dat",n+1,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_CFD_ProbeLine/REEF3D-CFD-probeline-%i-00%i.dat",n+1,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_CFD_ProbeLine/REEF3D-CFD-probeline-%i-0%i.dat",n+1,num);

			if(num>99999)
			sprintf(name,"./REEF3D_CFD_ProbeLine/REEF3D-CFD-probeline-%i-%i.dat",n+1,num);
			}
			
			lineout[n].open(name);

			lineout[n]<<"Line Probe ID:  "<<n+1<<endl<<endl;
			lineout[n]<<"x_start     x_end     y_start     y_end     z_start     z_end"<<endl;
			
			lineout[n]<<p->P62_xs[n]<<"\t "<<p->P62_xe[n]<<"\t "<<p->P62_ys[n]<<"\t "<<p->P62_ye[n]<<"\t "<<p->P62_zs[n]<<"\t "<<p->P62_ze[n]<<endl;

			lineout[n]<<endl<<endl;
			
			lineout[n]<<"X \t Y \t Z \t U \t V \t W \t P \t Kin \t Epsilon/Omega \t Eddy-viscosity \t LS"<<endl;
			
			lineout[n]<<endl<<endl;
		}
		
		for(q=0;q<elnum[n];++q)
		{
			if(flag[n][q]>0)
			{
			xp = p->P62_xs[n] + double(elid[n][q])*ds[n]*(p->P62_xe[n]-p->P62_xs[n])/(norm[n]>eps?norm[n]:1.0e20);
			yp = p->P62_ys[n] + double(elid[n][q])*ds[n]*(p->P62_ye[n]-p->P62_ys[n])/(norm[n]>eps?norm[n]:1.0e20);
			zp = p->P62_zs[n] + double(elid[n][q])*ds[n]*(p->P62_ze[n]-p->P62_zs[n])/(norm[n]>eps?norm[n]:1.0e20);
			
			U[n][q] = p->ccipol1(a->u, xp, yp, zp);
			V[n][q] = p->ccipol2(a->v, xp, yp, zp);
			W[n][q] = p->ccipol3(a->w, xp, yp, zp);
			P[n][q] = p->ccipol4_a(a->press, xp, yp, zp) - p->pressgage;
			K[n][q] = pturb->ccipol_kinval(p, pgc, xp, yp, zp);
			E[n][q] = pturb->ccipol_epsval(p, pgc, xp, yp, zp);
			VT[n][q] = p->ccipol4_a(a->eddyv, xp, yp, zp);
			LS[n][q] = p->ccipol4(a->phi, xp, yp, zp);
			}
			
			if(flag[n][q]==0)
			{
			U[n][q] = 0.0;
			V[n][q] = 0.0;
			W[n][q] = 0.0;
			P[n][q] = 0.0;
			K[n][q] = 0.0;
			E[n][q] = 0.0;
			VT[n][q] = 0.0;
			LS[n][q] = 0.0;
			}
		}
		
		pgc->gatherv_double(U[n],elnum[n],VAL[n],elnum_all[n],displ[n]); 
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		U_all[n][elid_all[n][q]] = VAL[n][q];
		

		pgc->gatherv_double(V[n],elnum[n],VAL[n],elnum_all[n],displ[n]);
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		V_all[n][elid_all[n][q]] = VAL[n][q];
		
		
		pgc->gatherv_double(W[n],elnum[n],VAL[n],elnum_all[n],displ[n]);
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		W_all[n][elid_all[n][q]] = VAL[n][q];
		
		
		pgc->gatherv_double(P[n],elnum[n],VAL[n],elnum_all[n],displ[n]);
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		P_all[n][elid_all[n][q]] = VAL[n][q];
		
		
		pgc->gatherv_double(K[n],elnum[n],VAL[n],elnum_all[n],displ[n]);
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		K_all[n][elid_all[n][q]] = VAL[n][q];
		
		
		pgc->gatherv_double(E[n],elnum[n],VAL[n],elnum_all[n],displ[n]);
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		E_all[n][elid_all[n][q]] = VAL[n][q];
		
		
		pgc->gatherv_double(VT[n],elnum[n],VAL[n],elnum_all[n],displ[n]);
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		VT_all[n][elid_all[n][q]] = VAL[n][q];
		
		
		pgc->gatherv_double(LS[n],elnum[n],VAL[n],elnum_all[n],displ[n]);
		
		if(p->mpirank==0)
		for(q=0;q<totelnum_all[n];++q)
		LS_all[n][elid_all[n][q]] = VAL[n][q];
        
        
        pgc->gatherv_int(active[n],elnum[n],active_all[n],elnum_all[n],displ[n]);
        pgc->gatherv_int(flag[n],elnum[n],flag_all[n],elnum_all[n],displ[n]);
        pgc->gatherv_int(iloc[n],elnum[n],iloc_all[n],elnum_all[n],displ[n]);
        pgc->gatherv_int(kloc[n],elnum[n],kloc_all[n],elnum_all[n],displ[n]);
        pgc->gatherv_double(xs[n],elnum[n],xs_all[n],elnum_all[n],displ[n]);
		


		if(p->mpirank==0)
		{
			for(q=0;q<totelnum_all[n];++q)
			{
				xp = p->P62_xs[n] + double(elid_all[n][q])*ds[n]*(p->P62_xe[n]-p->P62_xs[n])/(norm[n]>eps?norm[n]:1.0e20);
				yp = p->P62_ys[n] + double(elid_all[n][q])*ds[n]*(p->P62_ye[n]-p->P62_ys[n])/(norm[n]>eps?norm[n]:1.0e20);
				zp = p->P62_zs[n] + double(elid_all[n][q])*ds[n]*(p->P62_ze[n]-p->P62_zs[n])/(norm[n]>eps?norm[n]:1.0e20);
				
			lineout[n]<<xp<<" \t "<<yp<<" \t "<<zp<<" \t "<<U_all[n][q]<<" \t "<<V_all[n][q]<<" \t "<<W_all[n][q];
			lineout[n]<<" \t "<<P_all[n][q]<<" \t "<<K_all[n][q]<<" \t "<<E_all[n][q]<<" \t "<<VT_all[n][q]<<" \t "<<LS_all[n][q]<<endl;//" \t \t "<<active_all[n][q]<<" \t "<<flag_all[n][q]<<" \t "<<xs_all[n][q]<<" \t "<<iloc_all[n][q]<<" \t "<<kloc_all[n][q]<<endl;
			}
		lineout[n].close();
		}
	}
	
	++linecount;
}

void probe_line::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void probe_line::ini_global_location(lexer *p, fdm *a, ghostcell *pgc)
{
	
	for(n=0;n<p->P62;++n)
	for(q=0;q<totelnum[n];++q)
	active[n][q]=0;
	
	for(n=0;n<p->P62;++n)
	norm[n] = sqrt(pow(p->P62_xe[n]-p->P62_xs[n],2.0) + pow(p->P62_ye[n]-p->P62_ys[n],2.0) + pow(p->P62_ze[n]-p->P62_zs[n],2.0));
	
	for(n=0;n<p->P62;++n)
	{
		t=0.0;
		count=0;
		for(q=0;q<totelnum[n];++q)
		{
		xloc = p->P62_xs[n] + t*ds[n]*(p->P62_xe[n]-p->P62_xs[n])/(norm[n]>eps?norm[n]:1.0e20);
		yloc = p->P62_ys[n] + t*ds[n]*(p->P62_ye[n]-p->P62_ys[n])/(norm[n]>eps?norm[n]:1.0e20);
		zloc = p->P62_zs[n] + t*ds[n]*(p->P62_ze[n]-p->P62_zs[n])/(norm[n]>eps?norm[n]:1.0e20);
            
            if(p->j_dir==0)
			if(xloc>domain_xs+eps_xs && xloc<=domain_xe+eps_xe 
			&& zloc>domain_zs+eps_zs && zloc<=domain_ze+eps_ze)
			{
			active[n][q]=1;
			++count;
			}
            
            if(p->j_dir==1)
			if(xloc>domain_xs+eps_xs && xloc<=domain_xe+eps_xe 
			&& yloc>domain_ys+eps_ys && yloc<=domain_ye+eps_ye 
			&& zloc>domain_zs+eps_zs && zloc<=domain_ze+eps_ze)
			{
			active[n][q]=1;
			++count;
			}
		
		t+=1.0;	
		}
		elnum[n]=count;
	}
    
    // --------------------

	for(n=0;n<p->P62;++n)
	{	
		count=0;
		allcount=0;
		for(q=0;q<totelnum[n];++q)
		{
			if(active[n][q]==1)
			{
			elid[n][count] = allcount;
			++count;
			}
		++allcount;
		}
	}
	
	for(n=0;n<p->P62;++n)
	pgc->gather_int(&elnum[n],1,elnum_all[n],1);
	
	for(n=0;n<p->P62;++n)
	for(q=0;q<p->M10;++q)
	displ[n][q]=0;
	
	for(n=0;n<p->P62;++n)
	for(q=1;q<p->M10;++q)
	displ[n][q]=displ[n][q-1]+elnum_all[n][q-1];
	
	for(n=0;n<p->P62;++n)
	pgc->gatherv_int(elid[n],elnum[n],elid_all[n],elnum_all[n],displ[n]);
	
	if(p->mpirank==0)
	{
	for(n=0;n<p->P62;++n)
	totelnum_all[n]=0;
	
	for(n=0;n<p->P62;++n)
	for(q=0;q<p->M10;++q)
	totelnum_all[n] +=  elnum_all[n][q];
	}
}

void probe_line::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
	int iiloc,jjloc,kkloc;
	
	for(n=0;n<p->P62;++n)
	for(q=0;q<elnum[n];++q)
	flag[n][q]=0;
	
	
	for(n=0;n<p->P62;++n)
	{	
		count=0;
		t=0.0;
		for(q=0;q<totelnum[n];++q)
		{
			if(active[n][q]==1)
			{
            xs[n][q] = p->P62_xs[n] + t*ds[n]*(p->P62_xe[n]-p->P62_xs[n])/(norm[n]>eps?norm[n]:1.0e20);
			iiloc = p->posc_i(p->P62_xs[n] + t*ds[n]*(p->P62_xe[n]-p->P62_xs[n])/(norm[n]>eps?norm[n]:1.0e20));
            
            if(p->j_dir==0)
            jjloc=0;
            
            if(p->j_dir==1)
			jjloc = p->posc_j(p->P62_ys[n] + t*ds[n]*(p->P62_ye[n]-p->P62_ys[n])/(norm[n]>eps?norm[n]:1.0e20));
            
			kkloc = p->posc_k(p->P62_zs[n] + t*ds[n]*(p->P62_ze[n]-p->P62_zs[n])/(norm[n]>eps?norm[n]:1.0e20));
			
            if(p->j_dir==0)
			check=boundcheck_ik(p,a,iiloc,jjloc,kkloc,1);
            
            if(p->j_dir==1)
			check=boundcheck(p,a,iiloc,jjloc,kkloc,1);
            
            //if(n==3 && p->mpirank==3)
            //cout<<p->mpirank<<" n: "<<n<<" i: "<<iiloc<<" j: "<<jjloc<<" k: "<<kkloc<<" check: "<<check<<" p->originx: "<<p->originx<<" p->endx: "<<p->endx<<endl;
			
            iloc[n][q] = iiloc;
            kloc[n][q] = p->mpirank;
            
            if(check==1)
			flag[n][count]=1;
			
			++count;
			}
			
		t+=1.0;
		}
	}
    
    

    //cout<<p->mpirank<<"  POSC_I: "<<p->posc_i(22.0)<<" XN_0: "<<p->XN[0+marge]<<"  "<<p->XN[p->knox+marge]<<" p->originx: "<<p->originx<<" p->endx: "<<p->endx<<endl;
}

int probe_line::conv(double a)
{

int b,c;
double d,diff;

c= int( a);
d=double(c);
diff=a-d;

b=c;

if(diff>0.5)
b=c+1;

if(diff<=-0.5)
b=c-1;

return b;

}
