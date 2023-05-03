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

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_breaking::breaking_algorithm(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    int ii,jj;
    
    if(p->A350>=0)
    if(p->count>count_n)
    {
    SLICELOOP4
    c->breaking(i,j)=0;
    
    count_n=p->count;
    }
    
    if(p->A350>=0)
    SLICELOOP4
    {
    bx(i,j)=0;
    by(i,j)=0;
    }
    
    pgc->gcsl_start4int(p,c->breaking,50);
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    
    
    if((p->A351==2 || p->A351==3) && p->count>1)
    {
    SLICELOOP4
    {
            // x
            if(c->Ex(i,j)   < -p->A355)
            {
                ii=i;
                
                bx(i,j) = 10;
                bx(i+1,j) = 10;
                bx(i-1,j) = 10;
                bx(i-2,j) = 10;

                i=ii;
            }
            
            if(c->Ex(i,j)   > p->A355)
            {
                ii=i;
                
                bx(i,j) = 20;
                bx(i-1,j) = 20;
                bx(i+1,j) = 20;
                bx(i+2,j) = 20;

                i=ii;
            }
            
            // y
            if(p->j_dir==1)
            if( c->Ey(i,j)   < -p->A355)
            {
                jj=j;
                
                by(i,j) = 10;
                by(i,j+1) = 10;
                by(i,j-1) = 10;
                by(i,j-2) = 10;

                j=jj;
            }
            
            if(p->j_dir==1)
            if( c->Ey(i,j)   > p->A355)
            {
                jj=j;
                
                by(i,j) = 20;
                by(i,j-1) = 20;
                by(i,j+1) = 20;
                by(i,j+2) = 20;

                j=jj;
            }
            
    }
    
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    
    
        SLICELOOP4
        if(bx(i,j)>0 || by(i,j)>0)
        {
        c->breaking(i,j)=1;
        }
    }
    
    
    
    
    if((p->A351==1 || p->A351==3) && p->count>1)
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
            
            if((eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A354*sqrt(9.81*c->WL(i,j)))
            {

                c->breaking(i-1,j)=2;
                c->breaking(i-2,j)=2;
                c->breaking(i,j)=2;
                c->breaking(i+1,j)=2;
                c->breaking(i+2,j)=2;
                
                if(p->j_dir==1)
                {
                c->breaking(i,j-2)=2;
                c->breaking(i,j-1)=2;
                c->breaking(i,j+1)=2;
                c->breaking(i,j+2)=2;
                }
            }
    }
    
    
    
    // -------------------
    if(p->A350==1)
    {
        SLICELOOP4
        c->vb(i,j) = 0.0;
        
        // coastline
        SLICELOOP4
        {
            
            if(c->coastline(i,j)>=0.0 && p->A346>0.0)
            {
                db = c->coastline(i,j);
                
                if(db<dist3)
                {
                c->vb(i,j) = rb3(p,db)*p->A346;
            
                }
            }
        }
        
        if(p->j_dir==0)
        SLICELOOP4
        {   
            
            if(c->breaking(i,j)>=1 || c->breaking(i-1,j)>=1 || c->breaking(i+1,j)>=1)
            c->vb(i,j) = p->A365*double(c->breaking(i,j));
            
            
            if(c->breaking(i,j)==0 &&(c->breaking(i-2,j)>=1 || c->breaking(i+2,j)>=1))
            c->vb(i,j) = 0.5*p->A365;
        }

        if(p->j_dir==1)
        SLICELOOP4
        {   
            
            if(c->breaking(i,j)>=1 || c->breaking(i-1,j)>=1 || c->breaking(i+1,j)>=1 || c->breaking(i,j-1)>=1 || c->breaking(i,j+1)>=1)
            c->vb(i,j) = p->A365*double(c->breaking(i,j));
            
            if(c->breaking(i,j)==0 &&( c->breaking(i-1,j-1)>=1 || c->breaking(i-1,j+1)>=1 || c->breaking(i+1,j-1)>=1 || c->breaking(i+1,j+1)>=1
           || c->breaking(i-2,j)>=1 || c->breaking(i+2,j)>=1 || c->breaking(i,j-2)>=1 || c->breaking(i,j+2)>=1))
            c->vb(i,j) = 0.5*p->A365;
        }
        
        if(p->A352==1)
        SLICELOOP4
        if(c->breaking(i,j)==2)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
        
        if(p->A352==2)
        SLICELOOP4
        if(c->breaking(i,j)==1)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
        
        if(p->A352==3)
        SLICELOOP4
        if(c->breaking(i,j)>=1)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
        
    pgc->gcsl_start4(p,c->vb,1);
    }
    
    if(p->A350==2)
    SLICELOOP4
    {
        if(c->breaking(i,j)>=1 || c->breaking(i-1,j)>=1 || c->breaking(i+1,j)>=1 || c->breaking(i,j-1)>=1 || c->breaking(i,j+1)>=1)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
    }
    
    SLICELOOP4
    c->breaklog(i,j)=0;
    
    // breaklog
    int count=0;
    
    SLICELOOP4
    if(c->breaking(i,j)>0)
    {
    c->breaklog(i,j)=1;
    ++count;
    }
    
    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
}


double fnpf_breaking::rb3(lexer *p, double x)
{
    double r=0.0;

    x=(dist3-fabs(x))/(dist3);
    x=MAX(x,0.0);
    
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);

	return r;
}

double fnpf_breaking::rb4(lexer *p, double x)
{
    double r=0.0;

    x=(dist4-fabs(x))/(dist4);
    x=MAX(x,0.0);
    
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);

	return r;
}
