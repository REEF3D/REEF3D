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

#include"fnpf_discrete_weights.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"
#include"vec.h"

fnpf_discrete_weights::fnpf_discrete_weights(lexer* p)
{
}

fnpf_discrete_weights::~fnpf_discrete_weights()
{
}

void fnpf_discrete_weights::ck_weights(lexer *p, double **ck, double *pos, int numpt, int order, int accuracy, int id)
{    
    double c1, c2, c3, c4, c5;
    double z;
    double x[20],c[20][20];
    
    int nd,mn;
    int r,s,t;
    
    nd = accuracy+1; 
    
    // ini
    for(int qn=0;qn<numpt;++qn)
    for(r=0;r<nd;++r)
    ck[qn+marge][r] = 0.0;
        
    for(int qn=0;qn<numpt;++qn)
    {
        nd = accuracy+1; 
        z = pos[qn+marge];
        
        //fill x regular
        for(r=0;r<nd;++r)
        x[r] = pos[qn+marge-accuracy/2 + r];
        
        if(accuracy>2)
        {
            
            // fill beginning
            if((id==1||id==2||id==3) && qn==0)
            {
                nd = accuracy; 
                
                for(r=0;r<nd;++r)
                x[r] = pos[qn+marge - 1 + r];
                
            }
            
            // fill end
            
            if((id==1||id==2||id==3) && qn==numpt-1)
            {
                nd = accuracy;
                
                for(r=0;r<nd;++r)
                x[r] = pos[qn+marge - 2 + r];
                
            }
            
            if(id==3 && qn==numpt-1)
            {
                nd = accuracy;
                
                for(r=0;r<nd;++r)
                x[r] = pos[qn+marge - 2 + r];
                
            }
        }
        
        // fill vertical vel
        if(id==6)
        for(r=0;r<nd;++r)
        x[r] = pos[qn+marge-accuracy + r];
 
            
        //---------------------
        // Algorithm
        c1 = 1.0;
        c4 = x[0] - z;
        for(r=0;r<=nd; ++r)
        for(s=0;s<=nd; ++s)
            
        c[r][s] = 0.0;
        
        c[0][0] = 1.0;

        for(r=1;r<nd; ++r)
        {
        
          mn = MIN(r,order);
          c2 = 1.0;
          c5 = c4;
          c4 = x[r] - z;

          for(s=0; s<r; ++s)
          {
                c3 = x[r] - x[s];
                c2 = c2*c3;
                
                if(s==r-1) 
                {
                for(t=mn;t>=1;t--)
                c[r][t] = c1*(double(t)*c[r-1][t-1] - c5*c[r-1][t])/c2;
         
                c[r][0] = -c1*c5*c[r-1][0]/c2;
                }
                
                for(t=mn;t>=1;t--)
                c[s][t] = (c4*c[s][t] - double(t)*c[s][t-1])/c3;
           
                c[s][0] = c4*c[s][0]/c3;
          }
          c1 = c2;
        }

        //------------------------
        
  
        // write into ct array:
        for(r=0;r<nd;++r)
        ck[qn+marge][r] = -c[r][2];
        
        if(accuracy>2)
        {
            // write beginning stencils
            if((id==1||id==2||id==3) && qn==0)
            {
            ck[qn+marge][0] = 0.0;
            
            for(r=0;r<nd;++r)
            ck[qn+marge][r+1] = -c[r][2];
            }
            
            // write end stencils
            if((id==1||id==2||id==3) && qn==numpt-1)
            {
            
            for(r=0;r<nd;++r)
            ck[qn+marge][r] = -c[r][2];
            }
        }
        
        if(id==6)
        for(r=0;r<nd;++r)
        ck[qn+marge][r] = c[r][1];
        
        
        
        // Debug Print
        /*
        if(id==3)
        if(p->mpirank==0)
        {
        for(r=0;r<accuracy+1;++r)
        cout<<ck[qn+marge][r]<<" ";
        cout<<endl;
        }*/
        /*
        if(id==6)
        if(p->mpirank==0)
        {
        for(r=0;r<accuracy+1;++r)
        cout<<ck[qn+marge][r]<<" ";
        cout<<endl;
        
        //cout<<"z: "<<z<<endl;
        }*/
        
        
        
    }


}




