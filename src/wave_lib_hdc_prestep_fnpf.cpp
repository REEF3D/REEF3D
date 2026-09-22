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

#include"wave_lib_hdc.h"
#include"lexer.h"
#include"ghostcell.h"

void wave_lib_hdc::wave_prestep_fnpf(lexer *p, ghostcell *pgc)
{
    // only at startup
    if(startup==0)
    {
        deltaT = simtime[1]-simtime[0];
        
        deltaT = deltaT>0.0?deltaT:1.0e20;
        
        t1 = (simtime[1]-(p->wavetime+p->I241))/deltaT;
        t2 = ((p->wavetime+p->I241)-simtime[0])/deltaT;
        
        q1 = diter;
        q2 = diter+1;
        
        if(file_conti==2)
        {
        filename_continuous(p,pgc);
        result.open(name, ios::binary);
        }
    
        read_result_fnpf(p,pgc,E1,F1,q1);
        read_result_fnpf(p,pgc,E2,F2,q2);
        startup=1;
        
        
        // find q1
        while(simtime[q1+1-diter]<=p->wavetime+p->I241)
        {
        ++q1;
        
        if(file_conti==2)
        read_result_fnpf(p,pgc,E1,F1,q1);
        }
            
        // find q2
        while(simtime[q2-diter]<p->wavetime+p->I241)
        {
        ++q2;
        
        if(file_conti==2 )
        read_result_fnpf(p,pgc,E2,F2,q2);
        }
    }
    
    q1n=q1;
    q2n=q2;
    
    // check: open next timestep
           
    // find q1
    while(simtime[q1+1-diter]<=p->wavetime+p->I241)
    ++q1;
        
    // find q2
    while(simtime[q2-diter]<p->wavetime+p->I241)
    ++q2;
    
    if(q2>=numiter+diter)
    endseries=1;
        
    q1=MIN(q1,numiter+diter);
    q2=MIN(q2,numiter+diter);
    
    if(q1==q2)
    ++q2;
    
    
    // single file read 
    if(file_conti==1)
    {
        if(q1!=q1n)
        {
        // Open File 1
        filename_single(p,pgc,q1);
        read_result_fnpf(p,pgc,E1,F1,q1);
        }
        
        if(q2!=q2n)
        {
        // Open File 2
        filename_single(p,pgc,q2);
        read_result_fnpf(p,pgc,E2,F2,q2);
        }
    }
        
    // contiuous file read
    if(file_conti==2)
    {
        if(q1!=q1n)
        fill_conti_fnpf(p,pgc);
        
        if(q2!=q2n)
        read_result_fnpf(p,pgc,E2,F2,q2);
    }
        

        deltaT = simtime[q2-diter]-simtime[q1-diter];
        
        if(p->mpirank==0)
        cout<<"HDC  q1: "<<q1<<" q2: "<<q2<<" t1: "<<t1<<" t2: "<<t2<<" deltaT: "<<deltaT<<" simtime[q1]: "<<simtime[q1-diter]<<" simtime[q2]: "<<simtime[q2-diter]<<endl;

        t1 = (simtime[q2-diter]-(p->wavetime+p->I241))/deltaT;
        t2 = ((p->wavetime+p->I241)-simtime[q1-diter])/deltaT;
        
        
    // time interpolation
    //if(q1!=q1n || q2!=q2n)
    time_interpol_fnpf(p);
    
}
