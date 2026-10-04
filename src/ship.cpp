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

#include"ship.h"
#include"ship_hull.h"
#include"ship_models.h"
#include"6DOF_rigidbody.h"
#include"6DOF_geometry.h"
#include"6DOF_output_dir.h"
#include"lexer.h"
#include<iostream>
#include<sstream>
#include<cstdio>

ship::ship(lexer *p, int number) : id(number), initialized(false),
                                   lpp(0.0), S(0.0), k(0.0), nu(p->W2), B44(0.0), B44q(0.0), Cd(0.0),
                                   thrust(0.0), xthrust(0.0), zthrust(0.0),
                                   friction(1), nstrip(40), lpp_in(false), S_in(false),
                                   xa(0.0), xf(0.0), zw(0.0),
                                   ub(0.0), vb(0.0), wb(0.0), pb(0.0), qb(0.0), rb_(0.0),
                                   Re(0.0), CF(0.0), XF(0.0), Ycf(0.0), Ncf(0.0), Kroll(0.0)
{
    // CFD resolves the wall shear (viscous momentum equations): no correlation-line friction by default
    if(p->A10==6)
    friction = 0;
    
    read(p);
    
    // NHFLOW X 38 1: the hull loads already contain a local skin friction
    if(friction==1 && p->X38==1)
    {
        friction = 0;
        
        if(p->mpirank==0)
        cout<<"ship: X 38 1 (local skin friction in the hull loads), ITTC-1957 friction switched off"<<endl;
    }
}

ship::~ship()
{
    if(out.is_open())
    out.close();
}

void ship::read(lexer *p)
{
    ifstream f("ship.dat");
    
    if(!f.is_open())
    {
        if(p->mpirank==0)
        cout<<"ship: X 350 1 but no ship.dat, defaults are used"<<endl;
        
        return;
    }
    
    string line;
    
    while(getline(f,line))
    {
        const size_t c = line.find('#');
        
        if(c!=string::npos)
        line = line.substr(0,c);
        
        istringstream ls(line);
        string key;
        
        if(!(ls>>key))
        continue;
        
        if(key=="lpp")
        {
            ls>>lpp;
            lpp_in = true;
        }
        else if(key=="wetted_surface")
        {
            ls>>S;
            S_in = true;
        }
        else if(key=="form_factor")
        ls>>k;
        
        else if(key=="friction")
        ls>>friction;
        
        else if(key=="viscosity")
        ls>>nu;
        
        else if(key=="roll_damping")
        ls>>B44>>B44q;
        
        else if(key=="crossflow")
        ls>>Cd;
        
        else if(key=="strips")
        ls>>nstrip;
        
        else if(key=="thrust")
        {
            ls>>thrust;
            
            if(!(ls>>xthrust>>zthrust))
            xthrust = zthrust = 0.0;
        }
        
        else if(p->mpirank==0)
        cout<<"ship: unknown keyword in ship.dat: "<<key<<endl;
    }
}

void ship::ini(lexer *p, const sixdof_rigidbody &b, const sixdof_geometry &g)
{
    // still water level in the body frame of the initial position (level hull)
    zw = p->F60 - b.c(2);
    
    ship_hull::waterline_extent(g.tri_x0,g.tri_y0,g.tri_z0,g.tricount,zw,xa,xf);
    
    if(!lpp_in)
    lpp = xf - xa;
    
    if(!S_in)
    S = ship_hull::wetted_surface(g.tri_x0,g.tri_y0,g.tri_z0,g.tricount,zw,p->j_dir==0);
    
    if(Cd>0.0)
    ship_hull::draft_strips(g.tri_x0,g.tri_y0,g.tri_z0,g.tricount,zw,xa,xf,nstrip,xs,dx,T);
    
    if(p->mpirank==0)
    {
        cout<<"ship "<<id<<": L = "<<lpp<<" m, S = "<<S<<" m^2, waterline x = "<<xa<<" .. "<<xf<<" m (from the CoG)"<<endl;
        cout<<"ship "<<id<<": friction "<<friction<<" (1+k) = "<<1.0+k<<", nu = "<<nu
            <<", roll damping "<<B44<<" "<<B44q<<", cross-flow Cd = "<<Cd<<", thrust "<<thrust<<endl;
    }
    
    initialized = true;
}

void ship::add_load(lexer *p, const sixdof_rigidbody &b, const sixdof_geometry &g, double *F)
{
    if(!initialized)
    ini(p,b,g);
    
    // velocity of the CoG and angular velocity in the ship frame
    Eigen::Matrix<double,6,1> u6;
    b.velocity(u6);
    
    const Eigen::Vector3d uI(u6(0),u6(1),u6(2));
    const Eigen::Vector3d wI(u6(3),u6(4),u6(5));
    const Eigen::Vector3d u = b.R.transpose()*uI;
    const Eigen::Vector3d w = b.R.transpose()*wI;
    
    ub = u(0); vb = u(1); wb = u(2);
    pb = w(0); qb = w(1); rb_ = w(2);
    
    // loads in the ship frame
    Eigen::Vector3d Fs(0.0,0.0,0.0), Ms(0.0,0.0,0.0);
    
    XF = 0.0;
    if(friction==1 && S>0.0 && lpp>0.0)
    XF = ship_models::friction(p->W1,nu,S,lpp,k,ub,Re,CF);
    
    Ycf = Ncf = 0.0;
    if(Cd>0.0)
    ship_models::crossflow(p->W1,Cd,xs,dx,T,vb,rb_,Ycf,Ncf);
    
    Kroll = ship_models::roll_damping(B44,B44q,pb);
    
    Fs(0) = XF + thrust;
    Fs(1) = Ycf;
    
    Ms(0) = Kroll;
    Ms(1) = zthrust*thrust;
    Ms(2) = Ncf;
    
    // inertial frame
    const Eigen::Vector3d FI = b.R*Fs;
    const Eigen::Vector3d MI = b.R*Ms;
    
    for(int n=0; n<3; ++n)
    {
        F[n]   += FI(n);
        F[n+3] += MI(n);
    }
}

void ship::print(lexer *p)
{
    if(p->mpirank!=0 || p->count%p->X19!=0)
    return;
    
    if(!out.is_open())
    {
        char name[1000];
        snprintf(name,sizeof(name),"%s/REEF3D_ship_%i.dat",sixdof_output_dir(p),id);
        out.open(name);
        out<<"time \t u [m/s] \t v [m/s] \t r [rad/s] \t Re \t C_F \t X_F [N] \t Y_cf [N] \t N_cf [Nm] \t K_roll [Nm] \t T [N]"<<endl;
    }
    
    out<<p->simtime<<" \t "<<ub<<" \t "<<vb<<" \t "<<rb_<<" \t "<<Re<<" \t "<<CF<<" \t "<<XF<<" \t "<<Ycf<<" \t "<<Ncf<<" \t "<<Kroll<<" \t "<<thrust<<endl;
}
