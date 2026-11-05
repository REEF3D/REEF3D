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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"momentum.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_hullbody.h"
#include<sys/stat.h>
#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
#include"net_interface.h"
#include"nhflow_reinidisc_fsf.h"

sixdof_obj_nhflow::sixdof_obj_nhflow(lexer *p, ghostcell *pgc, int number) : sixdof_obj_2D(p,pgc,number)
{
    if(p->A10==5)
    {
    pnhfrdisc = new nhflow_reinidisc_fsf(p);
    
    p->Iarray(IO,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CL,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(CR,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(FRK1,p->imax*p->jmax*(p->kmax+2));
    p->Darray(DTT,p->imax*p->jmax*(p->kmax+2));
    p->Darray(LL,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(fsf,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(vert,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(nflag,p->imax*p->jmax*(p->kmax+2));
    }
}

sixdof_obj_nhflow::~sixdof_obj_nhflow()
{
}


void sixdof_obj_nhflow::initialize_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF_df_ini "<<endl;
    
    if(p->mpirank==0)
    mkdir("./REEF3D_NHFLOW_6DOF",0777);
    
    
    // mean horizontal spacing and band width of the ray cast
    ini_band(p,pgc);
    
    
    // Initialise folder structure
    if(p->X50==1)
	print_ini_vtp(p,pgc);
    
    if(p->X50==2)
    print_ini_stl(p,pgc);
 
    // Initialise processor boundaries
    ini_parallel(p,pgc);
    
    // Initialise objects
	objects_create(p,pgc);
    
    // Initialise fbvel
	ini_fbvel(p,pgc);
    
    // Level Set for floating body
    ray_cast(p,d,pgc);
	nhflow_reini_RK2(p,d,pgc,d->FB);
    
    // Calculate geometrical properties
    geometry_parameters_nhflow(p,d,pgc);
    
    // Initialise position of bodies
    iniPosition_RBM(p,pgc);
	
	// Recalculate distances
	ray_cast(p,d,pgc);
	nhflow_reini_RK2(p,d,pgc,d->FB);
    
    // Initialise global variables
	update_fbvel(p,pgc);
    
    // X 15 8: the hull as a sharp thin body (nhflow_hullbody), the thin-body service with it
    if(p->X15==8 && d->hullbody==nullptr)
    d->hullbody = new nhflow_hullbody(p,d,pgc,this);

    // Print initial body (not on a restart from a checkpoint, I 45 1: printed already)
    if(p->I45!=1)
    {
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
    }

	// Mooring
	if(p->X310==0)
	{
		pmooring.push_back(new mooring_void());
	}
	else
	{
		pgc->bcast_int(&p->mooring_count,1);	

		Xme.resize(p->mooring_count);
		Yme.resize(p->mooring_count);
		Zme.resize(p->mooring_count);
		Kme.resize(p->mooring_count);
		Mme.resize(p->mooring_count);
		Nme.resize(p->mooring_count);

		if(p->mpirank==0)
		mkdir("./REEF3D_NHFLOW_6DOF_Mooring",0777);	

		pmooring.reserve(p->mooring_count);
		X311_xen.resize(p->mooring_count,0.0);
		X311_yen.resize(p->mooring_count,0.0);
		X311_zen.resize(p->mooring_count,0.0);
			
		for (int i=0; i < p->mooring_count; i++)
		{
			if(p->X310==1)
			{
				pmooring.push_back(new mooring_Catenary(i));
			}	
			else if(p->X310==2)
			{
				pmooring.push_back(new mooring_barQuasiStatic(i)); 
			}	
			else if(p->X310==3)
			{
                pmooring.push_back(new mooring_dynamic(i));
			}
			else if(p->X310==4)
			{
				pmooring.push_back(new mooring_Spring(i));
			}
		
			// Fairlead offset in the body frame: mooringForces() applies R_*offset,
			// so remove the initial rotation (X101) here to avoid rotating twice
			Eigen::Vector3d fl(p->X311_xe[i] - p->xg, p->X311_ye[i] - p->yg, p->X311_ze[i] - p->zg);
			fl = R_.transpose()*fl;
			X311_xen[i] = fl(0);
			X311_yen[i] = fl(1);
			X311_zen[i] = fl(2);
		
			pmooring[i]->initialize(p,pgc);
		}
	}	

    // Net
    pnetinter->initialize_nhflow(p,d,pgc);
    
    // Membranes (X 330): node positions in the body frame
    if(p->X330>0)
    pnetinter->membrane_attach_nhflow(p,c_,R_);
    
    if(p->X320>0)
    {
    Xne.resize(p->net_count);
    Yne.resize(p->net_count);
    Zne.resize(p->net_count);
    Kne.resize(p->net_count);
    Mne.resize(p->net_count);
    Nne.resize(p->net_count);
    }
    
    // Initial hydrodynamic + gravity load. Xe..Ne are otherwise first set in reforcing(), i.e. after
    // the first pressure solve, so the first RK stage of the first time step would use
    // uninitialised values.
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xd=Yd=Zd=Kd=Md=Nd=0.0;
    
    // porous floating body (X 16)
    if(p->X16==1)
    {
    porosity_nhflow(p,d,pgc);
    
        if(p->mpirank==0)
        {
        cout<<"6DOF porous floating body: n="<<p->X16_n<<" d50="<<p->X16_d50
            <<" alpha="<<p->X16_alpha<<" beta="<<p->X16_beta<<endl;
        cout<<"6DOF porous floating body: X 21 / X 22 define the mass of the solid skeleton,"
            <<" i.e. bulk density rho_s*(1-n), buoyancy acts on (1-n)*V_sub"<<endl;
        }
    }
    
    if(p->X10==1)
    hydrodynamic_forces_nhflow(p,d,pgc,d->WL,false);
    
    // power take-off (X 500): explicit (X 500 1) or with the stiff-element Jacobians on the
    // left-hand side (X 500 2, apply_added_mass with A = 0); two-way bodies only (X 10 1)
    ini_pto(p,pgc);
    
    if(pto_on_ && p->X10!=1 && p->mpirank==0)
    cout<<"6DOF PTO: X 500 needs a two-way floating body (X 10 1) in NHFLOW to act on the motion"<<endl;
}



void sixdof_obj_nhflow::ini_band(lexer *p, ghostcell *pgc)
{
    // DSM: mean horizontal spacing (ray cast band NB*DSM, forcing widths); NB: band in cells.
    // Shared by the floating body (initialize_nhflow) and the wavemaker (initialize_wavemaker).
    int num=0;
    DSM=0.0;
    
    if(p->j_dir==0)
    SLICELOOP4
    {
    DSM += p->DXN[IP];
        
    ++num;
    }
    
    if(p->j_dir==1 && p->geo_curv==0)
    SLICELOOP4
    {
    DSM += 0.5*(p->DXN[IP] + p->DYN[JP]);
        
    ++num;
    }

    // curvilinear grid (G 50): the physical mean widths of the cells
    if(p->j_dir==1 && p->geo_curv==1)
    SLICELOOP4
    {
    DSM += p->geo_width(i,j);

    ++num;
    }

    DSM = pgc->globalsum(DSM);
    num = pgc->globalisum(num);
    
    DSM = DSM/double(MAX(num,1));
    
    NB = MAX(3.0, 2.0*p->A526 + 1.0);
}
