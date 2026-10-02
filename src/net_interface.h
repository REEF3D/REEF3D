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

#ifndef SIXDOF_CFD_H_
#define SIXDOF_CFD_H_
#include <Eigen/Dense>
#include"6DOF.h"
#include"6DOF_obj.h"
#include<vector>

class mooring;
class net;
class net_membrane;
class slice;
class ddweno_f_nug;

using namespace std;

class net_interface : public increment
{
public:
	net_interface(lexer*, ghostcell*);
	virtual ~net_interface();

    void start_cfd(lexer*, fdm*, ghostcell*, double, Eigen::Matrix3d, bool) {};
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*, double, Eigen::Matrix3d, bool) {};

    void netForces_cfd(lexer*, fdm*, ghostcell*, double, Eigen::Matrix3d, vector<double>&, vector<double>&, 
                        vector<double>&, vector<double>&, vector<double>&, vector<double>&, bool);
                        
    void netForces_nhflow(lexer*, fdm_nhf*, ghostcell*, double, Eigen::Matrix3d, vector<double>&, vector<double>&, 
                        vector<double>&, vector<double>&, vector<double>&, vector<double>&, bool);
    

    void dlm_cfd(lexer*, fdm*, ghostcell*, int);
    void dlm_nhflow(lexer*, fdm_nhf*, ghostcell*, int);
    
    void initialize_cfd(lexer*, fdm*, ghostcell*);
    void initialize_nhflow(lexer*, fdm_nhf*, ghostcell*);
    
    // impermeable membranes (X 330, membrane.dat), NHFLOW
    void membrane_ini_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void membrane_forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void membrane_reaction_nhflow(lexer*, fdm_nhf*, ghostcell*, double, slice&, bool);
    void membrane_pgrad(lexer*, fdm_nhf*, double, double*, double*, double*, slice&, int);
    bool membrane_iterated();
    void membrane_reforce_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    bool membrane_couple_nhflow(lexer*, fdm_nhf*, ghostcell*, int, double, slice&, int);
    void membrane_attach_nhflow(lexer*, const Eigen::Vector3d&, const Eigen::Matrix3d&);
    void membrane_body_nhflow(const Eigen::Vector3d&, const Eigen::Matrix3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
    void membraneForces_nhflow(lexer*, const Eigen::Vector3d&, const Eigen::Matrix3d&, double&, double&, double&, double&, double&, double&);
    double membrane_addedmass_nhflow(lexer*);
    Eigen::Matrix3d membrane_stiffness_nhflow(lexer*);
    
    typedef vector<Eigen::Vector3d> EigenMat;
    
    
private:
    vector<net*> pnet;
    vector<net_membrane*> pmem;
   
    field1 kernel_x;
    field2 kernel_y;
    field3 kernel_z;
    
    double *KX,*KY,*KZ;
   
   double kernel_peskin(const double&);
   
};

#endif
