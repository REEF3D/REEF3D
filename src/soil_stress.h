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
Architect: Hans Bihs
--------------------------------------------------------------------*/
// REEF3D::CFD soil stress model (S 200): wave-induced response of a poro-elastic seabed, after the
// seabed model of Titi Sui (Hohai University). Biot's consolidation equations (u-p form) in the
// cells below the sediment bed (topo < 0, outside solids), one-way coupled to the CFD: the dynamic
// fluid pressure at the bed (CFD pressure minus its value at the start time S 213) is the pore pressure
// at the seabed surface and the load on the soil skeleton. The bed may move (S 10): the soil region
// follows it.
//
//   equilibrium   div sigma' - grad p = rho d2u/dt2               (S 200 2; S 200 1 quasi-static: 0)
//   storage       n beta dp/dt + d(div u)/dt = div(k/gamma_w grad p)
//
//   sigma' cross-anisotropic (z the axis of symmetry), Jeng (1997):
//     sxx = Ev (C11 exx + C12 eyy + C13 ezz), syy = Ev (C12 exx + C11 eyy + C13 ezz),
//     szz = Ev (C13 exx + C13 eyy + C33 ezz), txy = Ev G gxy, txz = Ev L gxz, tyz = Ev L gyz
//   beta = 1/Kw + (1-Sr)/Pw0; extension and tension positive, p the excess pore pressure
//
// Discretisation: the staggered grid of the CFD (displacements at the cell faces, p at the cell
// centres). The elastic operator is the gradient of the discrete strain energy (symmetric), the
// coupling term its exact transpose, variable-step BDF2 in time. Both equations are solved together in
// every step: the symmetric saddle-point system
//
//      | A + c0 rho M     -B^T          | |u|   |f|
//      | -B               -(S + T/a0)   | |p| = |g|
//
// with MINRES and a block-diagonal preconditioner: per rank the Cholesky factors of A + c0 rho M and of
// S + T/a0 + V/(Ev C33) (for B A^-1 B^T), block Jacobi over the ranks. The model has its own fields and a
// one-layer halo of its own (edges and corners included); nothing of the fluid solver is changed.
//
// Boundaries: bed surface (cell face to a fluid cell): p = p_b, sigma'.n = 0; solids (structures, e.g. a
// pipeline) and the bottom of the domain: u = 0, impermeable; lateral domain edges S 211 (1 roller:
// u.n = 0, no shear; 2 fixed: u = 0), impermeable.

#ifndef SOIL_STRESS_H_
#define SOIL_STRESS_H_

#include<vector>
#include<fstream>
#include<mpi.h>

class lexer;
class fdm;
class ghostcell;

class soil_stress
{
public:
    soil_stress(lexer*, fdm*, ghostcell*);
    virtual ~soil_stress();

    // at the new time level of a CFD step (after simtime += dt)
    void step(lexer*, fdm*, ghostcell*);

    enum celltype : int {NOCELL=0, SOIL=1, FREE=2, WALL=3, ROLLER=4};
    enum dofkind : int {NODOF=0, UNK=1, ZERO=2};

private:
    // a vector of the system: ux, uy, uz at the faces, p at the cell centres (4 blocks of N);
    // tc >= 0: the unit vector of the unknown tc (block), tq (index), for the diagonal
    struct vget
    {
        const double *x;
        int tc, tq;
    };

    // set-up (soil_stress.cpp)
    void ini(lexer*, fdm*, ghostcell*);
    void classify(lexer*, fdm*, ghostcell*);
    void bed_pressure(lexer*, fdm*, ghostcell*);
    void halo(ghostcell*, double*);
    void halo_all(ghostcell*, std::vector<double>&);

    // solver (soil_stress_solve.cpp)
    void solve(lexer*, ghostcell*);
    void rhs(std::vector<double>&) const;
    void apply(ghostcell*, std::vector<double>&, std::vector<double>&);
    void precon(const std::vector<double>&, std::vector<double>&) const;
    double dot(ghostcell*, const std::vector<double>&, const std::vector<double>&) const;
    void diagonal();
    void factorize();
    double mech_row(int, int,int,int, const vget&) const;
    double flow_row(int,int,int, const vget&) const;
    double flow_trans(int,int,int, int&) const;
    double strain_vol(int,int,int, const vget&) const;
    void normal_stress(int,int,int, const vget&, double&, double&, double&) const;
    double edge_xz(int,int,int, const vget&, double&, double&) const;
    double edge_yz(int,int,int, const vget&, double&, double&) const;
    double edge_xy(int,int,int, const vget&, double&, double&) const;
    double face_mass(int, int,int,int) const;

    // output (soil_stress_print.cpp)
    void postprocess(lexer*);
    void print_vtu(lexer*);
    void print_probes(lexer*, ghostcell*);

    // the arrays have a halo of one cell: i = -1 .. knox, j = -1 .. knoy, k = -1 .. knoz
    inline int id(int i, int j, int k) const {return ((i+1)*NY + (j+1))*NZ + (k+1);}
    inline int type(int i, int j, int k) const {return int(ctype[id(i,j,k)]);}
    inline double U(const vget &g, int c, int i, int j, int k) const
    {
        const int q = id(i,j,k);
        if(int(kind[c][q])!=UNK)
        return 0.0;
        if(g.tc>=0)
        return (g.tc==c && g.tq==q) ? 1.0 : 0.0;
        return g.x[c*N+q];
    }
    inline double P(const vget &g, int i, int j, int k) const
    {
        const int q = id(i,j,k);
        if(int(ctype[q])!=SOIL)
        return 0.0;
        if(g.tc>=0)
        return (g.tc==3 && g.tq==q) ? 1.0 : 0.0;
        return g.x[3*N+q];
    }
    double DXN(int) const;
    double DYN(int) const;
    double DZN(int) const;
    double DXP(int) const;
    double DYP(int) const;
    double DZP(int) const;

    lexer *p;
    int knox,knoy,knoz,NY,NZ,N;
    bool dim3;
    bool active=false;
    int model;

    // material
    double Ev,C11,C12,C13,C33,Gh,Gv,Kx,Ky,Kz,rho,rhof,grav,gsub,K0,patm;

    // state
    std::vector<double> ctype;              // cell type (celltype)
    std::vector<double> kind[3];            // dof kind of the faces i+1/2, j+1/2, k+1/2 (dofkind)
    std::vector<double> X, Xn, Xnn;         // ux, uy, uz, p now, at n and n-1
    std::vector<double> pb;                 // pore pressure at the bed: at the fluid cells
    std::vector<double> pref;               // CFD pressure at the start time
    std::vector<double> stor;               // storage n beta of the cells
    std::vector<double> D;                  // diagonal (fallback preconditioner)
    struct factor;                          // the block preconditioner: Cholesky of A and of S + T/a0 + beta (Eigen)
    factor *fac = nullptr;
    std::vector<int> own[4];                // owned unknowns: faces and soil cells of the rank

    // MINRES vectors
    std::vector<double> b, r1, r2, y, v, w, w1, w2;

    // time
    double t_last, dt, dt_old;
    double a0,a1,a2;    // BDF2 first derivative
    double c0,c1,c2;    // second derivative (S 200 2)
    int nsteps, nsolve;

    // statistics
    int iter, iter_max, soilcells;
    double resid, soiltime;

    // output
    std::vector<double> post;               // cell outputs (postprocess)
    static constexpr int npost = 14;
    double printtime;
    int printcount;
    std::ofstream probefile;
    MPI_Comm colcomm;                       // the ranks of a column (liquefaction potential)
};

#endif
