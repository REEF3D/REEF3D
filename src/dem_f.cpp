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

#include"dem_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4a.h"
#include<fstream>
#include<sstream>
#include<map>
#include<random>
#include<sys/stat.h>
#include<mpi.h>

dem_f::dem_f(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"DEM startup ..."<<endl;

    solver = 0;
    coupling = p->E11;
    initialized = false;
    nb = 0;

    infile = "dem.txt";
    node_spacing_factor = 0.15;
    sdf_res = 32;
    quad_res = 6;
    gridwalls = 1;

    hybrid_ratio = p->E12;
    Ca = p->E17;
    hs_factor = std::max(p->E18,0.1);
    travel = std::max(p->E20,0.01);
    fluidacc = p->E22;
    kernel_cells = p->E23;
    minsub = std::max(1,p->E16);
    nsub = 1;

    core.maxiter = p->E13;
    core.tol = p->E14;
    core.beta = p->E19;
    core.manifold = p->E21;

    Sx = nullptr;
    Sy = nullptr;
    Sz = nullptr;
    ALPHA = nullptr;
    SX = SY = SZ = ALPHAV = nullptr;

    dt_old = 0.0;
    printtime = 0.0;
    printcount = 0;
    steptime = 0.0;
    maxiter_used = 0;

    read_input(p,pgc);
}

dem_f::~dem_f()
{
    delete Sx;
    delete Sy;
    delete Sz;
    delete ALPHA;
    delete [] SX;
    delete [] SY;
    delete [] SZ;
    delete [] ALPHAV;
}

// ---------------------------------------------------------------------------------------------
// input
// ---------------------------------------------------------------------------------------------

void dem_f::read_input(lexer *p, ghostcell *pgc)
{
    ifstream in(infile.c_str());

    if(!in.is_open())
    {
        if(p->mpirank==0)
        cout<<"DEM: input file "<<infile<<" not found"<<endl;
        pgc->final(true);
    }

    struct shapespec {string type; vector<double> par; string file;};
    map<int,shapespec> shapedef;
    map<int,int> shapeindex, matindex;
    vector<dem_body> bodies;
    vector<int> bshape, bmat;

    mt19937 rng(12345);
    uniform_real_distribution<double> uni(0.0,1.0);

    auto rotation = [](double rx, double ry, double rz)
    {
        const double deg = 3.14159265358979323846/180.0;
        return dem_quat(Eigen::AngleAxisd(rz*deg,dem_vec::UnitZ())*Eigen::AngleAxisd(ry*deg,dem_vec::UnitY())*Eigen::AngleAxisd(rx*deg,dem_vec::UnitX()));
    };

    auto random_rotation = [&]()
    {
        // uniform random unit quaternion (Shoemake)
        double u1=uni(rng), u2=uni(rng), u3=uni(rng);
        const double tp = 2.0*3.14159265358979323846;
        dem_quat q(sqrt(u1)*cos(tp*u3), sqrt(1.0-u1)*sin(tp*u2), sqrt(1.0-u1)*cos(tp*u2), sqrt(u1)*sin(tp*u3));
        q.normalize();
        return q;
    };

    string line;
    int lineno=0;
    while(getline(in,line))
    {
        ++lineno;
        size_t hash = line.find('#');
        if(hash!=string::npos)
        line = line.substr(0,hash);

        istringstream ss(line);
        string key;
        if(!(ss>>key))
        continue;

        if(key=="material")
        {
            int id;
            dem_material m;
            ss>>id>>m.rho>>m.friction>>m.restitution;
            matindex[id] = core.mats.size();
            core.mats.push_back(m);
        }
        else if(key=="wall_material")
        ss>>core.wallmat.friction>>core.wallmat.restitution;

        else if(key=="shape")
        {
            int id;
            shapespec s;
            ss>>id>>s.type;
            if(s.type=="stl")
            {
                double scale=1.0;
                ss>>s.file;
                if(!(ss>>scale)) scale=1.0;
                s.par.push_back(scale);
            }
            else
            {
                double v;
                while(ss>>v)
                s.par.push_back(v);
            }
            shapedef[id] = s;
        }
        else if(key=="node_spacing")
        ss>>node_spacing_factor;

        else if(key=="sdf_resolution")
        ss>>sdf_res;

        else if(key=="quadrature")
        ss>>quad_res;

        else if(key=="grid_walls")
        ss>>gridwalls;

        else if(key=="wall")
        {
            dem_vec n,x;
            ss>>n(0)>>n(1)>>n(2)>>x(0)>>x(1)>>x(2);
            dem_plane pl;
            pl.n = n.normalized();
            pl.d = pl.n.dot(x);
            core.planes.push_back(pl);
        }
        else if(key=="particle" || key=="fixed")
        {
            int sid, mid;
            dem_body B;
            double rx=0,ry=0,rz=0;
            ss>>sid>>mid>>B.x(0)>>B.x(1)>>B.x(2);
            if(ss>>rx>>ry>>rz)
            {
                if(key=="particle")
                ss>>B.v(0)>>B.v(1)>>B.v(2);
            }
            B.q = rotation(rx,ry,rz);
            B.fixed = (key=="fixed");
            bodies.push_back(B);
            bshape.push_back(sid);
            bmat.push_back(mid);
        }
        else if(key=="block")
        {
            int sid, mid, nx, ny, nz, randrot=0;
            unsigned int seed=1;
            dem_vec x0, dx;
            ss>>sid>>mid>>nx>>ny>>nz>>x0(0)>>x0(1)>>x0(2)>>dx(0)>>dx(1)>>dx(2);
            if(ss>>randrot)
            {
                if(ss>>seed)
                rng.seed(seed);
            }

            for(int k=0; k<nz; ++k)
            for(int j=0; j<ny; ++j)
            for(int i=0; i<nx; ++i)
            {
                dem_body B;
                B.x = x0 + dem_vec(i*dx(0),j*dx(1),k*dx(2));
                B.q = randrot ? random_rotation() : dem_quat::Identity();
                bodies.push_back(B);
                bshape.push_back(sid);
                bmat.push_back(mid);
            }
        }
        else if(p->mpirank==0)
        cout<<"DEM: unknown keyword '"<<key<<"' in line "<<lineno<<" of "<<infile<<endl;
    }

    if(core.mats.empty())
    {
        matindex[0]=0;
        core.mats.push_back(dem_material());
    }

    // build shapes
    for(auto &sd : shapedef)
    {
        const shapespec &s = sd.second;
        dem_shape S;
        bool ok=true;

        if(s.type=="sphere" && s.par.size()>=1)
        S.build_sphere(s.par[0]);

        else if(s.type=="box" && s.par.size()>=3)
        S.build_box(s.par[0],s.par[1],s.par[2]);

        else if(s.type=="cylinder" && s.par.size()>=2)
        S.build_cylinder(s.par[0],s.par[1]);

        else if(s.type=="ellipsoid" && s.par.size()>=3)
        S.build_ellipsoid(s.par[0],s.par[1],s.par[2],sdf_res);

        else if(s.type=="stl")
        ok = S.build_mesh(s.file,s.par[0],sdf_res);

        else
        ok=false;

        if(!ok)
        {
            if(p->mpirank==0)
            cout<<"DEM: invalid shape "<<sd.first<<" ("<<s.type<<")"<<endl;
            pgc->final(true);
        }

        S.make_nodes(node_spacing_factor*S.deq);
        S.make_quadrature(quad_res);

        shapeindex[sd.first] = core.shapes.size();
        core.shapes.push_back(S);
    }

    for(size_t n=0; n<bodies.size(); ++n)
    {
        if(shapeindex.find(bshape[n])==shapeindex.end() || matindex.find(bmat[n])==matindex.end())
        {
            if(p->mpirank==0)
            cout<<"DEM: particle "<<n<<" refers to an undefined shape or material"<<endl;
            pgc->final(true);
        }
        bodies[n].shape = shapeindex[bshape[n]];
        bodies[n].mat = matindex[bmat[n]];
        bodies[n].id = n;
    }

    core.bodies = bodies;
}

// ---------------------------------------------------------------------------------------------
// initialisation
// ---------------------------------------------------------------------------------------------

void dem_f::ini(lexer *p, ghostcell *pgc)
{
    // grid length scale: CFD mean cell size; NHFLOW horizontal cell size only
    if(solver==6)
    dxs = p->DXM;
    else
    {
        double sum=0.0, cnt=0.0;
        for(i=0; i<p->knox; ++i)
        {
            sum += p->DXN[IP];
            cnt += 1.0;
        }
        if(p->j_dir==1)
        for(j=0; j<p->knoy; ++j)
        {
            sum += p->DYN[JP];
            cnt += 1.0;
        }
        sum = pgc->globalsum(sum);
        cnt = pgc->globalsum(cnt);
        dxs = cnt>0.0 ? sum/cnt : p->DXM;
    }

    core.gravity = dem_vec(p->W20,p->W21,p->W22);
    core.plane2D = (p->j_dir==0);
    core.initialize();

    nb = core.bodies.size();

    rhof.assign(nb,0.0);
    nuf.assign(nb,0.0);
    epsf.assign(nb,1.0);
    vsub.assign(nb,0.0);
    dvol.assign(nb,0.0);
    ufl.assign(nb,dem_vec::Zero());
    ufl_old.assign(nb,dem_vec::Zero());
    Fb.assign(nb,dem_vec::Zero());
    Tb.assign(nb,dem_vec::Zero());
    fluidcount.assign(nb,0);
    ufl_valid.assign(nb,false);
    Fibm.assign(nb,dem_vec::Zero());
    Fstage.assign(3*nb,dem_vec::Zero());
    Tstage.assign(3*nb,dem_vec::Zero());
    Tibm.assign(nb,dem_vec::Zero());
    mfl.assign(nb,0.0);
    hvol.assign(nb,0.0);
    Ffp.assign(nb,dem_vec::Zero());
    Fhyd.assign(nb,dem_vec::Zero());
    vprev.assign(nb,dem_vec::Zero());
    Ifl.assign(nb,dem_vec::Zero());
    Lfl.assign(nb,dem_vec::Zero());
    Ifl_old.assign(nb,dem_vec::Zero());
    Lfl_old.assign(nb,dem_vec::Zero());
    Ifl_valid.assign(nb,false);
    wprev.assign(nb,dem_vec::Zero());

    for(int n=0; n<nb; ++n)
    {
        vprev[n] = core.bodies[n].v;
        wprev[n] = core.bodies[n].w;
    }

    // coupling mode per particle
    int nres=0;
    for(auto &B : core.bodies)
    {
        double ratio = core.shapes[B.shape].deq/dxs;

        if(coupling==2)
        B.mode = 1;
        else if(coupling==3)
        B.mode = ratio>=hybrid_ratio ? 1 : 0;
        else
        B.mode = 0;

        nres += B.mode;
    }

    basemode.assign(nb,0);
    for(int n=0; n<nb; ++n)
    basemode[n] = core.bodies[n].mode;

    if(p->mpirank==0)
    {
        mkdir("./REEF3D_DEM",0777);
        if(p->E15>0.0)
        mkdir("./REEF3D_DEM_VTP",0777);

        cout<<"DEM: "<<core.shapes.size()<<" shapes, "<<core.mats.size()<<" materials, "<<nb<<" particles, "<<core.planes.size()<<" walls"<<endl;
        for(size_t s=0; s<core.shapes.size(); ++s)
        {
            const dem_shape &S = core.shapes[s];
            cout<<"DEM: shape "<<s<<" type "<<S.type<<" V: "<<S.volume<<" d_eq: "<<S.deq<<" r_bound: "<<S.rbound
                <<" sphericity: "<<S.sphericity<<" nodes: "<<S.nodes.size()<<" quadrature points: "<<S.qp.size()<<endl;
        }
        cout<<"DEM: coupling "<<coupling<<", resolved particles: "<<nres<<", unresolved: "<<nb-nres<<endl;
    }

    printtime = p->simtime;
}

void dem_f::ini_cfd(lexer *p, ghostcell *pgc)
{
    Sx = new field1(p);
    Sy = new field2(p);
    Sz = new field3(p);
    ALPHA = new field4a(p);

    ULOOP
    (*Sx)(i,j,k) = 0.0;
    VLOOP
    (*Sy)(i,j,k) = 0.0;
    WLOOP
    (*Sz)(i,j,k) = 0.0;
    LOOP
    (*ALPHA)(i,j,k) = 0.0;
}

void dem_f::ini_nhflow(lexer *p, ghostcell *pgc)
{
    p->Darray(SX,p->imax*p->jmax*(p->kmax+2));
    p->Darray(SY,p->imax*p->jmax*(p->kmax+2));
    p->Darray(SZ,p->imax*p->jmax*(p->kmax+2));
    p->Darray(ALPHAV,p->imax*p->jmax*(p->kmax+2));
}

// ---------------------------------------------------------------------------------------------
// driver entry points
// ---------------------------------------------------------------------------------------------

void dem_f::start_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(!initialized)
    {
        solver = 6;
        ini(p,pgc);
        ini_cfd(p,pgc);
        initialized = true;
        print(p,pgc);
    }

    double starttime = pgc->timer();

    fluid_cfd(p,a,pgc);
    internal_cfd(p,a,pgc);
    set_loads(p,pgc);

    core.gravity = dem_vec(p->W20,p->W21,p->W22);
    dem_core::wallfunc wf = nullptr;
    if(gridwalls==1 && (p->toporead>0 || p->S10>0 || p->solidread==1))
    wf = [&](dem_core &c, double margin, vector<dem_contact> &cts) {walls_cfd(p,a,pgc,margin,cts);};

    nsub = substeps(p);
    maxiter_used = 0;
    for(int s=0; s<nsub; ++s)
    {
        core.step(p->dt/double(nsub),wf);
        maxiter_used = std::max(maxiter_used,core.iterations);
    }

    // hydrodynamic force over the step (drag and added mass evaluated with the new velocity)
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        dem_vec ap = (B.v - vprev[n])/p->dt;
        Fhyd[n] = B.F + B.K*(B.uf - B.v) + B.madd*(B.af - ap);
    }

    sync(p,pgc);
    deactivate(p);

    feedback_cfd(p,a,pgc);

    for(int n=0; n<nb; ++n)
    {
        Fibm[n].setZero();
        Tibm[n].setZero();
        mfl[n] = 0.0;
        hvol[n] = 0.0;
    }
    dt_old = p->dt;

    steptime = pgc->timer()-starttime;
    print(p,pgc);
}

void dem_f::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!initialized)
    {
        solver = 5;
        ini(p,pgc);
        ini_nhflow(p,pgc);
        initialized = true;
        print(p,pgc);
    }

    double starttime = pgc->timer();

    fluid_nhflow(p,d,pgc);
    internal_nhflow(p,d,pgc);
    set_loads(p,pgc);

    core.gravity = dem_vec(p->W20,p->W21,p->W22);
    dem_core::wallfunc wf = nullptr;
    if(gridwalls==1)
    wf = [&](dem_core &c, double margin, vector<dem_contact> &cts) {walls_nhflow(p,d,pgc,margin,cts);};

    nsub = substeps(p);
    maxiter_used = 0;
    for(int s=0; s<nsub; ++s)
    {
        core.step(p->dt/double(nsub),wf);
        maxiter_used = std::max(maxiter_used,core.iterations);
    }

    // hydrodynamic force over the step (drag and added mass evaluated with the new velocity)
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        dem_vec ap = (B.v - vprev[n])/p->dt;
        Fhyd[n] = B.F + B.K*(B.uf - B.v) + B.madd*(B.af - ap);
    }

    sync(p,pgc);
    deactivate(p);

    feedback_nhflow(p,d,pgc);

    for(int n=0; n<nb; ++n)
    {
        Fibm[n].setZero();
        Tibm[n].setZero();
        mfl[n] = 0.0;
        hvol[n] = 0.0;
    }
    dt_old = p->dt;

    steptime = pgc->timer()-starttime;
    print(p,pgc);
}

// ---------------------------------------------------------------------------------------------
// loads
// ---------------------------------------------------------------------------------------------

void dem_f::set_loads(lexer *p, ghostcell *pgc)
{
    dem_vec g(p->W20,p->W21,p->W22);

    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        const dem_shape &S = core.shapes[B.shape];

        B.F.setZero();
        B.T.setZero();
        B.K = 0.0;
        B.Kr = 0.0;
        B.madd = 0.0;
        B.uf.setZero();
        B.af.setZero();

        // lagged particle acceleration over the last fluid step
        vprev[n] = B.v;
        wprev[n] = B.w;

        if(!B.active || B.fixed || coupling==0)
        continue;

        if(B.mode==0)
        {
            // unresolved
            if(fluidcount[n]==0 || rhof[n]<=0.0)
            continue;

            double fsub = solver==5 ? std::min(1.0,vsub[n]/S.volume) : 1.0;

            B.F += Fb[n];
            B.T += Tb[n];

            if(fsub>0.0)
            {
                drag(p,n,rhof[n],nuf[n],epsf[n]);
                B.K *= fsub;
                B.Kr *= fsub;

                // added mass; the fluid acceleration terms use the local Eulerian acceleration, which
                // contains the particle's self-induced flow, hence optional (E 22)
                B.madd = Ca*rhof[n]*S.volume*fsub;
                if(fluidacc==1 && ufl_valid[n] && dt_old>0.0)
                {
                    dem_vec af = (ufl[n]-ufl_old[n])/dt_old;
                    B.af = af;
                    B.F += rhof[n]*S.volume*fsub*af;
                }
            }
        }
        else
        {
            // resolved: forcing integral, rate of change of the fluid momentum inside the particle
            // (evaluated from the fluid field, Kempe & Froehlich 2012), buoyancy
            B.F = Fibm[n];
            B.T = Tibm[n];

            if(Ifl_valid[n] && dt_old>0.0)
            {
                B.F += (Ifl[n] - Ifl_old[n])/dt_old;
                B.T += (Lfl[n] - Lfl_old[n])/dt_old;
            }

            if(solver==6)
            B.F -= mfl[n]*g;         // gravity acts on the forced fluid in CFD, Uhlmann (2005)


            if(solver==5)
            {
                B.F += Fb[n];        // hydrostatic pressure is not part of the NHFLOW pressure
                B.T += Tb[n];
            }
        }
    }
}

void dem_f::drag(lexer *p, int n, double rho, double nu, double eps)
{
    dem_body &B = core.bodies[n];
    const dem_shape &S = core.shapes[B.shape];

    const double pi = 3.14159265358979323846;
    double d = S.deq;
    double A = 0.25*pi*d*d;
    double phi = std::max(0.1,std::min(1.0,S.sphericity));
    eps = std::max(0.2,std::min(1.0,eps));
    nu = std::max(nu,1.0e-12);

    double ur = (ufl[n]-B.v).norm();
    double Re = eps*ur*d/nu;

    // Haider & Levenspiel (1989), non-spherical particles; Cd*|u_rel|
    double cdu;
    if(Re<1.0e-8)
    cdu = 24.0*nu/(eps*d);
    else
    {
        double a1 = exp(2.3288 - 6.4581*phi + 2.4486*phi*phi);
        double b1 = 0.0964 + 0.5565*phi;
        double c1 = exp(4.905 - 13.8944*phi + 18.4222*phi*phi - 10.2599*phi*phi*phi);
        double d1 = exp(1.4681 + 12.2584*phi - 20.7322*phi*phi + 15.8855*phi*phi*phi);
        double cd = 24.0/Re*(1.0 + a1*pow(Re,b1)) + c1/(1.0 + d1/Re);
        cdu = cd*ur;
    }

    // Di Felice (1994) voidage function
    double chi = 3.7;
    if(Re>1.0e-8)
    chi = 3.7 - 0.65*exp(-0.5*pow(1.5 - log10(Re),2.0));

    B.K  = 0.5*rho*A*pow(eps,2.0-chi)*cdu;
    B.Kr = pi*rho*nu*d*d*d;
    B.uf = ufl[n];
}

// ---------------------------------------------------------------------------------------------
// parallel helpers
// ---------------------------------------------------------------------------------------------

bool dem_f::owns(lexer *p, double x, double y, double z)
{
    auto in = [](double v, double lo, double hi, double gmax)
    {
        return v>=lo && (v<hi || (hi>=gmax-1.0e-12 && v<=hi));
    };

    if(!in(x,p->originx,p->endx,p->global_xmax))
    return false;

    if(p->j_dir==1 && !in(y,p->originy,p->endy,p->global_ymax))
    return false;

    if(solver==6 && !in(z,p->originz,p->endz,p->global_zmax))
    return false;

    return true;
}

void dem_f::sync(lexer *p, ghostcell *pgc)
{
    const int nv = 14;
    vector<double> buf(nb*nv);

    if(p->mpirank==0)
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        double *b = &buf[n*nv];
        b[0]=B.x(0); b[1]=B.x(1); b[2]=B.x(2);
        b[3]=B.q.w(); b[4]=B.q.x(); b[5]=B.q.y(); b[6]=B.q.z();
        b[7]=B.v(0); b[8]=B.v(1); b[9]=B.v(2);
        b[10]=B.w(0); b[11]=B.w(1); b[12]=B.w(2);
        b[13]=B.active ? 1.0 : 0.0;
    }

    MPI_Bcast(buf.data(),nb*nv,MPI_DOUBLE,0,pgc->mpi_comm);

    if(p->mpirank>0)
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        const double *b = &buf[n*nv];
        B.x = dem_vec(b[0],b[1],b[2]);
        B.q = dem_quat(b[3],b[4],b[5],b[6]);
        B.v = dem_vec(b[7],b[8],b[9]);
        B.w = dem_vec(b[10],b[11],b[12]);
        B.active = b[13]>0.5;
        B.R = B.q.toRotationMatrix();
    }
}

void dem_f::deactivate(lexer *p)
{
    for(auto &B : core.bodies)
    {
        if(!B.active)
        continue;

        double r = 2.0*core.shapes[B.shape].rbound;
        bool out = B.x(0)<p->global_xmin-r || B.x(0)>p->global_xmax+r;
        if(solver==6)
        out = out || B.x(2)<p->global_zmin-r || B.x(2)>p->global_zmax+r;
        if(p->j_dir==1)
        out = out || B.x(1)<p->global_ymin-r || B.x(1)>p->global_ymax+r;

        if(out)
        {
            B.active = false;
            B.v.setZero();
            B.w.setZero();
            if(p->mpirank==0)
            cout<<"DEM: particle "<<B.id<<" left the domain and is deactivated"<<endl;
        }
    }
}

void dem_f::gather_walls(lexer *p, ghostcell *pgc, const vector<double> &loc, vector<dem_contact> &cts)
{
    // packed: body, feature, x(3), n(3), gap
    const int nv = 9;
    int nloc = loc.size();
    vector<int> counts(p->mpi_size), displs(p->mpi_size,0);
    MPI_Allgather(&nloc,1,MPI_INT,counts.data(),1,MPI_INT,pgc->mpi_comm);

    int ntot=0;
    for(int q=0; q<p->mpi_size; ++q)
    {
        displs[q] = ntot;
        ntot += counts[q];
    }

    vector<double> all(std::max(ntot,1));
    MPI_Allgatherv(loc.data(),nloc,MPI_DOUBLE,all.data(),counts.data(),displs.data(),MPI_DOUBLE,pgc->mpi_comm);

    for(int c=0; c<ntot/nv; ++c)
    {
        const double *b = &all[c*nv];
        dem_contact C;
        C.a = int(b[0]);
        C.b = -1;
        int feature = int(b[1]);
        C.key = core.make_key(C.a,-1,feature);
        C.x = dem_vec(b[2],b[3],b[4]);
        C.n = dem_vec(b[5],b[6],b[7]);
        C.gap = b[8];
        const dem_material &M = core.mats[core.bodies[C.a].mat];
        C.mu = 0.5*(M.friction + core.wallmat.friction);
        C.e  = 0.5*(M.restitution + core.wallmat.restitution);
        cts.push_back(C);
    }
}

double dem_f::heaviside(double phi, double eps)
{
    if(phi<-eps)
    return 1.0;
    if(phi>eps)
    return 0.0;
    return 0.5*(1.0 - phi/eps - sin(3.14159265358979323846*phi/eps)/3.14159265358979323846);
}

double dem_f::kernel(double r, double R)
{
    if(r>=R)
    return 0.0;
    double s = 1.0 - (r*r)/(R*R);
    return s*s;
}

void dem_f::cellrange(lexer *p, int n, double R, int &i0, int &i1, int &j0, int &j1, int &k0, int &k1)
{
    const dem_vec &x = core.bodies[n].x;

    i0=0; i1=-1; j0=0; j1=-1; k0=0; k1=-1;

    if(x(0)+R<p->originx || x(0)-R>p->endx)
    return;
    if(p->j_dir==1 && (x(1)+R<p->originy || x(1)-R>p->endy))
    return;
    if(solver==6 && (x(2)+R<p->originz || x(2)-R>p->endz))
    return;

    i0 = std::max(0,std::min(p->knox-1,p->posc_i(x(0)-R)-1));
    i1 = std::max(0,std::min(p->knox-1,p->posc_i(x(0)+R)+1));

    if(p->j_dir==1)
    {
        j0 = std::max(0,std::min(p->knoy-1,p->posc_j(x(1)-R)-1));
        j1 = std::max(0,std::min(p->knoy-1,p->posc_j(x(1)+R)+1));
    }
    else
    {
        j0=0;
        j1=0;
    }

    if(solver==6)
    {
        k0 = std::max(0,std::min(p->knoz-1,p->posc_k(x(2)-R)-1));
        k1 = std::max(0,std::min(p->knoz-1,p->posc_k(x(2)+R)+1));
    }
    else
    {
        k0 = 0;
        k1 = p->knoz-1;
    }
}

void dem_f::combine_stages(lexer *p, ghostcell *pgc, int iter, double alpha)
{
    // time average of the stage forcing with the Runge-Kutta weights:
    // SSP-RK3 (alpha 1, 1/4, 2/3): 1/6, 1/6, 2/3; RK2 (alpha 1, 1/2): 1/2, 1/2; otherwise final stage only
    double wgt[3] = {0.0,0.0,0.0};

    if(iter==2 && fabs(alpha-2.0/3.0)<1.0e-8)
    {
        wgt[0]=1.0/6.0;
        wgt[1]=1.0/6.0;
        wgt[2]=2.0/3.0;
    }
    else if(iter==2 && fabs(alpha-1.0/3.0)<1.0e-8)
    {
        // low-storage RK3 (N 40 4, alpha 8/15, 2/15, 1/3): additive stages
        wgt[0]=8.0/15.0;
        wgt[1]=2.0/15.0;
        wgt[2]=1.0/3.0;
    }
    else if(iter==1 && fabs(alpha-0.5)<1.0e-8)
    {
        wgt[0]=0.5;
        wgt[1]=0.5;
    }
    else if(iter==0)
    wgt[0]=1.0;
    else
    {
        wgt[std::min(std::max(iter,0),2)]=1.0;
        if(p->mpirank==0 && !rkwarn)
        cout<<"DEM: warning, unknown Runge-Kutta stage weights (final stage "<<iter<<", alpha "<<alpha<<"), resolved forces use the final stage only"<<endl;
        rkwarn=true;
    }

    vector<double> buf(8*nb);
    for(int q=0; q<nb; ++q)
    {
        dem_vec F = dem_vec::Zero(), T = dem_vec::Zero();
        for(int s=0; s<3; ++s)
        {
            F += wgt[s]*Fstage[3*q+s];
            T += wgt[s]*Tstage[3*q+s];
            Fstage[3*q+s].setZero();
            Tstage[3*q+s].setZero();
        }
        buf[8*q+0]=F(0); buf[8*q+1]=F(1); buf[8*q+2]=F(2);
        buf[8*q+3]=T(0); buf[8*q+4]=T(1); buf[8*q+5]=T(2);
        buf[8*q+6]=mfl[q];
        buf[8*q+7]=hvol[q];
    }

    MPI_Allreduce(MPI_IN_PLACE,buf.data(),8*nb,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    for(int q=0; q<nb; ++q)
    {
        Fibm[q] = dem_vec(buf[8*q+0],buf[8*q+1],buf[8*q+2]);
        Tibm[q] = dem_vec(buf[8*q+3],buf[8*q+4],buf[8*q+5]);

        // fluid mass of the particle volume: mean density over the smoothed indicator times the exact
        // volume (the smoothed indicator overestimates the volume by O(eps^2/r)); NHFLOW: submerged volume
        const dem_body &B = core.bodies[q];
        double vol = solver==5 ? vsub[q] : core.shapes[B.shape].volume;
        mfl[q] = buf[8*q+7]>0.0 ? buf[8*q+6]/buf[8*q+7]*vol : 0.0;
    }
}

int dem_f::substeps(lexer *p)
{
    // particles must not travel more than a fraction of their size per DEM step (contact detection)
    double ns = ceil(core.max_velocity()*p->dt/(travel*core.min_rbound()));
    int n = std::max(minsub,int(std::min(ns,1.0e6)));

    if(n>1000)
    {
        if(p->mpirank==0)
        cout<<"DEM: warning, "<<n<<" substeps required, limited to 1000 (max particle velocity "<<core.max_velocity()<<" m/s)"<<endl;
        n=1000;
    }
    return n;
}

dem_vec dem_f::relpos(lexer *p, double x, double y, double z, const dem_vec &c)
{
    return dem_vec(x-c(0), p->j_dir==1 ? y-c(1) : 0.0, z-c(2));
}
