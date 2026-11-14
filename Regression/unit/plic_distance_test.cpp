// Standalone check of plic_distance (the printed phi of a VOF-PLIC run, F 80 4): no MPI, no
// REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src plic_distance_test.cpp -o plic_distance_test
// Run:    ./plic_distance_test
// Checks: the plane constant holds the cell's fraction under the plane (sampled, random normals,
//         stretched cells); with the vof of a flat surface, of a tilted plane and of a surface
//         exactly on a cell face, the distance at the cell centres is that to the plane in the
//         band and the vertical distance to the surface outside it, with the sign of the water.
// Architect: Hans Bihs

// plic_distance.cpp needs only these members of the lexer, a field, a ghost cell exchange and
// increment::marge: small stand-ins, so the test builds without the solver (the include guards
// keep the real headers out)
#include<cmath>
#include<vector>
#define LEXER_H_
#define GHOSTCELL_H_
#define INCREMENT_H_
#define FIELD_H_
#define FIELD4_H_

class increment
{
public:
    static int marge;
};
int increment::marge = 5;

class lexer
{
public:
    int knox=0, knoy=1, knoz=0, j_dir=0;
    int imin=-3, jmin=-3, kmin=-3, imax=0, jmax=0, kmax=0;
    double F93=0.001, F94=0.999;
    std::vector<int> flag;
    int *flag4=nullptr;
    std::vector<double> XP, YP, ZP, DXN, DYN, DZN;

    // a grid of nx x 1 x nz cells (2D) from (0, 0) with cells dx, dz (dz may vary: dzk)
    void grid(int nx, int nz, double dx, const std::vector<double> &dzk, int ny=1)
    {
        knox=nx; knoy=ny; knoz=nz; j_dir = ny>1 ? 1 : 0;
        imax=nx+6; jmax=ny+6; kmax=nz+6;
        flag.assign(imax*jmax*kmax, 1);
        flag4=flag.data();
        const int m=increment::marge;
        XP.assign(nx+2*m, 0.0); DXN.assign(nx+2*m, dx);
        YP.assign(ny+2*m, 0.5); DYN.assign(ny+2*m, ny>1 ? dx : 1.0);
        if(ny>1)
        for(int j=-m; j<ny+m; ++j)
        YP[j+m] = (j+0.5)*dx;
        ZP.assign(nz+2*m, 0.0); DZN.assign(nz+2*m, dzk[0]);
        for(int i=-m; i<nx+m; ++i)
        XP[i+m] = (i+0.5)*dx;
        for(int k=0; k<nz; ++k)
        DZN[k+m] = dzk[k];
        for(int k=-m; k<0; ++k)
        DZN[k+m] = dzk[0];
        for(int k=nz; k<nz+m; ++k)
        DZN[k+m] = dzk[nz-1];
        double z=0.0;
        for(int k=0; k<nz; ++k) { ZP[k+m] = z+0.5*DZN[k+m]; z += DZN[k+m]; }
        z=0.0;
        for(int k=-1; k>=-m; --k) { z -= DZN[k+m]; ZP[k+m] = z+0.5*DZN[k+m]; }
        z=ZP[nz-1+m]+0.5*DZN[nz-1+m];
        for(int k=nz; k<nz+m; ++k) { ZP[k+m] = z+0.5*DZN[k+m]; z += DZN[k+m]; }
    }
};

class field
{
public:
    explicit field(lexer *p) : p(p), V(p->imax*p->jmax*p->kmax, 0.0) {}
    virtual ~field() {}
    double& operator()(int i, int j, int k)
    {
        return V[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
    }
    lexer *p;
    std::vector<double> V;
};

class field4 : public field
{
public:
    explicit field4(lexer *p) : field(p) {}
};

class ghostcell
{
public:
    // one rank, 2D: the ghost cells in y are the cell itself (as REEF3D's in a 2D run)
    void start4(lexer *p, field &f, int)
    {
        if(p->knoy==1)
        {
        for(int i=-3; i<p->knox+3; ++i)
        for(int k=-3; k<p->knoz+3; ++k)
        for(int j=-3; j<4; ++j)
        if(j!=0)
        f(i,j,k) = f(std::min(std::max(i,0),p->knox-1),0,std::min(std::max(k,0),p->knoz-1));
        }
        // 3D: the tests fill the ghost cells of the vof themselves; the distance is checked inside
    }
};

#include"../../src/plic_distance.cpp"

#include<cstdio>
#include<random>
#include<functional>

static int failed = 0;

static void check(bool ok, const char *what, double got, double want)
{
    if(!ok)
    {
        ++failed;
        std::printf("FAIL %s: %.9g (want %.9g)\n", what, got, want);
    }
}

// the vof of every cell (ghosts too) under the surface z = h(x), sampled
static void fill_vof(lexer &p, field &f, const std::function<double(double)> &h)
{
    const int m=increment::marge, n=200;
    for(int i=-3; i<p.knox+3; ++i)
    for(int k=-3; k<p.knoz+3; ++k)
    {
        const double x0=p.XP[i+m]-0.5*p.DXN[i+m], z0=p.ZP[k+m]-0.5*p.DZN[k+m];
        double v=0.0;
        for(int a=0; a<n; ++a)
        {
            const double x = x0+(a+0.5)/n*p.DXN[i+m];
            v += std::min(std::max((h(x)-z0)/p.DZN[k+m],0.0),1.0);
        }
        for(int j=-3; j<4; ++j)
        f(i,j,k) = v/n;
    }
}

int main()
{
    // 1. the plane constant: the sampled fraction under the plane is V0
    {
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> u(-1.0,1.0), v(0.02,0.98);
        double worst=0.0;
        for(int t=0; t<200; ++t)
        {
            double nx=u(rng), ny=u(rng), nz=u(rng);
            const double s=std::sqrt(nx*nx+ny*ny+nz*nz);
            nx/=s; ny/=s; nz/=s;
            const double dx=0.5+v(rng), dy=0.5+v(rng), dz=0.5+v(rng), V0=v(rng);
            const double a = plic_distance::plane_constant(nx,ny,nz,dx,dy,dz,V0);
            const int n=40;
            int under=0;
            for(int i=0; i<n; ++i)
            for(int j=0; j<n; ++j)
            for(int k=0; k<n; ++k)
            if(nx*((i+0.5)/n-0.5)*dx + ny*((j+0.5)/n-0.5)*dy + nz*((k+0.5)/n-0.5)*dz <= a)
            ++under;
            worst = std::max(worst, std::fabs(double(under)/(n*n*n)-V0));
        }
        check(worst<0.01, "plane constant: largest sampled fraction error", worst, 0.0);
        check(std::fabs(plic_distance::plane_constant(0,0,1,1,1,1,0.3)+0.2)<1e-12, "flat, V0 0.3", plic_distance::plane_constant(0,0,1,1,1,1,0.3), -0.2);
        check(std::fabs(plic_distance::plane_constant(1,0,0,2,1,1,0.7)-0.4)<1e-12, "x normal, dx 2, V0 0.7", plic_distance::plane_constant(1,0,0,2,1,1,0.7), 0.4);
    }

    // 2. a flat surface in a cell, stretched cells in z: the exact distance everywhere
    {
        lexer p;
        std::vector<double> dz(24);
        for(int k=0; k<24; ++k)
        dz[k] = 0.05*(1.0+0.03*k);
        p.grid(16,24,0.05,dz);
        field4 vof(&p);
        const double zs=0.6137;
        fill_vof(p,vof,[zs](double){return zs;});
        ghostcell g;
        plic_distance d(&p);
        field4 &phi = d.start(&p,&g,vof);
        double worst=0.0;
        for(int i=0; i<p.knox; ++i)
        for(int k=0; k<p.knoz; ++k)
        worst = std::max(worst, std::fabs(phi(i,0,k)-(zs-p.ZP[k+5])));
        check(worst<1e-9, "flat surface: largest distance error", worst, 0.0);
    }

    // 3. a tilted plane: the distance to it in the band (Youngs is not exact for a tilted plane:
    //    within 3 % of a cell), the vertical distance to the surface outside it, the sign of the water
    {
        lexer p;
        p.grid(40,30,0.05,std::vector<double>(30,0.05));
        field4 vof(&p);
        const double s=0.3, c=std::cos(std::atan(s));
        auto h = [s](double x){return 0.7+s*(x-1.0);};
        fill_vof(p,vof,h);
        ghostcell g;
        plic_distance d(&p);
        field4 &phi = d.start(&p,&g,vof);
        double worst_band=0.0, worst_far=0.0;
        int sign=0;
        for(int i=3; i<p.knox-3; ++i)
        for(int k=0; k<p.knoz; ++k)
        {
            const double x=p.XP[i+5], z=p.ZP[k+5], vertical=h(x)-z, normal=vertical*c;
            if(std::fabs(normal)<1.5*0.05)
            worst_band = std::max(worst_band, std::fabs(phi(i,0,k)-normal));
            else if(std::fabs(vertical)>3.0*0.05)
            worst_far = std::max(worst_far, std::fabs(phi(i,0,k)-vertical));
            if((phi(i,0,k)>0.0) != (vof(i,0,k)>=0.5))
            ++sign;
        }
        // where phi = 0 lies on each column (linear between the cell centres, as an iso-surface)
        double worst_zero=0.0;
        for(int i=3; i<p.knox-3; ++i)
        for(int k=0; k<p.knoz-1; ++k)
        if((phi(i,0,k)>=0.0) != (phi(i,0,k+1)>=0.0))
        {
            const double z0=p.ZP[k+5], z1=p.ZP[k+6], t=phi(i,0,k)/(phi(i,0,k)-phi(i,0,k+1));
            worst_zero = std::max(worst_zero, std::fabs(z0+t*(z1-z0)-h(p.XP[i+5])));
        }
        std::printf("tilted plane: band %.4f, outside %.4f, zero %.4f cells\n", worst_band/0.05, worst_far/0.05, worst_zero/0.05);
        check(worst_zero<0.05*0.05, "tilted plane: largest error of phi = 0 on the columns", worst_zero, 0.0);
        check(worst_band<0.15*0.05, "tilted plane: largest distance error in the band", worst_band, 0.0);
        check(worst_far<0.25*0.05, "tilted plane: largest vertical distance error outside", worst_far, 0.0);
        check(sign==0, "tilted plane: cells with the wrong sign", sign, 0);
    }

    // 4. the surface exactly on a cell face (no interface cell): the plane of the face
    {
        lexer p;
        p.grid(8,20,0.1,std::vector<double>(20,0.1));
        field4 vof(&p);
        fill_vof(p,vof,[](double){return 1.0;});
        for(int i=-3; i<p.knox+3; ++i)          // no interface cell: only full and empty
        for(int k=-3; k<p.knoz+3; ++k)
        for(int j=-3; j<4; ++j)
        vof(i,j,k) = vof(i,j,k)>0.5 ? 1.0 : 0.0;
        ghostcell g;
        plic_distance d(&p);
        field4 &phi = d.start(&p,&g,vof);
        double worst=0.0;
        for(int i=0; i<p.knox; ++i)
        for(int k=0; k<p.knoz; ++k)
        worst = std::max(worst, std::fabs(phi(i,0,k)-(1.0-p.ZP[k+5])));
        check(worst<1e-9, "surface on a face: largest distance error", worst, 0.0);
        // a step on the face: one column of water higher, the distance to its side face
        for(int k=10; k<12; ++k)
        for(int j=-3; j<4; ++j)
        vof(4,j,k) = 1.0;
        field4 &step = d.start(&p,&g,vof);
        const double side = std::fabs(step(5,0,10));      // the air cell beside the column: half a cell
        check(std::fabs(side-0.05)<1e-9, "step: distance to the side face", side, 0.05);
    }

    // 5. 3D: a plane tilted in x and y, where phi = 0 lies on each column
    {
        lexer p;
        const int nx=24, ny=16, nz=20;
        const double dx=0.05;
        p.grid(nx,nz,dx,std::vector<double>(nz,dx),ny);
        field4 vof(&p);
        auto h = [](double x, double y){return 0.5+0.3*(x-0.6)+0.2*(y-0.4);};
        const int m=increment::marge, n=24;
        for(int i=-3; i<nx+3; ++i)
        for(int j=-3; j<ny+3; ++j)
        for(int k=-3; k<nz+3; ++k)
        {
            const double x0=p.XP[i+m]-0.5*dx, y0=p.YP[j+m]-0.5*dx, z0=p.ZP[k+m]-0.5*dx;
            double v=0.0;
            for(int a=0; a<n; ++a)
            for(int b=0; b<n; ++b)
            v += std::min(std::max((h(x0+(a+0.5)/n*dx,y0+(b+0.5)/n*dx)-z0)/dx,0.0),1.0);
            vof(i,j,k) = v/(n*n);
        }
        ghostcell g;
        plic_distance d(&p);
        field4 &phi = d.start(&p,&g,vof);
        double worst_zero=0.0;
        int sign=0;
        for(int i=3; i<nx-3; ++i)
        for(int j=3; j<ny-3; ++j)
        for(int k=0; k<nz-1; ++k)
        {
            if((phi(i,j,k)>=0.0) != (phi(i,j,k+1)>=0.0))
            {
                const double z0=p.ZP[k+m], z1=p.ZP[k+1+m], t=phi(i,j,k)/(phi(i,j,k)-phi(i,j,k+1));
                worst_zero = std::max(worst_zero, std::fabs(z0+t*(z1-z0)-h(p.XP[i+m],p.YP[j+m])));
            }
            if((phi(i,j,k)>0.0) != (vof(i,j,k)>=0.5))
            ++sign;
        }
        std::printf("3D tilted plane: zero %.4f cells\n", worst_zero/dx);
        check(worst_zero<0.05*dx, "3D tilted plane: largest error of phi = 0 on the columns", worst_zero, 0.0);
        check(sign==0, "3D tilted plane: cells with the wrong sign", sign, 0);
    }

    if(failed==0)
    std::printf("plic_distance_test: all checks passed\n");
    return failed==0 ? 0 : 1;
}
