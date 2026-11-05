// Standalone verification of the LAGOON body writer (lagoon_bodies): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../src lagoon_bodies_test.cpp ../../src/lagoon_store.cpp -lz -o lagoon_bodies_test
// Run:    ./lagoon_bodies_test     (writes lagoon_bodies_test.lagoon in the current folder)
// The files follow LAGOON's lagoon-format SPEC.md 0.2, section 9; LAGOON reads them
// with lagoon_format.store (Store(...).bodies["nhflow_body"]).
// Architect: Hans Bihs
#include"lagoon_store.h"
#include<cmath>
#include<cstdio>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<sys/stat.h>
#include<vector>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static std::string slurp(const std::string &file)
{
    std::ifstream in(file.c_str(), std::ios::binary);
    std::stringstream s;
    s<<in.rdbuf();
    return s.str();
}

static bool exists(const std::string &file)
{
    struct stat info;
    return stat(file.c_str(), &info)==0;
}

// rotation matrix (row major) of a unit quaternion (w, x, y, z), as REEF3D's E G^T
static void matrix(const double e[4], double R[9])
{
    const double E[3][4] = {{-e[1], e[0], -e[3], e[2]}, {-e[2], e[3], e[0], -e[1]}, {-e[3], -e[2], e[1], e[0]}};
    const double G[3][4] = {{-e[1], e[0], e[3], -e[2]}, {-e[2], -e[3], e[0], e[1]}, {-e[3], e[2], -e[1], e[0]}};
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
        {
            R[3*i+j] = 0.0;
            for(int k=0; k<4; ++k)
                R[3*i+j] += E[i][k]*G[j][k];
        }
}

// a box of 12 triangles, each with its own three points (as REEF3D keeps them)
static std::vector<double> box(double lx, double ly, double lz)
{
    const int faces[6][4] = {{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}};
    const double corner[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    std::vector<double> points;
    for(const int *f : faces)
        for(int tri=0; tri<2; ++tri)
        {
            const int k[3] = {f[0], f[1+tri], f[2+tri]};
            for(int q=0; q<3; ++q)
            {
                points.push_back((corner[k[q]][0]-0.5)*lx);
                points.push_back((corner[k[q]][1]-0.5)*ly);
                points.push_back((corner[k[q]][2]-0.5)*lz);
            }
        }
    return points;
}

int main()
{
    std::cout<<"LAGOON body writer"<<std::endl;
    const std::string path = "./lagoon_bodies_test.lagoon";
    std::remove((path + "/zarr.json").c_str());

    // quaternions: REEF3D's matrix of e and back
    {
        double e[4] = {0.3, -0.5, 0.7, 0.2};
        double norm = std::sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]+e[3]*e[3]);
        for(double &v : e) v /= norm;
        double R[9], q[4];
        matrix(e, R);
        lagoon_bodies::quaternion(R, q);
        double worst = 0.0;
        for(int i=0; i<4; ++i) worst = std::max(worst, std::fabs(q[i]-e[i]));
        check(worst < 1e-14, "quaternion of REEF3D's rotation matrix (w, x, y, z)");
        const double turned[4] = {-e[0], -e[1], -e[2], -e[3]};  // the same rotation
        matrix(turned, R);
        lagoon_bodies::quaternion(R, q);
        check(q[0] >= 0.0 && std::fabs(q[1]-e[1]) < 1e-14, "w >= 0");
    }

    // two bodies, 8 outputs: they heave, drift and turn about a tilted axis
    lagoon_bodies writer(path, "NHFLOW", "nhflow_body", "REEF3D_NHFLOW_6DOF_VTP",
                         "{\"type\": \"run\", \"run\": \"test\"}");
    const std::vector<double> x0[2] = {box(2.0, 1.0, 0.5), box(0.6, 0.6, 0.6)};
    bool all = true;
    for(int t=0; t<8; ++t)
        for(int b=0; b<2; ++b)
        {
            const double angle = 0.1*t*(b+1);
            double e[4] = {std::cos(angle/2), 0.3*std::sin(angle/2), 0.0, std::sqrt(0.91)*std::sin(angle/2)};
            double R[9];
            matrix(e, R);
            const double c[3] = {700.0 + 0.5*t + 10*b, 40.0, 1.5 + 0.2*std::sin(t)};
            const int points = int(x0[b].size()/3);
            std::vector<double> x(x0[b].size());
            for(int i=0; i<points; ++i)
                for(int r=0; r<3; ++r)
                    x[3*i+r] = R[3*r]*x0[b][3*i] + R[3*r+1]*x0[b][3*i+1] + R[3*r+2]*x0[b][3*i+2] + c[r];
            all = writer.output(b, points, x0[b].data(), x.data(), R, c, 0.5*t, t) && all;
        }
    check(all, "8 outputs of 2 rigid bodies in the store");
    const std::string set = path + "/bodies/nhflow_body";
    check(exists(path + "/zarr.json") && exists(path + "/bodies/zarr.json"), "store root and bodies group");
    const std::string attrs = slurp(set + "/zarr.json");
    check(attrs.find("\"bodies\": [0, 1]")!=std::string::npos && attrs.find("\"dataset\": \"nhflow_body\"")!=std::string::npos,
          "set attributes: both bodies, the dataset key");
    check(slurp(set + "/time/zarr.json").find("\"shape\": [8]")!=std::string::npos, "time counts 8 outputs");
    check(slurp(set + "/body_1/rotation/zarr.json").find("\"shape\": [8, 4]")!=std::string::npos, "rotation (8, 4)");
    check(slurp(set + "/body_0/vertices/zarr.json").find("\"shape\": [36, 3]")!=std::string::npos, "vertices once (36, 3)");
    check(exists(set + "/body_0/translation/c/0/0") && exists(set + "/body_0/triangles/c/0/0"), "chunks written");

    // a body that changes shape: refused, and nothing more is written
    {
        std::vector<double> x(x0[0]);
        for(double &v : x) v *= 1.1;
        const double R[9] = {1,0,0, 0,1,0, 0,0,1}, c[3] = {0,0,0};
        check(!writer.output(0, int(x.size()/3), x0[0].data(), x.data(), R, c, 4.0, 8) && !writer.usable(),
              "a body off its rigid motion is refused (its VTP is written instead)");
        check(slurp(set + "/time/zarr.json").find("\"shape\": [8]")!=std::string::npos, "still 8 outputs");
    }

    std::cout<<(nfail ? "FAILED" : "all passed")<<std::endl;
    return nfail ? 1 : 0;
}
