#include "6DOF_forcesext_file.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

void sixdof_forcesext_file::read_format(lexer *p, ghostcell *pgc)
{
    char name[100];
    double val,val0,val1;
    double sign,beta,s;
    int count;

    std::cout << "in read_format\n";

    sprintf(name,"6DOF_forces.dat");

    // open file and determine number of columns
    ifstream file(name, ios_base::in);

    if (!file)
    {
        std::cerr
            << "no '6DOF_forces.dat' file found"
            << " rank=" << p->mpirank
            << std::endl;

        std::abort();
    }

    std::string line;
    std::getline(file, line);

    std::istringstream iss(line);

    colnum = 0;

    while (iss >> val)
        ++colnum;

    std::cerr
        << "rank=" << p->mpirank
        << " colnum=" << colnum
        << std::endl;

    // rewind file
    file.clear();
    file.seekg(0);

    // count rows
    count = 0;

    while (true)
    {
        bool complete_row = true;

        for (qn = 0; qn < colnum; ++qn)
        {
            if (!(file >> val))
            {
                complete_row = false;
                break;
            }
        }

        if (!complete_row)
            break;

        ++count;
    }

    ptnum = count;

    file.close();

    // allocate
    p->Darray(data, ptnum, colnum);

    // re-open file
    file.open(name, ios_base::in);

    if (!file)
    {
        std::cout
            << std::endl
            << "no '6DOF_forces.dat' file found"
            << std::endl
            << std::endl;

        std::abort();
    }

    // read file
    rowcount = 0;

    while (rowcount < ptnum)
    {
        for (qn = 0; qn < colnum; ++qn)
        {
            if (!(file >> data[rowcount][qn]))
            {
                std::abort();
            }
        }

        ++rowcount;
    }

    ts = data[0][0];
    te = data[ptnum-1][0];
}
void sixdof_forcesext_file::forcesext_trans(
    lexer *p,
    ghostcell *pgc
)
{
    Xext = 0.0;
    Yext = 0.0;
    Zext = 0.0;
    std::cout << "simtime: " << p->simtime << "\n";
    if (p->simtime < ts || p->simtime > te){
        std::cout << "returned here\n";
        return;

    }
    int rank = p->mpirank;



    if (!data)
    {
        return;
    }

    if (!data[timecount])
    {
        return;
    }

    while (
        timecount < ptnum - 1 &&
        p->simtime > data[timecount][0]
    )
    {
        ++timecount;

    }


    if (timecount == 0)
    {
        Xext = data[0][1];
        Yext = data[0][2];
        Zext = data[0][3];
        return;
    }

    const int i0 = timecount - 1;
    const int i1 = timecount;

    const double t0 = data[i0][0];
    const double t1 = data[i1][0];

    const double alpha = (p->simtime - t0) / (t1 - t0);
    Yext = data[i0][1] + alpha * ( data[i1][1] - data[i0][1] );

    Xext = data[i0][2] + alpha * ( data[i1][2] - data[i0][2] );

    Zext = data[i0][3] + alpha * ( data[i1][3] - data[i0][3] );
}
void sixdof_forcesext_file::forcesext_rot(    lexer *p,
ghostcell *pgc){
    Kext = 0.0;
    Mext = 0.0;
    Next = 0.0;
    if (colnum !=7)
        return;
    

    std::cout << "simtime: " << p->simtime << "\n";
    if (p->simtime < ts || p->simtime > te){
        std::cout << "returned here\n";
        return;

    }
    int rank = p->mpirank;



    if (!data)
    {
        return;
    }

    if (!data[timecount])
    {
        return;
    }

    while (
        timecount < ptnum - 1 &&
        p->simtime > data[timecount][0]
    )
    {
        ++timecount;

    }


    if (timecount == 0)
    {
        Kext = data[0][4];
        Mext = data[0][5];
        Next = data[0][6];
        return;
    }

    const int i0 = timecount - 1;
    const int i1 = timecount;

    const double t0 = data[i0][0];
    const double t1 = data[i1][0];

    const double alpha = (p->simtime - t0) / (t1 - t0);
    Kext = data[i0][4] + alpha * ( data[i1][4] - data[i0][4] );

    Mext = data[i0][5] + alpha * ( data[i1][5] - data[i0][5] );

    Next = data[i0][6] + alpha * ( data[i1][6] - data[i0][6] );
    std::cout << " Kext " << Kext << " Mext " << Mext<< " Next " << Next << "\n";
}
sixdof_forcesext_file::sixdof_forcesext_file(lexer* p, ghostcell* g){
    ini(p,g);
    read_format(p,g);
}

void sixdof_forcesext_file::ini(lexer* p, ghostcell* g){
    return;
}
