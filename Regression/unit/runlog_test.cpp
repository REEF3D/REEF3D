// Standalone verification of the run log (runlog) for outputs that go into the LAGOON
// store only (P 18 1, no VTU/VTP files): no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src runlog_test.cpp -o runlog_test
// Run:    ./runlog_test     (writes runlog_test_case/REEF3D_Case/REEF3D_CFD_run.jsonl)
// The lines follow LAGOON's lagoon-format RUN.md: a "run" line, a "stream" line per
// stream (format "zarr" with the store folder and the group in it), an "output" line
// per output with the time and iteration of the output, an "end" line.
// Architect: Hans Bihs

// runlog.cpp needs only these members of the lexer: a small stand-in, so the test
// builds without the solver (the include guard keeps the real lexer.h out)
#define LEXER_H_
class runlog;
class lexer
{
public:
    int mpirank = 0, M10 = 4, A10 = 6, count = 0;
    double simtime = 0.0, phimean = 0.5;
    double global_xmin = 0.0, global_ymin = 0.0, global_zmin = 0.0;
    double global_xmax = 10.0, global_ymax = 1.0, global_zmax = 1.0;
    runlog *plog = nullptr;
};
#include"../../src/runlog.cpp"

#include<fstream>
#include<iostream>
#include<string>
#include<sys/stat.h>
#include<unistd.h>
#include<vector>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static std::vector<std::string> lines(const std::string &file)
{
    std::ifstream in(file.c_str());
    std::vector<std::string> all;
    for(std::string line; std::getline(in, line);)
        all.push_back(line);
    return all;
}

static bool has(const std::string &line, const std::string &part)
{
    return line.find(part)!=std::string::npos;
}

int main()
{
    mkdir("runlog_test_case", 0777);
    if(chdir("runlog_test_case")!=0)
        return 1;
    std::remove("REEF3D_Case/REEF3D_CFD_run.jsonl");
    {
        std::ofstream ctrl("ctrl.txt");
        ctrl<<"A 10 6\nP 18 1\n";
    }

    lexer p;
    runlog log(&p, "test");
    const std::string file = "REEF3D_Case/REEF3D_CFD_run.jsonl";
    check(lines(file).empty(), "nothing is written before the first output (the run line is lazy)");

    // the volume store counts output 0 only at the next output: by then the solver is
    // at a later time and iteration, the log takes the ones of the output
    p.simtime = 0.25;
    p.count = 40;
    log.stored(&p, 0, 0.0, 0, "lagoon_volume", "volume", "./REEF3D_CFD.lagoon", "volume", p.M10);
    p.simtime = 0.5;
    p.count = 80;
    log.stored(&p, 1, 0.25, 40, "lagoon_volume", "volume", "./REEF3D_CFD.lagoon", "volume", p.M10);
    log.stored(&p, 3, 0.5, 80, "lagoon_bed", "bed", "./REEF3D_CFD.lagoon", "bed", p.M10);
    log.end(&p, "finished");

    const std::vector<std::string> all = lines(file);
    check(all.size()==7, "run, stream, output, output, stream, output, end: " + std::to_string(all.size()) + " lines");
    if(all.size()==7)
    {
        check(has(all[0], "\"type\":\"run\"") && has(all[0], "\"solver\":\"CFD\"") && has(all[0], "\"ranks\":4"),
              "the run line comes with the first stored output");
        check(has(all[1], "\"type\":\"stream\"") && has(all[1], "\"name\":\"lagoon_volume\""), "the volume stream");
        check(has(all[1], "\"folder\":\"REEF3D_CFD.lagoon\""), "the store folder, without ./");
        check(has(all[1], "\"group\":\"volume\"") && has(all[1], "\"format\":\"zarr\""), "group and format zarr");
        check(has(all[1], "\"pieces\":4") && !has(all[1], "first_step") && !has(all[1], "\"files\""),
              "pieces, no first_step for step 0, no file pattern");
        check(has(all[2], "\"step\":0") && has(all[2], "\"time\":0,") && has(all[2], "\"iteration\":0,"),
              "output 0 with its own time and iteration, not the solver's current ones");
        check(has(all[3], "\"step\":1") && has(all[3], "\"time\":0.25,") && has(all[3], "\"iteration\":40,")
              && has(all[3], "\"streams\":[\"lagoon_volume\"]"), "output 1");
        check(has(all[4], "\"name\":\"lagoon_bed\"") && has(all[4], "\"first_step\":3") && has(all[4], "\"group\":\"bed\""),
              "a second stream, registered once, with its first step");
        check(has(all[5], "\"step\":3") && has(all[5], "\"streams\":[\"lagoon_bed\"]"), "the bed output");
        check(has(all[6], "\"type\":\"end\"") && has(all[6], "\"status\":\"finished\""), "the end line");
    }

    std::cout<<(nfail==0 ? "ALL PASSED" : std::to_string(nfail) + " FAILED")<<std::endl;
    return nfail==0 ? 0 : 1;
}
