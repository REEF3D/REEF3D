// Architect: Hans Bihs
// seastate_test: seastate_implicit::run_threads serial (seastate_store::fill shares out its tiles to the threads of
// the solver since Phase 13; the unit test links without the solver)
#include<functional>
class seastate_implicit
{
public:
    static void run_threads(int nt, const std::function<void(int)> &f);
};
void seastate_implicit::run_threads(int nt, const std::function<void(int)> &f)
{
    for(int t=0; t<nt; ++t)
    f(t);
}
