#!/bin/sh
# Builds and runs the REEF3D::CORAL unit tests (no MPI, no REEF3D binary), from Regression/unit:
#   ./coral_tests.sh            all tests (CXX=g++, JOBS=4 by default)
# The CORAL sources are compiled once into a temporary folder, each test links against them.
# Architect: Hans Bihs
CXX=${CXX:-g++}
JOBS=${JOBS:-4}
FLAGS="-O2 -DEIGEN_NO_DEBUG -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src"
# the user PTO controllers are loaded with dlopen (libdl on Linux; part of the C library on macOS)
LIBS=""
[ "$(uname -s)" = "Linux" ] && LIBS="-ldl"
CC=${CC:-cc}
OBJ=$(mktemp -d)
trap 'rm -rf "$OBJ"' EXIT
fail=0
# CORAL objects (all folders except reef3d/ and app/), in parallel
for f in $(find ../../src/coral -name '*.cpp' -not -path '*/reef3d/*' -not -path '*/app/*'); do
    echo "$f"
done | xargs -P "$JOBS" -I{} sh -c 'o="$0/$(basename {} .cpp).o"; '"$CXX $FLAGS"' -c {} -o "$o" || exit 255' "$OBJ" || fail=1
# the example user controller of coral_wec_test
$CC -shared -fPIC -O2 coral_pto_user_example.c -o "$OBJ/coral_pto_user_example.so" || fail=1
for t in coral_test coral_hydro_test coral_pto_test coral_mooring_test coral_io_test coral_beam_test coral_linear_test \
         coral_superelement_test coral_soil_test coral_offshore_test coral_fluid_test coral_wind_test coral_aero_test coral_rotor_test coral_control_test coral_floater_test coral_unsteady_aero_test \
         coral_moordyn_test coral_wec_test coral_slender_test coral_slam_test; do
    if ! $CXX $FLAGS $t.cpp "$OBJ"/*.o $LIBS -o "$OBJ/$t"; then echo "$t: build failed"; fail=1; continue; fi
    if [ "$t" = coral_wec_test ]; then res=$("$OBJ/$t" "$OBJ/coral_pto_user_example.so" | tail -1); else res=$("$OBJ/$t" | tail -1); fi
    echo "$t: $res"
    [ "$res" = "all passed" ] || fail=1
done
python3 coral_layering_check.py | tail -1
exit $fail
