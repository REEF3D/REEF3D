#!/bin/sh
# Builds MoorDyn v2 (MoorDyn-C, the shared library) and the comparison program coral_moordyn_compare
# (REEF3D::CORAL against MoorDyn), from Regression/moordyn:
#   ./build_moordyn.sh [<work folder>]      (default ./moordyn_build; MD_TAG=v2.7.2, CXX=g++, JOBS=4)
# MoorDyn is cloned at the tag into <work folder>/MoorDyn and built without VTK and wrappers; the program
# is linked with the CORAL sources (all folders except reef3d/ and app/) and an rpath to the library.
# Architect: Hans Bihs
set -e
W=${1:-./moordyn_build}
TAG=${MD_TAG:-v2.7.2}
CXX=${CXX:-g++}
JOBS=${JOBS:-4}
mkdir -p "$W"
W=$(cd "$W" && pwd)
if [ ! -d "$W/MoorDyn" ]; then
    git clone --depth 1 --branch "$TAG" https://github.com/FloatingArrayDesign/MoorDyn.git "$W/MoorDyn"
fi
mkdir -p "$W/build"
(cd "$W/build" && cmake ../MoorDyn -DCMAKE_BUILD_TYPE=Release -DFORTRAN_WRAPPER=OFF -DPYTHON_WRAPPER=OFF \
     -DMATLAB_WRAPPER=OFF -DBUILD_TESTING=OFF -DUSE_VTK=OFF > cmake.log && make -j "$JOBS" moordyn > make.log)
FLAGS="-O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src"
OBJ="$W/coral_obj"
mkdir -p "$OBJ"
for f in $(find ../../src/coral -name '*.cpp' -not -path '*/reef3d/*' -not -path '*/app/*'); do
    echo "$f"
done | xargs -P "$JOBS" -I{} sh -c 'o="$0/$(basename {} .cpp).o"; '"$CXX $FLAGS"' -c {} -o "$o" || exit 255' "$OBJ"
$CXX $FLAGS -I"$W/MoorDyn/source" -I"$W/build/source" coral_moordyn_compare.cpp "$OBJ"/*.o \
    -L"$W/build/source" -lmoordyn -Wl,-rpath,"$W/build/source" -o coral_moordyn_compare
echo "built coral_moordyn_compare (MoorDyn $TAG in $W)"
