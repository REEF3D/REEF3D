# Build:  make -j<N> [release|all|dev|debug|coral]      (release is the default; N = number of cores, e.g. -j18)
#
# Every configuration has its own object directory (build/<config>), so switching between
# release/all/dev/debug never mixes objects and never forces a full rebuild of the other one.
#
# Build-time switches (all optional):
#   HYPRE=1   build the hypre solvers N 10 10-39 (benchmarking only, needs hypre in HYPRE_DIR);
#             the default build has no hypre dependency; N 10 10-39 then fall back to REEFMG (N 10 1)
#   PCH=0     no precompiled header (src/reef3d_pch.h: standard library, MPI, Eigen); PCH=1 is the default
#             and halves the compile time, the frontend (header parsing) is ~80% of it without
#   LTO=0     release without link-time optimisation (faster link); clang uses ThinLTO (parallel,
#             incremental with a cache in build/release/lto-cache), LTO=full the previous full LTO
#             (single-threaded link); GCC -flto=auto
#   OMP=0     no OpenMP threads for the FEM solid (Z 30); OMP=auto (default) turns them on when available
#   The REEFMG preconditioner on the Apple Silicon GPU (N 16 1, src/reefmg_device_metal.cpp) is part of
#   every macOS build: nothing to switch on, no framework to link (Metal is loaded at run time)
#
# Layout of this file:
#   1. configuration       goal -> CONFIG, build directories, compiler, git version
#   2. sources and flags   sources, include paths, libraries, flags per configuration, OpenMP
#   3. special objects     git version, SIMD wave objects, precompiled header (variables only)
#   4. rules               compile, link, precompiled header, CORAL, clean, info

# ── 1. configuration ─────────────────────────────────────────────────────────
HYPRE        ?= 0
PCH          ?= 1
LTO          ?= 1
OMP          ?= auto

CONFIG       := $(or $(firstword $(filter release all dev debug coral,$(MAKECMDGOALS))),release)
ifeq ($(HYPRE),1)
  BUILD_ROOT   := ./build_hypre
else
  BUILD_ROOT   := ./build
endif
OBJ_DIR      := $(BUILD_ROOT)/$(CONFIG)
APP_DIR      := ./bin
TARGET       := REEF3D
APP          := $(APP_DIR)/$(TARGET)
CXX          := mpicxx
IS_CLANG     := $(shell $(CXX) --version 2>/dev/null | grep -c -i clang)

GIT_BRANCH   := $(shell git rev-parse --abbrev-ref HEAD)
GIT_COMMIT   := $(shell git rev-parse --short=7 HEAD)
GIT_DIRTY    := $(shell git diff --quiet --ignore-submodules HEAD -- || echo -dirty)
GIT_VERSION  := $(GIT_COMMIT)$(GIT_DIRTY)

# ── 2. sources, include paths, libraries ─────────────────────────────────────
CXXFLAGS     := -std=c++20 -pthread
EIGEN_DIR    := ThirdParty/eigen-5.0.0
INCLUDE      := -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY -I src
# dlopen of user PTO controllers (CORAL pto ... user <library>): libdl on Linux (glibc < 2.34), in libSystem on macOS
DL_LIBS      :=
ifeq ($(shell uname -s),Linux)
  DL_LIBS    := -ldl
endif
# zlib: the LAGOON store output (P 18) compresses its chunks with gzip, in a thread
LDFLAGS      := -lz -pthread $(DL_LIBS)
SRC          := $(filter-out src/hypre_%.cpp,$(wildcard src/*.cpp))

ifeq ($(HYPRE),1)
  CXXFLAGS   += -DREEF3D_USE_HYPRE
  HYPRE_DIR  := /usr/local/hypre
  INCLUDE    += -I ${HYPRE_DIR}/include
  LDFLAGS    += -L ${HYPRE_DIR}/lib/ -lHYPRE
  SRC        := $(wildcard src/*.cpp)
endif

# REEF3D::CORAL (src/coral/<folder>/): the solver build takes every folder except app/ (the standalone
# CORAL executable), bin/CORAL every folder except reef3d/; its includes start at src/
# (#include"coral/core/coral_types.h"). Appended after the HYPRE block, which replaces SRC.
CORAL_ALL    := $(shell find src/coral -name '*.cpp' 2>/dev/null | sort)
SRC          += $(filter-out src/coral/app/%,$(CORAL_ALL))
CORAL_SRC    := $(filter-out src/coral/reef3d/%,$(CORAL_ALL))

# ── compiler flags per configuration ─────────────────────────────────────────
# Eigen types returned by value make GCC note the C++14/C++17 ABI change (-Wpsabi)
ifeq ($(IS_CLANG),0)
  CXXFLAGS += -Wno-psabi
endif

ifeq ($(CONFIG),all)
  CONFIG_FLAGS := -O3 -w
else ifeq ($(CONFIG),coral)
  CONFIG_FLAGS := -O3 -DNDEBUG -DEIGEN_NO_DEBUG -w
else ifeq ($(CONFIG),debug)
  CONFIG_FLAGS := -O0 -g -g3 -Wall
else ifeq ($(CONFIG),dev)
  CONFIG_FLAGS := -O3 -Wall -pedantic -Wpedantic -Wextra -Wshadow -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2
else ifeq ($(CONFIG),release)
  CONFIG_FLAGS := -O3 -DNDEBUG -DEIGEN_NO_DEBUG -march=native -w
endif

# Link-time optimisation (release only). GCC: -flto=auto. clang: ThinLTO, whose optimisation runs in
# parallel and, with the cache, only for the changed objects; full LTO (LTO=full) runs single-threaded
# over all ~1500 objects on every link
ifeq ($(CONFIG),release)
  ifneq ($(LTO),0)
    LTO_MODE     := $(if $(filter 0,$(IS_CLANG)),auto,$(if $(filter full,$(LTO)),full,thin))
    CONFIG_FLAGS += -flto=$(LTO_MODE)
    LTO_LDFLAGS  := -flto=$(LTO_MODE)
    ifeq ($(LTO_MODE)$(shell uname -s),thinDarwin)
      LTO_LDFLAGS += -Wl,-cache_path_lto,$(OBJ_DIR)/lto-cache
    endif
  endif
endif

CXXFLAGS += $(CONFIG_FLAGS) -DBUILD=\"$(CONFIG)\"

# OpenMP threads for the FEM solid (Z 30): on when the compiler can build OpenMP code (OMP=auto,
# the default; OMP=0 builds without). GCC: -fopenmp; Apple clang: libomp (brew install libomp),
# found with brew --prefix libomp or OMP_PREFIX=<dir>
ifeq ($(shell uname -s)$(IS_CLANG),Darwin1)
  OMP_PREFIX   ?= $(shell brew --prefix libomp 2>/dev/null)
  OMP_CFLAGS   := -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include
  OMP_LIBS     := -L$(OMP_PREFIX)/lib -lomp
else
  OMP_CFLAGS   := -fopenmp
  OMP_LIBS     := -fopenmp
endif
ifeq ($(OMP),auto)
  USE_OMP      := $(shell printf '\043include <omp.h>\nint main(){return omp_get_max_threads()>0 ? 0 : 1;}\n' | \
                $(CXX) $(OMP_CFLAGS) -x c++ - -o /dev/null $(OMP_LIBS) >/dev/null 2>&1 && echo 1 || echo 0)
else
  USE_OMP      := $(OMP)
endif
ifeq ($(USE_OMP),1)
  OMPFLAGS     := $(OMP_CFLAGS)
  OMPLIBS      := $(OMP_LIBS)
endif

OBJECTS      := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

# ── 3. special objects ───────────────────────────────────────────────────────
# git version: only driver.cpp and driver_log.cpp print it. A stamp file that changes with the
# branch/commit/dirty state recompiles just these two files, so a new commit neither rebuilds
# everything nor leaves a stale version in the executable
VERSION_STAMP := $(OBJ_DIR)/git_version.txt
VERSION_OBJS  := $(OBJ_DIR)/src/driver.o $(OBJ_DIR)/src/driver_log.o
ifeq ($(filter clean,$(MAKECMDGOALS)),)
  _VERSION_LINE := $(GIT_BRANCH) $(GIT_VERSION)
  _VERSION_OLD  := $(shell cat $(VERSION_STAMP) 2>/dev/null)
  ifneq ($(_VERSION_LINE),$(_VERSION_OLD))
    _DUMMY        := $(shell mkdir -p $(OBJ_DIR) && echo '$(_VERSION_LINE)' > $(VERSION_STAMP))
  endif
endif
$(VERSION_OBJS): CXXFLAGS += -DVERSION=\"$(GIT_VERSION)\" -DBRANCH=\"$(GIT_BRANCH)\"

# Irregular-wave component sums: vectorised libm (glibc libmvec) on Linux/GCC, SIMD reductions over
# the second-order terms (#pragma omp simd, GCC and Clang)
WAVE_OBJ       := $(OBJ_DIR)/src/wave_lib_irregular
WAVE_SIMD_OBJS := $(WAVE_OBJ)_1st.o $(WAVE_OBJ)_2nd_a.o $(WAVE_OBJ)_2nd_b.o
ifeq ($(shell uname -s),Linux)
  $(WAVE_OBJ)_1st.o: CXXFLAGS += -fopenmp-simd -fno-math-errno -DREEF3D_SIMD_MATH
  LDFLAGS += -lmvec
else
  $(WAVE_OBJ)_1st.o: CXXFLAGS += -fopenmp-simd
endif
$(WAVE_OBJ)_2nd_a.o $(WAVE_OBJ)_2nd_b.o: CXXFLAGS += -fopenmp-simd

# Precompiled header: src/reef3d_pch.h (standard library, MPI, Eigen) is compiled once per configuration
# and included into every source through PCH_FLAGS. GCC needs the copy in the build directory (it reads
# <header>.gch next to the header); clang reads the .pch file directly. PCH=0 disables it.
USE_PCH      := $(if $(filter coral,$(CONFIG)),0,$(PCH))

# The VERSION/BRANCH objects and the SIMD wave objects have their own flags: they are built without
# the PCH and do not depend on it, so the PCH never inherits their target-specific flags
NO_PCH_OBJS  := $(VERSION_OBJS) $(WAVE_SIMD_OBJS)
$(NO_PCH_OBJS): PCH_FLAGS :=

ifeq ($(USE_PCH),1)
  PCH_DIR      := $(OBJ_DIR)/pch
  PCH_COPY     := $(PCH_DIR)/reef3d_pch.h
  ifeq ($(IS_CLANG),0)
    PCH_OUT    := $(PCH_COPY).gch
    PCH_FLAGS  := -include $(PCH_COPY)
  else
    PCH_OUT    := $(PCH_COPY).pch
    PCH_FLAGS  := -include-pch $(PCH_OUT)
  endif
endif

# ── 4. rules ─────────────────────────────────────────────────────────────────
.PHONY: all clean coral debug dev info release

.DEFAULT_GOAL := release

release all dev debug: $(OBJ_DIR)/$(TARGET)
	@mkdir -p $(APP_DIR)
	@cp -p $(OBJ_DIR)/$(TARGET) $(APP)
	@echo "[*] $(CONFIG): $(APP)"

$(OBJ_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(OBJ_DIR)/lto-cache
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LTO_LDFLAGS) $(LDFLAGS) $(OMPLIBS)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDE) $(PCH_FLAGS) -MMD -MP -c $< -o $@

$(VERSION_OBJS): $(VERSION_STAMP)

ifeq ($(USE_PCH),1)
$(PCH_COPY): src/reef3d_pch.h
	@mkdir -p $(@D)
	@cp $< $@

$(PCH_OUT): $(PCH_COPY)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDE) -x c++-header -MMD -MP -MT $@ -MF $@.d -c $< -o $@

$(filter-out $(NO_PCH_OBJS),$(OBJECTS)): $(PCH_OUT)

-include $(PCH_OUT).d
endif

# REEF3D::CORAL standalone executable (bin/CORAL): src/coral without reef3d/, with app/; no MPI needed
CORAL_APP     := $(APP_DIR)/CORAL
CORAL_OBJECTS := $(CORAL_SRC:%.cpp=$(OBJ_DIR)/%.o)

coral: $(CORAL_APP)

$(CORAL_APP): $(CORAL_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OMPLIBS) $(DL_LIBS)

ifeq ($(CONFIG),coral)
  -include $(CORAL_OBJECTS:.o=.d)
else
  -include $(DEPENDENCIES)
endif

clean:
	-@rm -rf $(APP_DIR) ./build ./build_hypre

info:
	@echo "[*] Configuration:   $(CONFIG)      "
	@echo "[*] Compiler:        $(CXX) (clang: $(IS_CLANG))"
	@echo "[*] Application dir: $(APP_DIR)     "
	@echo "[*] Object dir:      $(OBJ_DIR)     "
	@echo "[*] PCH / LTO / OMP: $(USE_PCH) / $(LTO) / $(USE_OMP)"
	@echo "[*] Sources:         $(SRC)         "
