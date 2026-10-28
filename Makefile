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
HYPRE        ?= 0
PCH          ?= 1
LTO          ?= 1

# ── configuration from the goal ───────────────────────────────────────────────
CONFIG       := $(firstword $(filter release all dev debug coral,$(MAKECMDGOALS)))
ifeq ($(CONFIG),)
CONFIG       := release
endif
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
HYPRE_DIR    := /usr/local/hypre
EIGEN_DIR    := ThirdParty/eigen-5.0.0
CXXFLAGS     := -std=c++20
ifeq ($(HYPRE),1)
CXXFLAGS     += -DREEF3D_USE_HYPRE
LDFLAGS      := -L ${HYPRE_DIR}/lib/ -lHYPRE
INCLUDE      := -I ${HYPRE_DIR}/include -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY 
SRC          := $(wildcard src/*.cpp)
else
LDFLAGS      :=
INCLUDE      := -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY 
SRC          := $(filter-out src/hypre_%.cpp,$(wildcard src/*.cpp))
endif
# REEF3D::CORAL (src/coral/<folder>/): every folder except app/ (the standalone CORAL executable);
# its includes start at src/ (#include"coral/core/coral_types.h")
SRC          += $(shell find src/coral -name '*.cpp' -not -path 'src/coral/app/*' 2>/dev/null | sort)
INCLUDE      += -I src
# zlib: the LAGOON store output (P 18) compresses its chunks with gzip, in a thread
LDFLAGS      += -lz -pthread

# ── compiler flags per configuration ─────────────────────────────────────────
ifeq ($(CONFIG),release)
CXXFLAGS     += -O3 -DNDEBUG -DEIGEN_NO_DEBUG -march=native -w -DBUILD=\"release\"
ifneq ($(LTO),0)
ifeq ($(IS_CLANG),0)
CXXFLAGS     += -flto=auto
LTO_LDFLAGS  := -flto=auto
else ifeq ($(LTO),full)
CXXFLAGS     += -flto=full
LTO_LDFLAGS  := -flto=full
else
# ThinLTO: the link-time optimisation runs in parallel and, with the cache, only for the changed
# objects; full LTO (-flto with clang) runs single-threaded over all ~1500 objects on every link
CXXFLAGS     += -flto=thin
LTO_LDFLAGS  := -flto=thin
ifeq ($(shell uname -s),Darwin)
LTO_LDFLAGS  += -Wl,-cache_path_lto,$(OBJ_DIR)/lto-cache
endif
endif
endif
else ifeq ($(CONFIG),all)
CXXFLAGS     += -O3 -w -DBUILD=\"all\"
else ifeq ($(CONFIG),dev)
CXXFLAGS     += -O3 -Wall -pedantic -Wpedantic -Wextra -Wshadow -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 #-Wold-style-cast 
CXXFLAGS     += -DBUILD=\"dev\"
else ifeq ($(CONFIG),debug)
CXXFLAGS     += -O0 -g -g3 -Wall -DBUILD=\"debug\"
else ifeq ($(CONFIG),coral)
CXXFLAGS     += -O3 -DNDEBUG -DEIGEN_NO_DEBUG -w -DBUILD=\"coral\"
endif

# OpenMP threads for the FEM solid (Z 30): on when the compiler can build OpenMP code (OMP=auto,
# the default; OMP=0 builds without). GCC: -fopenmp; Apple clang: libomp (brew install libomp),
# found with brew --prefix libomp or OMP_PREFIX=<dir>
OMP          ?= auto
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

# ── precompiled header ───────────────────────────────────────────────────────
# src/reef3d_pch.h is built once per configuration (with the same flags as the sources) and
# force-included into every source; clang picks <header>.pch, GCC <header>.gch next to the copy
USE_PCH      := $(PCH)
ifeq ($(CONFIG),coral)
USE_PCH      := 0
endif
ifeq ($(USE_PCH),1)
PCH_DIR      := $(OBJ_DIR)/pch
ifeq ($(IS_CLANG),0)
PCH_OUT      := $(PCH_DIR)/reef3d_pch.h.gch
else
PCH_OUT      := $(PCH_DIR)/reef3d_pch.h.pch
endif
PCH_FLAGS    := -include $(PCH_DIR)/reef3d_pch.h
# snapshot of the configuration flags: target-specific additions (driver*.o, wave_lib_irregular*.o)
# must not leak into the PCH when make happens to build it as their prerequisite
PCH_CXXFLAGS := $(CXXFLAGS) $(OMPFLAGS) $(INCLUDE)
endif

# ── git version: only driver.cpp and driver_log.cpp print it ─────────────────
# A stamp file that changes with the branch/commit/dirty state recompiles just these two files,
# so a new commit neither rebuilds everything nor leaves a stale version in the executable
VERSION_STAMP := $(OBJ_DIR)/git_version.txt
VERSION_OBJS  := $(OBJ_DIR)/src/driver.o $(OBJ_DIR)/src/driver_log.o
ifeq ($(filter clean,$(MAKECMDGOALS)),)
_VERSION_LINE := $(GIT_BRANCH) $(GIT_VERSION)
_VERSION_OLD  := $(shell cat $(VERSION_STAMP) 2>/dev/null)
ifneq ($(_VERSION_LINE),$(_VERSION_OLD))
_DUMMY        := $(shell mkdir -p $(OBJ_DIR) && echo '$(_VERSION_LINE)' > $(VERSION_STAMP))
endif
endif
$(VERSION_OBJS): $(VERSION_STAMP)
$(VERSION_OBJS): CXXFLAGS += -DVERSION=\"$(GIT_VERSION)\" -DBRANCH=\"$(GIT_BRANCH)\"

.PHONY: all clean debug dev info release coral

.DEFAULT_GOAL := release

release all dev debug: $(OBJ_DIR)/$(TARGET)
	@mkdir -p $(APP_DIR)
	@cp -p $(OBJ_DIR)/$(TARGET) $(APP)
	@echo "[*] $(CONFIG): $(APP)"

# Vectorised libm (glibc libmvec) for the irregular-wave component sums (Linux/GCC)
# These files have their own language flags and therefore do not use the precompiled header.
WAVE_SIMD_OBJS := $(OBJ_DIR)/src/wave_lib_irregular_1st.o $(OBJ_DIR)/src/wave_lib_irregular_2nd_a.o $(OBJ_DIR)/src/wave_lib_irregular_2nd_b.o
$(WAVE_SIMD_OBJS): PCH_FLAGS :=
ifeq ($(shell uname -s),Linux)
$(OBJ_DIR)/src/wave_lib_irregular_1st.o: CXXFLAGS += -fopenmp-simd -fno-math-errno -DREEF3D_SIMD_MATH
LDFLAGS += -lmvec
else
$(OBJ_DIR)/src/wave_lib_irregular_1st.o: CXXFLAGS += -fopenmp-simd
endif

# SIMD reductions over the second-order irregular-wave terms (#pragma omp simd, GCC and Clang)
$(OBJ_DIR)/src/wave_lib_irregular_2nd_a.o $(OBJ_DIR)/src/wave_lib_irregular_2nd_b.o: CXXFLAGS += -fopenmp-simd

ifeq ($(USE_PCH),1)
$(PCH_OUT): src/reef3d_pch.h
	@mkdir -p $(@D)
	@cp src/reef3d_pch.h $(PCH_DIR)/reef3d_pch.h
	$(CXX) $(PCH_CXXFLAGS) -x c++-header -MMD -MP -MT $@ -MF $@.d -c $(PCH_DIR)/reef3d_pch.h -o $@

$(OBJECTS): $(PCH_OUT)

-include $(PCH_OUT).d
endif

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDE) $(PCH_FLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(OBJ_DIR)/lto-cache
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LTO_LDFLAGS) $(LDFLAGS) $(OMPLIBS)

# REEF3D::CORAL standalone executable (bin/CORAL): src/coral without reef3d/, with app/; no MPI needed
CORAL_APP     := $(APP_DIR)/CORAL
CORAL_SRC     := $(shell find src/coral -name '*.cpp' -not -path 'src/coral/reef3d/*' 2>/dev/null | sort)
CORAL_OBJECTS := $(CORAL_SRC:%.cpp=$(OBJ_DIR)/%.o)

coral: $(CORAL_APP)

$(CORAL_APP): $(CORAL_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OMPLIBS)

ifeq ($(CONFIG),coral)
-include $(CORAL_OBJECTS:.o=.d)
else
-include $(DEPENDENCIES)
endif

clean:
	-@rm -rf $(APP_DIR) ./build ./build_hypre

info:
	@echo "[*] Configuration:   ${CONFIG}      "
	@echo "[*] Compiler:        $(CXX) (clang: $(IS_CLANG))"
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] PCH / LTO / OMP: $(USE_PCH) / $(LTO) / $(USE_OMP)"
	@echo "[*] Sources:         ${SRC}         "
