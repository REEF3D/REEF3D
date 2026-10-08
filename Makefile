HYPRE        ?= 0
OMP          ?= auto

OBJ_DIR      := ./build
APP_DIR      := ./bin
TARGET       := REEF3D
APP          := $(APP_DIR)/$(TARGET)
CXX          := mpicxx
CXXFLAGS     := -std=c++20 -pthread
GIT_BRANCH   := $(shell git rev-parse --abbrev-ref HEAD)
GIT_COMMIT   := $(shell git rev-parse --short=7 HEAD)
GIT_DIRTY    := $(shell git diff --quiet --ignore-submodules HEAD -- || echo -dirty)
GIT_VERSION  := $(GIT_COMMIT)$(GIT_DIRTY)
CXXFLAGS     += -DVERSION=\"$(GIT_VERSION)\" -DBRANCH=\"$(GIT_BRANCH)\"
EIGEN_DIR    := ThirdParty/eigen-5.0.0
INCLUDE      := -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY
LDFLAGS      := -lz -pthread
SRC          := $(filter-out src/hypre_%.cpp,$(wildcard src/*.cpp))

# hypre is opt-in, for benchmarking only: make HYPRE=1 builds the hypre solvers N 10 10-39 (needs hypre in HYPRE_DIR).
# The default build has no hypre dependency; N 10 10-39 then fall back to REEFMG (N 10 1).
ifeq ($(HYPRE),1)
OBJ_DIR      := ./build_hypre
HYPRE_DIR    := /usr/local/hypre
CXXFLAGS     += -DREEF3D_USE_HYPRE
INCLUDE      += -I ${HYPRE_DIR}/include
LDFLAGS      += -L ${HYPRE_DIR}/lib/ -lHYPRE
SRC          := $(wildcard src/*.cpp)
endif

# OpenMP threads for the FEM solid (Z 30): on when the compiler can build OpenMP code (OMP=auto,
# the default; OMP=0 builds without). GCC: -fopenmp; Apple clang: libomp (brew install libomp),
# found with brew --prefix libomp or OMP_PREFIX=<dir>
ifeq ($(shell uname -s)$(shell $(CXX) --version 2>/dev/null | grep -c -i clang),Darwin1)
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

.PHONY: all clean debug dev info release

.DEFAULT_GOAL := release

all: CXXFLAGS += -O3 -w
all: CXXFLAGS += -DBUILD=\"all\"
all: $(APP)

release: CXXFLAGS += -O3 -DNDEBUG -DEIGEN_NO_DEBUG -march=native -flto -w
release: CXXFLAGS += -DBUILD=\"release\"
release: LDFLAGS += -flto=auto
release: $(APP)

dev: CXXFLAGS += -O3 -Wall -pedantic -Wpedantic -Wextra -Wshadow -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 #-Wold-style-cast 
dev: CXXFLAGS += -DBUILD=\"dev\"
dev: $(APP)

debug: CXXFLAGS += -O0 -g -g3 -Wall
debug: CXXFLAGS += -DBUILD=\"debug\"
debug: $(APP)

# Vectorised libm (glibc libmvec) for the irregular-wave component sums (Linux/GCC)
ifeq ($(shell uname -s),Linux)
$(OBJ_DIR)/src/wave_lib_irregular_1st.o: CXXFLAGS += -fopenmp-simd -fno-math-errno -DREEF3D_SIMD_MATH
LDFLAGS += -lmvec
endif

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(INCLUDE) -MMD -MP -c $< -o $@

$(APP): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(OMPLIBS)

-include $(DEPENDENCIES)

clean:
	-@rm -rvf $(APP_DIR) $(OBJ_DIR)

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
