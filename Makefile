OBJ_DIR      := ./build
APP_DIR      := ./bin
TARGET       := REEF3D
APP		     := $(APP_DIR)/$(TARGET)
CXX          := mpicxx
GIT_BRANCH   := $(shell git rev-parse --abbrev-ref HEAD)
GIT_VERSION  := "$(shell git describe --dirty --always --tags)"
HYPRE_DIR    := /usr/local/hypre
EIGEN_DIR    := ThirdParty/eigen-3.3.8 
CXXFLAGS     := -std=c++11 -DVERSION=\"$(GIT_VERSION)\" -DBRANCH=\"$(GIT_BRANCH)\"
LDFLAGS      := -L ${HYPRE_DIR}/lib/ -lHYPRE
INCLUDE      := -I ${HYPRE_DIR}/include -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY 
SRC          := $(wildcard src/*.cpp)
OBJECTS      := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

.PHONY: all clean debug dev info release

.DEFAULT_GOAL := all

all: CXXFLAGS += -O3 -w
all: CXXFLAGS += -DBUILD=\"all\"
all: $(APP)

release: CXXFLAGS += -O3 -DNDEBUG -DEIGEN_NO_DEBUG -march=native -flto -w
release: CXXFLAGS += -DBUILD=\"release\"
release: LDFLAGS += -flto
release: $(APP)

dev: CXXFLAGS += -O3 -Wall -pedantic -Wpedantic -Wextra -Wshadow -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 #-Wold-style-cast 
dev: CXXFLAGS += -DBUILD=\"dev\"
dev: $(APP)

debug: CXXFLAGS += -O0 -g -g3 -Wall
debug: CXXFLAGS += -DBUILD=\"debug\"
debug: $(APP)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -MMD -c $< -o $@

$(APP): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

-include $(DEPENDENCIES)

clean:
	-@rm -rvf $(APP_DIR) $(OBJ_DIR)

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
