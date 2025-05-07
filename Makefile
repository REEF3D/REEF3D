BUILD        := ./build
BIN          := ./bin
TARGET       := REEF3D
CXX          := mpicxx
GIT_BRANCH   := $(shell git rev-parse --abbrev-ref HEAD)
GIT_VERSION  := "$(shell git describe --abbrev=8 --dirty --always --tags)"
OBJ_DIR      := $(BUILD)
APP_DIR      := $(BIN)
HYPRE_DIR    := /usr/local/hypre
EIGEN_DIR    := ThirdParty/eigen-3.3.8 
CXXFLAGS     := -std=c++11 -DVERSION=\"$(GIT_VERSION)\" -DBRANCH=\"$(GIT_BRANCH)\"
LDFLAGS      := -L ${HYPRE_DIR}/lib/ -lHYPRE
INCLUDE      := -I ${HYPRE_DIR}/include -I ${EIGEN_DIR} -DEIGEN_MPL2_ONLY 
SRC          := $(wildcard src/*.cpp)
OBJECTS      := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

.PHONY: all build clean debug info

all: CXXFLAGS += -O3 -w
all: build

debug: CXXFLAGS += -O0 -g -g3 -Wall
debug: build

build: $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
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
