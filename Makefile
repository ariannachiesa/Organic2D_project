BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

all: $(BUILD_DIR)/test1_mesh $(BUILD_DIR)/test2_mesh $(BUILD_DIR)/test3_mesh $(BUILD_DIR)/test_Laplace \
$(BUILD_DIR)/test_LinPoisson $(BUILD_DIR)/test1_NLPoisson $(BUILD_DIR)/test2_NLPoisson $(BUILD_DIR)/test_CVcurve

# Directory where Bim++ is installed.
#BIMPP_PREFIX = $(shell pwd)/../bimpp
#BIMPP_PREFIX = ./bimpp_install
BIMPP_PREFIX = /vagrant/bimpp_install

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

mkOctaveInc = ${mkOctavePrefix}/include/octave-4.2.0/
mkOctaveLib = ${mkOctavePrefix}/lib/octave/4.2.0/
mkP4estInc = /vagrant/p4est/include
mkP4estLib = /vagrant/p4est/lib

CPPFLAGS = -I./src/ -I$(BIMPP_PREFIX)/include \
-I$(mkP4estInc) -I$(mkOctaveInc) -I$(mkOctaveInc)/octave -I${mkMumpsInc}

CXXFLAGS = -Wall -std=c++14
CXX = mpic++

LIBS = $(BIMPP_PREFIX)/lib ${mkOpenblasLib} \
$(mkP4estLib) $(mkOctaveLib) ${mkOpenblasLib} \
${mkMumpsLib} ${mkScotchLib}

comma := ,

LDFLAGS = $(addprefix -Wl$(comma)-rpath$(comma), $(LIBS)) \
$(addprefix -L, $(LIBS))

LDLIBS = -lbim -lbimio -lbimlis \
-lp4est -lsc \
-lbimmumps -lbimp4est \
-lopenblas -loctave -loctinterp \
-lopenblas -lptscotcherr

$(BUILD_DIR)/test1_mesh: $(OBJS) $(BUILD_DIR)/test/test1_mesh.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)
	
$(BUILD_DIR)/test2_mesh: $(OBJS) $(BUILD_DIR)/test/test2_mesh.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)
	
$(BUILD_DIR)/test3_mesh: $(OBJS) $(BUILD_DIR)/test/test3_mesh.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/test_Laplace: $(OBJS) $(BUILD_DIR)/test/test_Laplace.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)
	
$(BUILD_DIR)/test_LinPoisson: $(OBJS) $(BUILD_DIR)/test/test_LinPoisson.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/test1_NLPoisson: $(OBJS) $(BUILD_DIR)/test/test1_NLPoisson.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)
	
$(BUILD_DIR)/test2_NLPoisson: $(OBJS) $(BUILD_DIR)/test/test2_NLPoisson.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)
	
$(BUILD_DIR)/test_CVcurve: $(OBJS) $(BUILD_DIR)/test/test_CVcurve.cpp.o
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)	
	
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

MKDIR_P ?= mkdir -p
