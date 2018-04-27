BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

# Directory where Bim++ is installed.
#BIMPP_PREFIX = ../bimpp
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

CXXFLAGS = -Wall -std=c++11
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

TARGET_EXEC ?= mis_main2D

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

MKDIR_P ?= mkdir -p
