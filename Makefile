.PHONY: all lib programs doxy clean delobj dellib delbin deldoxy

OBJS_LIB = $(shell ls src/*.cpp | sed 's/\.cpp/.o/')
PROGS    = $(shell ls src/*.C | sed 's/\.C//' | sed 's/src\///')
# export $DK2NU to dk2nu installation directory
DK2NU_INC = ${DK2NU}/include
DK2NU_LIB = ${DK2NU}/lib
BOOST_INC = ${BOOST_ROOT}/include
INCLUDES  = -I./include -I$(shell root-config --incdir) -I${BOOST_INC} -I${DK2NU_INC}
ROOT_LIB  = $(shell root-config --libdir)
ROOT_LIBS  =$(shell root-config --libs) -lEG

# for dynamic linking
# LDFLAGS = '-Wl,-rpath,$$ORIGIN/../lib' -L./lib -L${DK2NU_LIB}
# LDLIBS  = -lppfx -ldk2nuTree

# for static linking
LDFLAGS = -L./lib -L${DK2NU_LIB}
LDLIBS  = -l:libppfx.a -l:libdk2nuTree.a

CC = g++
CFLAGS = -fPIC -DLINUX -O3 -g -pthread -std=c++17 -m64 -Wall -Wextra -pedantic
# AW (2025-02-13): flag to toggle ICARUS-specific code
ICARUSFLAG = -DUH_ICARUS
# AW (2025-02-13): This flag was just a test to mimic what would happen if we had data in the x_F < 0 region. Don't use it.
# ICARUSFLAG += -DNUA_XF_MIRRORING
CFLAGS += $(ICARUSFLAG)

all: lib programs doxy

lib: dirs libppfx.so libppfx.a

libppfx.so: $(OBJS_LIB)
	$(CC) -shared -o lib/$@ $^ -Wl,-rpath,${DK2NU_LIB},-rpath,${ROOT_LIB} -L${DK2NU_LIB} -ldk2nuTree

libppfx.a: $(OBJS_LIB)
	ar rcs lib/$@ $^

programs: dirs $(PROGS)
	echo making $(PROGS)

$(PROGS): % : src/%.o $(OBJS_LIB) lib
	$(CC) -Wall -o bin/$@ $< $(PPFX_OBJS) $(ROOT_LIBS) $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

%.o: %.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

%.o: %.C
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

%.o: %.cxx
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

dirs:
	if [ ! -d bin ]; then mkdir -p bin; fi
	if [ ! -d lib ]; then mkdir -p lib; fi

doxy:
	doxygen doxygen/config_doxygen

clean:  deldoxy delobj dellib delbin

delobj:
	-rm src/*.o

dellib:
	if [ -d lib ]; then rm -rf lib; fi

delbin:
	if [ -d bin ]; then rm -rf bin; fi

deldoxy:
	if [ -d html ]; then rm -rf html; fi
	if [ -d latex ]; then rm -rf latex; fi
