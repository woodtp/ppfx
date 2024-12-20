OBJS_LIB = $(shell ls src/*.cpp | sed 's/\.cpp/.o/')
PROGS = $(shell ls src/*.C | sed 's/\.C//' | sed 's/src\///')
INCLUDES = -I./include -I$(shell root-config --incdir) -I${BOOST_INC} -I${DK2NU}/include
DEPLIBS=$(shell root-config --libs) -lEG
LDFLAGS= '-Wl,-rpath,$$ORIGIN/../lib' -L./lib
LDLIBS = -lppfx -L${DK2NU}/lib -ldk2nuTree

CC = g++
CFLAGS = -fPIC -DLINUX -O3 -g $(shell root-config --cflags) -Wall -Wextra -pedantic
ICARUSFLAG = -DUH_ICARUS
# ICARUSFLAG += -DNUA_XF_MIRRORING
CFLAGS += $(ICARUSFLAG)

all:    lib programs doxy

lib: libppfx.so

libppfx.so: $(OBJS_LIB)
	if [ ! -d lib ]; then mkdir -p lib; fi

	$(CC) -shared -o lib/$@ $^ -L${DK2NU}/lib -ldk2nuTree

programs: $(PROGS)
	echo making $(PROGS)

$(PROGS): % : src/%.o $(OBJS_LIB)  libppfx.so
	if [ ! -d bin ]; then mkdir -p bin; fi

	$(CC) -Wall -o bin/$@ $< $(PPFX_OBJS) $(DEPLIBS) $(LDFLAGS) $(LDLIBS)
%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

%.o: %.cc
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

%.o: %.C
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

%.o: %.cxx
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

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
