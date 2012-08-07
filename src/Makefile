# $Id: Makefile,v 1.4 2008/04/24 22:15:48 rsz Exp $
# For posix systems.

VTKINC = -I/usr/include/vtk
VTKLIB = -L/usr/lib64 -lvtkCommon -lvtkFiltering -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO  -lvtkRendering
MAYA = /usr/autodesk/maya2012-x64
MAYALIB = -L$(MAYA)/lib -lFoundation -lOpenMaya -lOpenMayaFX

INCLUDES = -I$(MAYA)/include/ $(VTKINC)
CFLAGS = -O3 -DVTK_USE_ANSI_STDLIB -fPIC -pthread -pipe -DLINUX -D_BOOL -DREQUIRE_IOSTREAM -Wno-deprecated -fno-gnu-keywords
LIBS = $(MAYALIB) $(VTKLIB)
CXXFLAGS = $(CFLAGS)

CC = /usr/local/gcc-4.1.2/bin/gcc412
CXX = /usr/local/gcc-4.1.2/bin/g++412

OBJS = importvtk.o pluginMain.o

.SUFFIXES: .cpp .cc .o .so .c

.c.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) $<

.cc.o:
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) $<

.cpp.o:
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) $<

.cc.i:
	$(CXX) -E $(INCLUDES) $(CXXFLAGS) $*.cc > $*.i

.cc.so:
	-rm -f $@
	$(LD) -o $@ $(INCLUDES) $< $(LIBS)

.cpp.so:
	-rm -f $@
	$(LD) -o $@ $(INCLUDES) $< $(LIBS)

.o.so:
	-rm -f $@
	$(LD) -o $@ $< $(LIBS)

all: importvtk.so

importvtk.so: $(OBJS)
	$(CXX) -shared -o $@ $(OBJS) $(LIBS)

clean:
	rm -f *.o *.so
