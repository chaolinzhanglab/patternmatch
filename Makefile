CC=			gcc
CXX=		g++
#CFLAGS=		-g -Wall -pg #-O2 -m64 -pg
CFLAGS=		-g -Wall -O2 -m64 #-pg
CXXFLAGS=	$(CFLAGS)

OBJS=	

PROG=		PatternMatch

INCLUDES= 	-I./gzstream

LIBS=		-lpopt -lm -lz
LIBS2=		./gzstream/libgzstream.a

SUBDIRS=	. gzstream

.SUFFIXES:.c .o .cpp

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)"  \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

PatternMatch:lib-recur $(OBJS) main.o
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS) $(LIBS2)

PatternMatch.o:FastaFile.hpp

cleanlocal:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

clean:cleanlocal-recur
