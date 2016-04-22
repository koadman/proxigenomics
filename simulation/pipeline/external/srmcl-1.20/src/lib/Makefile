include ../../Makefile.in

INCLUDES = -I. 

CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES) -g
LD = $(CC) -L. 


OBJS = coarsen.o match.o ccgraph.o memory.o \
       util.o metrics.o merge.o\
       graph.o io.o \
       subgraphs.o \
       myqsort.o \
       mlmcl.o mclbase.o mclutils.o rmcl.o \
       srmcl.o 

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp

../../libmlrmcl.a: $(OBJS) *.h 
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f ../../libmlrmcl.a
