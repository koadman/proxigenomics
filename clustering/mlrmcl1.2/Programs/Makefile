include ../Makefile.in

BINDIR = ..

INCLUDES = -I../Lib 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmlrmcl -lm
LD = $(CC) $(LDOPTIONS) -L. -L..
MLRMCLLIB = ../libmlrmcl.a

MLRMCLOBJ = mlrmcl.o 
NCUTOBJ = ncut.o 
MERGEOBJ = mergeClusters.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/mlrmcl $(BINDIR)/ncut $(BINDIR)/mergeClusters

$(BINDIR)/mlrmcl: $(MLRMCLOBJ) $(MLRMCLLIB)
	$(LD) -o $@ $(MLRMCLOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/mergeClusters: $(MERGEOBJ) $(MLRMCLLIB)
	$(LD) -o $@ $(MERGEOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/ncut: $(NCUTOBJ) $(MLRMCLLIB)
	$(LD) -o $@ $(NCUTOBJ) $(LIBS)
	chmod 755 $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
        rm -f $(BINDIR)/mlrmcl ;\
        rm -f $(BINDIR)/ncut ;\
		rm -f $(BINDIR)/mergeClusters ;\
