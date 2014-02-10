include ../Makefile.in

BINDIR = ..

INCLUDES = -I../Lib 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmlrmcl -lm
LD = $(CC) $(LDOPTIONS) -L. -L..
MLRMCLLIB = ../libmlrmcl.a

MLRMCLOBJ = mlrmcl.o 
RMCLOBJ = rmcl.o 
SRMCLOBJ = srmcl.o
NCUTOBJ = ncut.o 
MERGEOBJ = mergeClusters.o

.c.o:
	
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/rmcl $(BINDIR)/srmcl  
#$(BINDIR)/mlrmcl $(BINDIR)/mergeClusters  $(BINDIR)/ncut
	

#$(BINDIR)/mlrmcl: $(MLRMCLOBJ) $(MLRMCLLIB)
#	$(LD) -o $@ $(MLRMCLOBJ) $(LIBS)
#	chmod 755 $@

$(BINDIR)/rmcl: $(RMCLOBJ) $(MLRMCLLIB)
	$(LD) -o $@ $(RMCLOBJ) $(LIBS)
	chmod 755 $@

#$(BINDIR)/softrmcl: $(SOFTRMCLOBJ) $(MLRMCLLIB)
#	$(LD) -o $@ $(SOFTRMCLOBJ) $(LIBS)
#	chmod 755 $@

$(BINDIR)/srmcl: $(SRMCLOBJ) $(MLRMCLLIB)
	$(LD) -o $@ $(SRMCLOBJ) $(LIBS)
	chmod 755 $@

#$(BINDIR)/mergeClusters: $(MERGEOBJ) $(MLRMCLLIB)
#	$(LD) -o $@ $(MERGEOBJ) $(LIBS)
#	chmod 755 $@

#$(BINDIR)/ncut: $(NCUTOBJ) $(MLRMCLLIB)
#	$(LD) -o $@ $(NCUTOBJ) $(LIBS)
#	chmod 755 $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
        rm -f $(BINDIR)/mlrmcl ;\
        rm -f $(BINDIR)/ncut ;\
		rm -f $(BINDIR)/mergeClusters ;\
