override CFLAGS += -I../include -I/usr/include/openblas
#LDLIBS += -lgsl -lopenblas -lrt -lm
override LDLIBS += -lgsl -lopenblas -lm -lpthread
# Some Linux distributions require separate lapacke library
# override LDLIBS += -llapacke
# Mac additions using homebrew installations
#override CFLAGS += -I/usr/local/opt/openblas/include -I/usr/local/opt/gsl/include
#override LDFLAGS += -L/usr/local/opt/openblas/lib -L/usr/local/opt/gsl/lib
# Harvard Medical School O2 cluster additions
#override CFLAGS += -I/n/app/openblas/0.2.19/include -I/n/app/gsl/2.3/include
#override LDFLAGS += -L/n/app/openblas/0.2.19/lib -L/n/app/gsl/2.3/lib/

ifeq ($(OPTIMIZE), 1)
	CFLAGS += -O2
endif

ifeq ($(DEBUG), 1)
	CFLAGS += -g # enable debugging
endif

ifeq ($(PROFILING), 1)
	CFLAGS += -pg # enable profiling
endif

ND=nicksrc
ED=eigensrc
KD=ksrc

NLIB = $(ND)/libnick.a

# ----- phony targets
.PHONY: all clean clobber install

EXE = baseprog convertf mergeit pca smshrink \
	$(ED)/pcatoy $(ED)/smartrel $(ED)/smarteigenstrat \
	$(ED)/twstats $(ED)/eigenstrat $(ED)/eigenstratQTL $(ED)/smartpca

all:    $(EXE)

install:  all
	mv $(EXE) ../bin

clobber:
	rm -f *.o */*.o */*.a $(ND)/*.o $(ED)/*.o $(KD)/*.o
	rm -f $(EXE)
	cd ../bin/ ; rm -f $(notdir $(EXE)) ; cd ../src

clean:
	rm -f *.o core core.* *.o $(ND)/*.o $(ED)/*.o $(KD)/*.o $(EXE)

# ----- build nicksrc/libnick.a
$(NLIB):
	$(MAKE) -C $(ND)

baseprog:	baseprog.o mcio.o egsubs.o admutils.o h2d.o $(ED)/exclude.o $(NLIB)

convertf:	convertf.o mcio.o egsubs.o admutils.o h2d.o $(ED)/exclude.o $(NLIB)

mergeit:	mergeit.o mcio.o admutils.o $(NLIB)

pca:	pca.o $(ED)/eigsubs.o eigensrc/eigx.o $(NLIB)

$(ED)/pcatoy: 	$(ED)/pcatoy.o eigensrc/eigsubs.o eigensrc/eigx.o $(NLIB)

$(ED)/smartrel:	$(ED)/smartrel.o twsubs.o mcio.o qpsubs.o admutils.o egsubs.o regsubs.o \
		$(ED)/eigsubs.o $(ED)/eigx.o $(ED)/smartsubs.o $(NLIB)    

$(ED)/smarteigenstrat:	$(ED)/smarteigenstrat.o mcio.o admutils.o $(NLIB)

$(ED)/twstats:	$(ED)/twstats.o $(NLIB) 

$(ED)/smartpca:	$(ED)/smartpca.o $(ED)/eigsubs.o $(ED)/exclude.o $(ED)/smartsubs.o $(ED)/eigx.o \
		mcio.o qpsubs.o admutils.o egsubs.o regsubs.o gval.o \
		$(NLIB) \
          	$(KD)/kjg_fpca.o $(KD)/kjg_gsl.o 
          	#-lpthread

smshrink:		smshrink.o $(ED)/eigsubs.o $(ED)/exclude.o $(ED)/smartsubs.o $(ED)/eigx.o \
		twsubs.o mcio.o qpsubs.o admutils.o egsubs.o regsubs.o gval.o \
		$(NLIB) \
          	$(KD)/kjg_fpca.o $(KD)/kjg_gsl.o

