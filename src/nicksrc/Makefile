CFLAGS = -I../../include -D_GNU_SOURCE

ifeq ($(DEBUG), 1)
	CFLAGS += -g # enable debugging
endif

ifeq ($(PROFILING), 1)
	CFLAGS += -pg # enable profiling
endif

NLIB = libnick.a
NLIBOBJS = gauss.o gds.o getpars.o linsubs.o sortit.o statsubs.o strsubs.o vsubs.o xsearch.o twtable.o

.PHONY: all clean

all: $(NLIB)

clean:
	rm *.a *.o

# ----- build nicksrc/libnick.a
$(NLIB): $(NLIBOBJS)
	ar -r $@ $^
