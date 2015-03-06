CFLAGS = -I../../include -I/opt/openblas/include

ifeq ($(DEBUG), 1)
	CFLAGS += -g # enable debugging
endif

ifeq ($(PROFILING), 1)
	CFLAGS += -pg # enable profiling
endif

.PHONY: all clean

all: kjg_gsl.o kjg_fpca.o

clean:
	rm *.o
