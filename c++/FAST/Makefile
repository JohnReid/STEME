CPP = g++
CPPFLAGS = -g -O3 -Winline --param large-function-growth=10000 --param max-inline-insns-single=10000 --param inline-unit-growth=10000 -Wno-deprecated
NR = 0

FFTW = /usr/u/niranjan/fftw-3.0.1/lib 

ifeq ($(NR), 1)
LIBS = -lm
else
LIBS = -lm -L $(FFTW) -lfftw3
endif

BINDIR = .

ME_OBJS = motif_evaluator.o theta_fns.o minimize.o convolution.o fft.o utils.o llr_score.o
ENTROPY_OBJS = entropy-test.cpp entropy-combined.o entropy.o convolution.o HS.o theta_fns.o minimize.o fft.o utils.o llr_score.o
EXES = HS_sFFT HS_csFFT HS_refine_csFFT MEtest

.cpp.o:
	$(CPP) -D NR=$(NR) $(CPPFLAGS) -c $<

depend:	
	$(CPP) $(CPPFLAGS) -MM *.cpp > depend 

clean:
	-rm *~ *.o

MEtest:	me-test.cpp $(ME_OBJS)
	$(CPP) $(CPPFLAGS) me-test.cpp $(ME_OBJS) $(LIBS) -o $(BINDIR)/$@

HS_sFFT:	$(ENTROPY_OBJS)
	$(CPP) $(CPPFLAGS) -D HS_SFFT $(ENTROPY_OBJS) $(LIBS) -o $(BINDIR)/$@

HS_csFFT:	$(ENTROPY_OBJS)
	$(CPP) $(CPPFLAGS) -D HS_CSFFT $(ENTROPY_OBJS) $(LIBS) -o $(BINDIR)/$@

HS_refine_csFFT:	$(ENTROPY_OBJS)
	$(CPP) $(CPPFLAGS) -D HS_REFINE_CSFFT $(ENTROPY_OBJS) $(LIBS) -o $(BINDIR)/$@

lib:	$(ENTROPY_OBJS)
	$(CPP) $(CPPFLAGS) $(ENTROPY_OBJS) $(LIBS) -shared -o libpvalue.a

HS_sFFT_lib:	entropy-test.cpp libpvalue.a
	$(CPP) $(CPPFLAGS) -D HS_SFFT entropy-test.cpp ./libpvalue.a $(LIBS) -o $(BINDIR)/$@

motifeval:	$(ME_OBJS)
	$(CPP) $(CPPFLAGS) $(ME_OBJS) $(LIBS) -shared -o libmotifeval.a

all:	
	make $(EXES)

include depend
