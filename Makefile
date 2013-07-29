CXX=g++
MPICXX=mpiCC
CPPFLAGS=-c
OPTFLAGS=-O3
LIBPATH=-ICommon
LDFLAGS=-fopenmp

EXEC1=dida
EXEC2=dsp

all: $(EXEC1) $(EXEC2)

SRCS1=dida.cpp partition.cpp merge.cpp

SRCS2=dsp.cpp Common/HashManager.cpp Common/BloomFilter.cpp Common/city.cpp

$(EXEC1):$(SRCS1)
	$(MPICXX) $(OPTFLAGS) -o $@ $^

$(EXEC2):$(SRCS2)
	$(CXX) $(OPTFLAGS) $(LDFLAGS) $(LIBPATH) -o $@ $^
    
clean:
	rm -rf *.o

cleanest: clean
	rm -rf $(EXEC1) $(EXEC2)
