all: dida

dida: dida.cpp
	g++-49 -I/home/hmohamadi/openmpi-1.6.5/include -pthread -L/home/hmohamadi/openmpi-1.6.5/lib -lmpi_cxx -lmpi -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ICommon -O3 -fopenmp -o bin/dida dida-3step.cpp Common/prt.cpp Common/Uncompress.cpp Common/SignalHandler.cpp Common/Fcontrol.cpp -Wall

dida-3step: dida-3step.cpp
	g++-49 -I/home/hmohamadi/openmpi-1.6.5/include -pthread -L/home/hmohamadi/openmpi-1.6.5/lib -lmpi_cxx -lmpi -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ICommon -O3 -fopenmp -o bin/dida dida-3step.cpp Common/prt.cpp Common/Uncompress.cpp Common/SignalHandler.cpp Common/Fcontrol.cpp -Wall

dida-working: dida-working.cpp
	g++-49 -I/home/hmohamadi/openmpi-1.6.5/include -pthread -L/home/hmohamadi/openmpi-1.6.5/lib -lmpi_cxx -lmpi -ldl -lm -Wl,--export-dynamic -lrt -lnsl -lutil -lm -ICommon -O3 -fopenmp -o bin/dida dida-working.cpp -Wall
#	mpic++ -ICommon -O3 -fopenmp -o bin/dida dida.cpp -ldl -Wall

intel: dida.cpp
	mpiicpc -mt_mpi -ICommon -O3 -openmp -o bin/dida dida.cpp -ldl -Wall

#dida: dida.cpp
#	mpic++ -ICommon -O3 -fopenmp -o bin/dida dida.cpp Common/Uncompress.cpp Common/SignalHandler.cpp Common/Fcontrol.cpp -ldl -Wall

#idida: dida.cpp
#	mpiicpc -mt_mpi -ICommon -O3 -openmp -o bin/dida dida.cpp Common/Uncompress.cpp Common/SignalHandler.cpp Common/Fcontrol.cpp -ldl -Wall

clean:
	rm bin/dida
