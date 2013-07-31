#include <mpi.h>
#include <cstdlib>
#include <getopt.h>
#include <algorithm>
#include <iterator>
#include "partition.h"
#include "merge.h"


#define PROGRAM "dida"

static const char VERSION_MESSAGE[] =
PROGRAM " Version 0-1.1 \n"
"Written by Hamid Mohamadi.\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY... TARGET\n"
"Map the sequences of the files QUERY to those of the file TARGET.\n"
"The index files TARGET.fai and TARGET.fm will be used if present.\n"
"\n"
" Options:\n"
"\n"
"  -p, --partition=N       divide reference to N partitions\n"
"  -l, --min-align=N       find matches at least N bp [1]\n"
"  -j, --threads=N         use N parallel threads [1]\n"
"  -b, --bfl=N         use N bp for Bloom filter [1]\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to hmohamadi@bcgsc.ca.\n";

static const char shortopts[] = "p:l:b:j:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "partition", required_argument, NULL, 'p' },
	{ "min-align", required_argument, NULL, 'l' },
    { "bfl", required_argument, NULL, 'b' },
	{ "threads", required_argument, NULL, 'j' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

int main(int argc, char **argv){
    
    int pNum, alnLmer,bflKmer,numThr;
    
    
    std::string commandLine;
    {
        std::ostringstream ss;
		char** last = argv + argc - 1;
		copy(argv, last, std::ostream_iterator<const char *>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}
    
    bool die = false;
	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'j': arg >> numThr; break;
			case 'l': arg >> alnLmer;break;
			case 'b': arg >> bflKmer; break;
			case 'p': arg >> pNum; break;
			case OPT_HELP:
                std::cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
                std::cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
            << (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}
	if (argc - optind != 2) {
        std::cerr << PROGRAM ": missing arguments\n";
		die = true;
	}
    
	if (die) {
        std::cerr << "Try `" << PROGRAM
        << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
    
    const char *readName(argv[argc-2]);
    const char *fasName(argv[argc-1]);
    
	int procSize, procRank, prcrNlen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Get_processor_name(processor_name,&prcrNlen);
	   
	printf("Starting process %d out of %d on %s\n", procRank, procSize, processor_name);

	std::string ubName, ueName;
	getFname(fasName, ubName, ueName);
	std::stringstream adjstm;
	adjstm << ubName << ".adj";

	std::stringstream indstm;
	std::stringstream indstm_tmp;
	std::stringstream alnstm;
	std::stringstream dspstm;

	int token=0;

	// 1. Partition
	if(procRank==0){
		getPartition((adjstm.str()).c_str(), pNum);
		for(int i=1;i<procSize;++i)MPI_Send(&token,1,MPI_INT,i,0,MPI_COMM_WORLD);
	}
	else MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);

	// 2. Dispatch
	if(procRank==0){
		for(int i=1;i<procSize;++i)MPI_Recv(&token,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);
		for(int i=1;i<procSize;++i)MPI_Send(&token,1,MPI_INT,i,0,MPI_COMM_WORLD);
	}
	else if (procRank==procSize-2){
		indstm<<"abyss-index "<<ubName<<"-"<<procRank<<".fa";
		std::cout<<"abyss-index "<<ubName<<"-"<<procRank<<".fa\n";
		std::system(indstm.str().c_str());
		indstm_tmp<<"abyss-index "<<ubName<<"-"<<procRank+1<<".fa";
		std::cout<<"abyss-index "<<ubName<<"-"<<procRank+1<<".fa\n";
		std::system(indstm_tmp.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
	}
	else if(procRank==procSize-1){
		dspstm<<"dsp "<<fasName<<" "<<bflKmer<<" "<<pNum<<" "<<readName;
		std::system(dspstm.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
	}

	// 3. Index
	else{
		indstm<<"abyss-index "<<ubName<<"-"<<procRank<<".fa";
		std::cout<<"abyss-index "<<ubName<<"-"<<procRank<<".fa\n";
		std::system(indstm.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
	}

	// 4. Align
	if(procRank==0)
		for(int i=1;i<procSize;++i)MPI_Recv(&token,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);
	else{
		alnstm<<"abyss-map --order -l"<<alnLmer<<" -j"<<numThr<< " mreads-"<<procRank<<".fastq "<<ubName<<"-"<<procRank<<".fa > aln-"<<procRank<<".sam";
		std::cout<<"abyss-map --order -l"<<alnLmer<<" -j"<<numThr<< " mreads-"<<procRank<<".fastq "<<ubName<<"-"<<procRank<<".fa > aln-"<<procRank<<".sam\n";
		std::system(alnstm.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
	}
		  
	// 5. Merge
	if(procRank==0)
		getMerge(pNum);

	MPI_Finalize();
	return 0; 
}
