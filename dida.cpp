#include <mpi.h>
#include <cstdlib>

#include "partition.h"
#include "merge.h"

int main(int argc, char *argv[]){
	    
	int procSize, procRank, prcrNlen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Get_processor_name(processor_name,&prcrNlen);
	   
	printf("Starting process %d out of %d on %s\n", procRank, procSize, processor_name);

	int pNum = atoi(argv[1]); // Number of Partitions
	int alnLmer = atoi(argv[2]); // Alignment censecutive match length
	int bflKmer=atoi(argv[3]); // Bloom filter kmer length
	const char *readName = argv[4]; // Name of read file in fastq format
	const char *fasName = argv[5]; // Name of Refrence, Contig file in fasta format

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
		indstm<<"/home/hmohamadi/abyss/bin/abyss-index "<<ubName<<"-"<<procRank<<".fa";
		std::cout<<"/home/hmohamadi/abyss/bin/abyss-index "<<ubName<<"-"<<procRank<<".fa\n";
		std::system(indstm.str().c_str());
		indstm_tmp<<"/home/hmohamadi/abyss/bin/abyss-index "<<ubName<<"-"<<procRank+1<<".fa";
		std::cout<<"/home/hmohamadi/abyss/bin/abyss-index "<<ubName<<"-"<<procRank+1<<".fa\n";
		std::system(indstm_tmp.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
	}
	else if(procRank==procSize-1){
		dspstm<<"/scratch/hmohamadi/DIDA/bin/dsp "<<fasName<<" "<<bflKmer<<" "<<pNum<<" "<<readName;
		std::system(dspstm.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
	}

	// 3. Index
	else{
		indstm<<"/home/hmohamadi/abyss/bin/abyss-index "<<ubName<<"-"<<procRank<<".fa";
		std::cout<<"/home/hmohamadi/abyss/bin/abyss-index "<<ubName<<"-"<<procRank<<".fa\n";
		std::system(indstm.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
	}

	// 4. Align
	if(procRank==0)
		for(int i=1;i<procSize;++i)MPI_Recv(&token,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);
	else{
		alnstm<<"/home/hmohamadi/abyss/bin/abyss-map --order -l"<<alnLmer<<" -j12 mreads-"<<procRank<<".fastq "<<ubName<<"-"<<procRank<<".fa > aln-"<<procRank<<".sam";
		std::cout<<"/home/hmohamadi/abyss/bin/abyss-map --order -l"<<alnLmer<<" -j12 mreads-"<<procRank<<".fastq "<<ubName<<"-"<<procRank<<".fa > aln-"<<procRank<<".sam\n";
		std::system(alnstm.str().c_str());
		MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
	}
		  
	// 5. Merge
	if(procRank==0)
		getMerge(pNum);

	MPI_Finalize();
	return 0; 
}
