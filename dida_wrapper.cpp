#include "config.h"
#include "Uncompress.h"
#include "FileUtil.h"
#include "prt.h"
#include "mrg.h"

#include "mpi.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <getopt.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <algorithm>
#include <queue>
#include <omp.h>

#define PROGRAM "dida-wrapper"


static const char VERSION_MESSAGE[] =
PROGRAM " Version 0-1.3 \n"
"Written by Hamid Mohamadi.\n"
"Copyright 2014 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY... TARGET\n"
"Dispatch the sequences of the files QUERY based on the Bloom filter of the file TARGET.\n"
"The index files TARGET.fai and TARGET.fm will be used if present.\n"
"\n"
" Options:\n"
"\n"
"  -p, --partition=N       divide reference to N partitions\n"
"  -j, --threads=N         use N parallel threads [partitions]\n"
"  -b, --bfl=N             use N bp for Bloom filter windows[20]\n"
"  -h, --hash=N            use N hash functions for Bloom filter[5]\n"
"  -a, --aligner=STR       use STR as the base aligner\n"
"  -m, --merge=STR         use STR modefor merging partial sam files[best]\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to hmohamadi@bcgsc.ca.\n";

namespace opt {
	/** Number of bits per item. */
	unsigned ibits = 8;
    
	/** The number of parallel threads. */
	static unsigned threads = 0;
    
	/** The number of partitions. */
	static int pnum = -1;
    
	/** The number of hash functions. */
	int nhash = 5;
    
	/** Minimum alignment length. */
	int alen = 20;
    
	/** The size of a b-mer. */
	int bmer = -1;
    
	/** The step size when breaking a read into b-mers. */
	int bmer_step = -1;

    /** The mapper name. */
    std::string mapper = "abyss-map";

    /** The merge mode. */
    std::string merge = "best";

    
	/** single-end library. */
	static int se;
    
	/** fastq mode dispatch. */
	static int fq;
}

static const char shortopts[] = "s:l:b:j:a:m:h:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "threads",	required_argument, NULL, 'j' },
	{ "bmer",	required_argument, NULL, 'b' },
	{ "alen",	required_argument, NULL, 'l' },
	{ "step",	required_argument, NULL, 's' },
	{ "hash",	required_argument, NULL, 'h' },
	{ "aligner",	required_argument, NULL, 'a' },
	{ "merge",	required_argument, NULL, 'm' },
	{ "se",	no_argument, &opt::se, 1 },
	{ "fq",	no_argument, &opt::fq, 1 },
	{ "help",	no_argument, NULL, OPT_HELP },
	{ "version",	no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static const char b2p[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //0
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //1
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //2
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //3
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //4   'A' 'C' 'G'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //5   'T'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //6   'a' 'c' 'g'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //7   't'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //8
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //9
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //10
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //11
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //12
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //13
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //14
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

struct faqRec {
    std::string readHead;
    std::string readSeq;
    std::string readQual;
};

size_t getInfo(const char *aName, unsigned k) {
	std::string line;
	std::ifstream faFile(aName);
    
    getline(faFile, line);
    if (line[0]!='>') {
        std::cerr << "Target file is not in correct format!\n";
        exit(EXIT_FAILURE);
    }
    size_t totItm=0, uLen=0;
    while (getline(faFile, line)) {
		if (line[0] != '>')
			uLen += line.length();
		else {
            if (uLen>=k)
                totItm+=uLen-k+1;
			uLen = 0;
		}
	}
    if (uLen>=k)
        totItm+=uLen-k+1;
    
	std::cerr << "|totLen|=" << totItm << "\n";
	faFile.close();
	return totItm;
}

// MurmurHash2, 64-bit versions, by Austin Appleby
// https://sites.google.com/site/murmurhash/MurmurHash2_64.cpp?attredirects=0
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed ){
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;
    
	uint64_t h = seed ^ (len * m);
    
	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);
    
	while(data != end)
	{
		uint64_t k = *data++;
        
		k *= m;
		k ^= k >> r;
		k *= m;
		
		h ^= k;
		h *= m;
	}
    
	const unsigned char * data2 = (const unsigned char*)data;
    
	switch(len & 7)
	{
        case 7: h ^= uint64_t(data2[6]) << 48;
        case 6: h ^= uint64_t(data2[5]) << 40;
        case 5: h ^= uint64_t(data2[4]) << 32;
        case 4: h ^= uint64_t(data2[3]) << 24;
        case 3: h ^= uint64_t(data2[2]) << 16;
        case 2: h ^= uint64_t(data2[1]) << 8;
        case 1: h ^= uint64_t(data2[0]);
	        h *= m;
	};
    
	h ^= h >> r;
	h *= m;
	h ^= h >> r;
    
	return h;
}

void filInsert(std::vector< std::vector<bool> > &myFilters, const unsigned pn, const std::string &bMer) {
	for (int i=0; i < opt::nhash; ++i)
		myFilters[pn][MurmurHash64A(bMer.c_str(), opt::bmer, i) % myFilters[pn].size()]=true;
}

bool filContain(const std::vector< std::vector<bool> > &myFilters, const unsigned pn, const std::string &bMer) {
	for (int i=0; i < opt::nhash; ++i)
		if (!myFilters[pn][MurmurHash64A(bMer.c_str(), opt::bmer, i) % myFilters[pn].size()])
			return false;
	return true;
}

static inline char rc(char c) {
    switch (c) {
        case 'A':
            return 'T';
            break;
        case 'C':
            return 'G';
            break;
        case 'G':
            return 'C';
            break;
        case 'T':
            return 'A';
            break;
        default:
            break;
    }
    return c;
}

bool is_empty(std::ifstream& pFile) {
    return pFile.peek() == std::ifstream::traits_type::eof();
}

void dida_partition(const int procRank, const char *refName) {
    if (procRank==0) {
        wgetPrt(refName, opt::pnum);
    }
}

void dida_index(const int procRank, const int procSize, const char *refName) {
    if (procRank < procSize && procRank > 0) {
        std::string refPartName = getPrtFilename(refName, procRank);
        if (opt::mapper== "abyss-map") {
            std::ostringstream amap_ind_stm;
            amap_ind_stm << "abyss-index " << refPartName;
            std::cerr<<amap_ind_stm.str()<<"\n";
            std::system(amap_ind_stm.str().c_str());
        }
        else if (opt::mapper== "bwa-mem"){
            std::ostringstream bwa_ind_stm;
            bwa_ind_stm << "bwa index " << refPartName;
            std::cerr<<bwa_ind_stm.str()<<"\n";
            std::system(bwa_ind_stm.str().c_str());
        }
        else if (opt::mapper== "bowtie2"){
            std::ostringstream bow_ind_stm;
            std::string bowRef = refPartName.substr(0,refPartName.rfind("."));
            bow_ind_stm << "bowtie2-build " << refPartName << " " << bowRef;
            std::cerr<<bow_ind_stm.str()<<"\n";
            std::system(bow_ind_stm.str().c_str());
        }
        else if (opt::mapper== "novoalign"){
            std::ostringstream novo_ind_stm;
            std::string novoRef = refPartName.substr(0,refPartName.rfind("."));
            novo_ind_stm << "novoindex -k " << opt::alen << " " << novoRef << " " << refPartName;
            std::cerr<<novo_ind_stm.str()<<"\n";
            std::system(novo_ind_stm.str().c_str());
        }
        else {
            std::cerr << "Error in aligner name;\n";
            exit(EXIT_FAILURE);
        }
    }
}

void getCanon(std::string &bMer) {
    int p=0, hLen=(opt::bmer-1)/2;
    while (bMer[p] == b2p[(unsigned char)bMer[opt::bmer-1-p]]) {
        ++p;
        if(p>=hLen) break;
    }
    if (bMer[p] > b2p[(unsigned char)bMer[opt::bmer-1-p]]) {
        for (int lIndex = p, rIndex = opt::bmer-1-p; lIndex<=rIndex; ++lIndex,--rIndex) {
            char tmp = b2p[(unsigned char)bMer[rIndex]];
            bMer[rIndex] = b2p[(unsigned char)bMer[lIndex]];
            bMer[lIndex] = tmp;
        }
    }
}

std::vector< std::vector<bool> > loadFilter(const char *refName) {
    
#ifdef _OPENMP
	double start = omp_get_wtime();
#else
	clock_t sTime = clock();
#endif
    
#ifdef _OPENMP
	unsigned tNum = omp_get_max_threads()>opt::pnum?opt::pnum:omp_get_max_threads();
    if (opt::threads < tNum && opt::threads > 0)
        tNum = opt::threads;
	std::cerr << "Number of threads=" << tNum << std::endl;
	omp_set_num_threads(tNum);
#endif
	
	int pIndex,chunk=1;
	//begin create filters
	std::vector< std::vector<bool> > myFilters(opt::pnum);
    
    std::cerr << "Loading filters ...\n";
#pragma omp parallel for shared(myFilters) private(pIndex) schedule(static,chunk)
    for (pIndex=0; pIndex<opt::pnum; ++pIndex) {
        std::string refPartName = getPrtFilename(refName, pIndex+1);
        size_t filterSize = opt::ibits*getInfo(refPartName.c_str(), opt::bmer);
        myFilters[pIndex].resize(filterSize);
        std::ifstream uFile(refPartName.c_str());
        
        std::string pline,line;
        getline(uFile, pline);
        while (getline(uFile, pline)) {
            if (pline[0] != '>')
                line += pline;
            else {
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                long uL= line.length();
                for (long j = 0; j < uL -opt::bmer+1; ++j) {
                    std::string bMer = line.substr(j,opt::bmer);
                    getCanon(bMer);
                    filInsert(myFilters,pIndex,bMer);
                }
                line.clear();
            }
        }
        std::transform(line.begin(), line.end(), line.begin(), ::toupper);
        long uL= line.length();
        for (long j = 0; j < uL -opt::bmer+1; ++j) {
            std::string bMer = line.substr(j,opt::bmer);
            getCanon(bMer);
            filInsert(myFilters,pIndex,bMer);
        }        
        uFile.close();
    }
    
	std::cerr<< "Loading BF done!\n";
#ifdef _OPENMP
	std::cerr << "Loading in sec: " << omp_get_wtime() - start << "\n";
#else
	std::cerr << "Running time of loading in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif
    
    return myFilters;
}

void dispatchRead(const char *libName, const std::vector< std::vector<bool> > &myFilters) {
#ifdef _OPENMP
    double start = omp_get_wtime();
#else
    clock_t sTime = clock();
#endif

    size_t buffSize = 4000000;
    std::ofstream rdFiles[opt::pnum];
	for (int i = 0; i < opt::pnum; ++i) {
		std::stringstream rstm;
		if (!opt::fq) rstm << "mreads-" << i+1 << ".fa";
		else rstm << "mreads-" << i+1 << ".fastq";
		rdFiles[i].open((rstm.str()).c_str());
	}
	std::ofstream msFile("lreads.sam");
    size_t fileNo=0, readId=0;
    std::string readHead, readSeq, readDir, readQual, rName;
    std::ifstream libFile(libName);
    while (getline(libFile, rName)) {
        std::ifstream readFile[2];
        readFile[0].open(rName.c_str());
        if (!opt::se) {
			getline(libFile, rName);
			readFile[1].open(rName.c_str());
		}
        bool readValid=true;
        while(readValid) {
            readValid=false;
            // set up readBuff
            std::vector<faqRec> readBuffer; // fixed-size to improve performance
            while (getline(readFile[fileNo], readHead)) {
                getline(readFile[fileNo], readSeq);
                std::transform(readSeq.begin(), readSeq.end(), readSeq.begin(), ::toupper);
                getline(readFile[fileNo], readDir);
                getline(readFile[fileNo], readQual);
                readHead[0]=':';
                faqRec rRec;
                std::ostringstream hstm;
                if(!opt::fq) hstm<<">"<<readId<<readHead;
                else hstm<<"@"<<readId<<readHead;
                rRec.readHead = hstm.str();
                rRec.readSeq = readSeq;
                rRec.readQual = readQual;
                readBuffer.push_back(rRec);
                if(!opt::se) fileNo=(fileNo + 1)%2;
                ++readId;
                if(readBuffer.size()==buffSize) break;
            }
            if(readBuffer.size()==buffSize) readValid=true;
            
            //dispatch buffer
            int pIndex;
            std::vector<bool> dspRead(buffSize,false);
#pragma omp parallel for shared(readBuffer,myFilters,rdFiles,dspRead) private(pIndex)
            for (pIndex=0; pIndex<opt::pnum; ++pIndex) {
                for(size_t bIndex = 0; bIndex<readBuffer.size(); ++bIndex) {
                    faqRec bRead = readBuffer[bIndex];
                    size_t readLen = bRead.readSeq.length();
                    for (size_t j=0; j <= readLen-opt::bmer; j+=opt::bmer_step) {
                        std::string bMer = bRead.readSeq.substr(j,opt::bmer);
                        getCanon(bMer);
                        if (filContain(myFilters, pIndex, bMer)) {
#pragma omp critical
                            dspRead[bIndex] = true;
                            if (!opt::fq)
                                rdFiles[pIndex] << bRead.readHead << "\n" << bRead.readSeq << "\n";
                            else
                                rdFiles[pIndex] << bRead.readHead << "\n" << bRead.readSeq << "\n+\n"<< bRead.readQual << "\n";
                            break;
                        }
                    }
                    
                }
            } // end dispatch buffer
            for(size_t bIndex = 0; bIndex<readBuffer.size(); ++bIndex) {
                if(!dspRead[bIndex]) msFile << readBuffer[bIndex].readHead.substr(1,std::string::npos) << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
            }
        }
        readFile[0].close();
        if (!opt::se)
            readFile[1].close();
        
    }
    libFile.close();
    msFile.close();
    for (int pIndex=0; pIndex<opt::pnum; ++pIndex)
    	rdFiles[pIndex].close();
    std::ofstream imdFile("maxinf", std::ios_base::app);
	imdFile<<readId<<"\n";
	imdFile.close();
#ifdef _OPENMP
    std::cerr << "Running time for dispatcher in sec: " << omp_get_wtime() - start << "\n";
#else
    std::cerr << "Running time for dispatcher in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif
}

void dida_dispatch(const int procRank, const char *libName, const char *refName) {
   	if (procRank==0) {
        std::vector< std::vector<bool> > myFilters = loadFilter(refName);
        dispatchRead(libName, myFilters);
    }
}

void dida_align(const int procRank, const int procSize, const char *refName) {
    if (procRank < procSize && procRank > 0) {
        std::string rfType;
        std::ostringstream readfile_stm;
        (!opt::fq)?rfType = ".fa":rfType = ".fastq";
        readfile_stm << "mreads-"<<procRank<<rfType;
        std::string refPartName = getPrtFilename(refName, procRank);
        
        if (opt::mapper== "abyss-map") {
            std::ostringstream amap_aln_stm;
            amap_aln_stm<<"abyss-map --order -j"<<omp_get_max_threads()<<" -l"<<opt::alen<<
            " "<<readfile_stm.str()<<" "<<refPartName<<" > aln-"<<procRank<<".sam";
            std::cerr<<amap_aln_stm.str()<<"\n";
            std::system(amap_aln_stm.str().c_str());
        }
        else if (opt::mapper== "bwa-mem"){
            std::ostringstream bwa_aln_stm;
            bwa_aln_stm << "bwa mem -t " <<omp_get_max_threads()<<" -k"<<opt::alen<<" "<<
            refPartName<<" "<<readfile_stm.str()<<" > aln-"<<procRank<<".sam";
            std::cerr<<bwa_aln_stm.str()<<"\n";
            std::system(bwa_aln_stm.str().c_str());
        }
        else if (opt::mapper== "bowtie2"){
            std::ostringstream bow_aln_stm;
            std::string bowRef = refPartName.substr(0,refPartName.rfind("."));
            bow_aln_stm << "bowtie2-align -f -p " <<omp_get_max_threads()<<
            " -x "<<bowRef<<" -U "<<readfile_stm.str()<< " -S aln-"<<procRank<<".sam";
            std::cerr<<bow_aln_stm.str()<<"\n";
            std::system(bow_aln_stm.str().c_str());
        }
        else if (opt::mapper== "novoalign"){
            std::ostringstream novo_aln_stm;
            std::string novoRef = refPartName.substr(0,refPartName.rfind("."));
            novo_aln_stm << "novoalign -o Sync -o SAM -f " << readfile_stm.str() <<
            " -d "<<novoRef<<" > aln-"<<procRank<<".sam";
            std::cerr<<novo_aln_stm.str()<<"\n";
            std::system(novo_aln_stm.str().c_str());
        }
        else {
            std::cerr << "Error in aligner name;\n";
            exit(EXIT_FAILURE);
        }
    }
}

void dida_merge(const int procRank, const char *refName) {
    if (procRank==0) {
        call_merger(opt::pnum, opt::mapper, opt::merge);
        for (int i=1; i <= opt::pnum; ++i) {
            std::ostringstream alnfile_stm;
            alnfile_stm << "aln-"<<i<<".sam";
            remove(alnfile_stm.str().c_str());
        }
        remove("lreads.sam");
        remove("maxinf");
        std::cout << "dida finished successfully!\n";
    }
    else {
        std::string rfType;
        std::ostringstream readfile_stm;
        (!opt::fq)?rfType = ".fa":rfType = ".fastq";
        readfile_stm << "mreads-"<<procRank<<rfType;
        std::string refPartName = getPrtFilename(refName, procRank);
        remove(readfile_stm.str().c_str());
        remove(refPartName.c_str());
        
        if (opt::mapper=="abyss-map") {
            std::ostringstream ind_stm;
            ind_stm << refPartName << ".fai";
            remove(ind_stm.str().c_str());
           
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << refPartName << ".fm";
            remove(ind_stm.str().c_str());
        }
        else if(opt::mapper=="bwa-mem") {
            std::ostringstream ind_stm;
            ind_stm << refPartName << ".amb";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << refPartName << ".ann";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << refPartName << ".bwt";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << refPartName << ".pac";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << refPartName << ".fa";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << refPartName << ".sa";
            remove(ind_stm.str().c_str());
        }
        else if(opt::mapper=="bowtie2") {
            std::string bowRef = refPartName.substr(0,refPartName.rfind("."));
            
            std::ostringstream ind_stm;
            ind_stm << bowRef << ".1.bt2";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << bowRef << ".2.bt2";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << bowRef << ".3.bt2";
            remove(ind_stm.str().c_str());
            
            ind_stm.str("");
            ind_stm.clear();
            ind_stm << bowRef << ".4.bt2";
            remove(ind_stm.str().c_str());

            ind_stm.str("");
            ind_stm.clear();
            ind_stm << bowRef << ".rev.1.bt2";
            remove(ind_stm.str().c_str());

            ind_stm.str("");
            ind_stm.clear();
            ind_stm << bowRef << ".rev.2.bt2";
            remove(ind_stm.str().c_str());
        }
        else if(opt::mapper=="novoalign") {
            std::string novoRef = refPartName.substr(0,refPartName.rfind("."));
            remove(novoRef.c_str());
        }
        
    }
}

void dida_barrier(const int procRank, const int procSize) {
    int token = 0;
    MPI_Status status;
    if (procRank==0) {
        for (int i=1; i<=procSize-1; ++i) MPI_Send(&token,1,MPI_INT,i,0,MPI_COMM_WORLD);
        for (int i=1; i<=procSize-1; ++i) MPI_Recv(&token,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);
    }
    else {
        MPI_Recv(&token,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        MPI_Send(&token,1,MPI_INT,0,0,MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    bool die = false;
    std::string blPath;
	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?':
				die = true; break;
			case 'j':
				arg >> opt::threads; break;
			case 'b':
				arg >> opt::bmer; break;
			case 'l':
				arg >> opt::alen; break;
			case 's':
				arg >> opt::bmer_step; break;
			case 'a':
				arg >> opt::mapper; break;
			case 'm':
				arg >> opt::merge; break;
			case 'h':
				arg >> opt::nhash; break;
                
			case OPT_HELP:
				std::cerr << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cerr << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			std::cerr << PROGRAM ": invalid option: `-"
			<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}
    
	if (opt::alen <= 1) {
		std::cerr << PROGRAM ": alignment length must at least 2.\n";
		die = true;
	}
    
	if (argc - optind < 2) {
		std::cerr << PROGRAM ": missing arguments\n";
		die = true;
	}
	if (die) {
		std::cerr << "Try `" << PROGRAM
        << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
    
	if (opt::bmer <= 0)
		opt::bmer = 3 * opt::alen / 4;
    
	if (opt::bmer_step <= 0)
		opt::bmer_step = opt::alen - opt::bmer + 1;
    
	const char *libName(argv[argc-2]);
	const char *refName(argv[argc-1]);
    
	int procSize, procRank, provided;
    
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    
    if (opt::pnum < 1)
		opt::pnum = procSize - 1;
    
	if (procRank == 0) {
		if (provided < MPI_THREAD_FUNNELED) {
			std::cerr << PROGRAM ": MPI_THREAD_SERIALIZED support required.\n"
            << "Install an MPI library which supports this level of thread safety.";
			exit(EXIT_FAILURE);
		}
		if (procSize < 3) {
			std::cerr << PROGRAM ": must specify at least 3 processes with mpirun\n";
			exit(EXIT_FAILURE);
		}
        
		std::cerr << "pnum=" << opt::pnum << " bmer=" << opt::bmer
        << " alen=" << opt::alen
        << " bmer_step=" << opt::bmer_step << std::endl;
	}
    
    double wTime = MPI_Wtime(); /* Wallclock time*/
    
    dida_partition(procRank, refName);
    dida_barrier(procRank, procSize);
    //MPI_Barrier(MPI_COMM_WORLD);
	dida_index(procRank, procSize, refName);
	dida_dispatch(procRank, libName, refName);
    dida_barrier(procRank, procSize);
    //MPI_Barrier(MPI_COMM_WORLD);
	dida_align(procRank, procSize, refName);
    dida_barrier(procRank, procSize);
    //MPI_Barrier(MPI_COMM_WORLD);
	dida_merge(procRank, refName);
    
    wTime = MPI_Wtime() - wTime;
    if (procRank==0)
		fprintf(stderr, "Wallclock time is %.2f seconds\n", wTime);
    
    MPI_Finalize();
	return 0;
}

