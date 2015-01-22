#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <getopt.h>
#include "Uncompress.h"
#ifdef _OPENMP
# include <omp.h>
#endif

#define PROGRAM "dsp"

static const char VERSION_MESSAGE[] =
PROGRAM " Version 1.0.0 \n"
"Written by Hamid Mohamadi.\n"
"Copyright 2014 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY\n"
"Dispatch the sequences of the files QUERY based on the Bloom filter of the file TARGET.\n"
"\n"
" Options:\n"
"\n"
"  -p, --partition=N       divide reference to N partitions\n"
"  -j, --threads=N         use N parallel threads [partitions]\n"
"  -l, --alen=N            the minimum alignment length [20]\n"
"  -b, --bmer=N            size of a bmer [alen/2]\n"
"  -s, --step=N            step size used when breaking a query sequence into bmers [bmer]\n"
"  -h, --hash=N            use N hash functions for Bloom filter [6]\n"
"      --se                single-end library\n"
"      --fq                dispatch reads in fastq format\n"
"      --store             store filters on disk\n"
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
	static int pnum = 1;
    
	/** The number of hash functions. */
	int nhash = 5;
    
	/** Minimum alignment length. */
	int alen = 20;
    
	/** The size of a b-mer. */
	int bmer = -1;
    
	/** The step size when breaking a read into b-mers. */
	int bmer_step = -1;
    
	/** single-end library. */
	static int se;

	/** fastq mode dispatch. */
	static int fq;
}

static const char shortopts[] = "s:l:b:p:j:d:h:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "threads",	required_argument, NULL, 'j' },
	{ "partition",	required_argument, NULL, 'p' },
    { "bmer",	required_argument, NULL, 'b' },
	{ "alen",	required_argument, NULL, 'l' },
	{ "step",	required_argument, NULL, 's' },
	{ "hash",	required_argument, NULL, 'h' },
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
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed ) {
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while (data != end)
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

std::vector< std::vector<bool> > loadFilter() {

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
        std::stringstream sstm;
        sstm << "mref-" << pIndex+1 << ".fa";
        size_t filterSize = opt::ibits*getInfo((sstm.str()).c_str(), opt::bmer);
        myFilters[pIndex].resize(filterSize);
        std::ifstream uFile((sstm.str()).c_str());
        
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
                    //size_t j=0;
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
}

int main(int argc, char** argv) {

#ifdef _OPENMP
	double start = omp_get_wtime();
#else
	clock_t sTime = clock();
#endif
    
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
			case 'p':
				arg >> opt::pnum; break;
			case 'l':
				arg >> opt::alen; break;
			case 's':
				arg >> opt::bmer_step; break;
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
    
	if (argc - optind < 1) {
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
		opt::bmer_step = opt::alen-opt::bmer+1;
    
    std::cerr<<"bmer-step="<<opt::bmer_step<<"\n";
    std::cerr<<"bmer="<<opt::bmer<<"\n";
    std::cerr<<"alen="<<opt::alen<<"\n";
    
	const char *libName(argv[argc-1]);
    
    std::vector< std::vector<bool> > myFilters = loadFilter();
    dispatchRead(libName, myFilters);
    
	

#ifdef _OPENMP
	std::cerr << "Running time in sec: " << omp_get_wtime() - start << "\n";
#else
	std::cerr << "Running time in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif
	return 0;
}
