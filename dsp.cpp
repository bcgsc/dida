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
"  -b, --bfl=N             use N bp for Bloom filter windows[1]\n"
"  -h, --hash=N            use N hash functions for Bloom filter [6]\n"
"  -l, --load=path         load filters from path on disk\n"
"  -d, --refdir=path       dir of reference on disk[current dir]\n"
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
    
	/** The size of a k-mer. */
	int bmer = 20;

	/** dir of subtargets. */
	std::string rdir = "./";

	/** single-end library. */
	static int se;

	/** fastq mode dispatch. */
	static int fq;
}


static const char shortopts[] = "p:b:j:d:h:";

enum { OPT_HELP = 1, OPT_VERSION };


static const struct option longopts[] = {
	{ "threads",	required_argument, NULL, 'j' },
	{ "partition",	required_argument, NULL, 'p' },
	{ "bfl",	required_argument, NULL, 'b' },
	{ "hash",	required_argument, NULL, 'h' },
	{ "refdir",	required_argument, NULL, 'd' },
	{ "se",	no_argument, &opt::se, 1 },
	{ "fq",	no_argument, &opt::fq, 1 },
	{ "help",	no_argument, NULL, OPT_HELP },
	{ "version",	no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

size_t getInfo(const char *aName, unsigned k) {
	std::string line;
	std::ifstream faFile(aName);
	size_t totItm=0;
	while (getline(faFile, line)) {
		getline(faFile, line);
		size_t uLen = line.length();
		totItm+=uLen-k+1;
	}
	std::cerr << "|totLen|=" << totItm << std::endl;
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
			case 'p':
				arg >> opt::pnum; break;
			case 'd':
				arg >> opt::rdir; break;
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
    
	if (opt::bmer == 0) {
		std::cerr << PROGRAM ": missing mandatory option `-b'\n";
		die = true;
	}
    
	if (argc - optind != 1) {
		std::cerr << PROGRAM ": missing arguments\n";
		die = true;
	}
    
	if (die) {
		std::cerr << "Try `" << PROGRAM
		<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
    
	const char *libName(argv[argc-1]);
    
    
#ifdef _OPENMP
	double start = omp_get_wtime();
#else
	clock_t sTime = clock();
#endif
    
    
	std::cerr << "Number of hash functions=" << opt::nhash << "\n";
    
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
        sstm << opt::rdir << "mref-" << pIndex+1 << ".fa";
        size_t filterSize = opt::ibits*getInfo((sstm.str()).c_str(), opt::bmer);
        myFilters[pIndex].resize(filterSize);
        std::ifstream uFile((sstm.str()).c_str());
        
        std::string line;
        while (getline(uFile, line)) {
            getline(uFile, line);
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            long uL= line.length();
            for (long j = 0; j < uL -opt::bmer+1; ++j) {
                std::string bMer = line.substr(j,opt::bmer);
                //Begin
                int p=0,hLen=(opt::bmer-1)/2;
                while (bMer[p] == rc(bMer[opt::bmer-1-p])) {
                    ++p;
                    if(p>hLen)break;
                }
                if (bMer[p] > rc(bMer[opt::bmer-1-p])) {
                    for (int lIndex = p, rIndex = opt::bmer-1-p; lIndex<rIndex; ++lIndex,--rIndex) {
                        char tmp = rc(bMer[rIndex]);
                        bMer[rIndex] = rc(bMer[lIndex]);
                        bMer[lIndex] = tmp;
                    }
                }
                //End
                filInsert(myFilters,pIndex,bMer);
            }
        }
        uFile.close();
    }
    
	std::cerr<< "Loading BF done!\n";
#ifdef _OPENMP
	std::cerr << "Loading in sec: " << omp_get_wtime() - start << "\n";
#else
	std::cerr << "Running time of loading in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif
    
	//************************************
	//whole read file
	std::ofstream rdFiles[opt::pnum];
	for (int i = 0; i < opt::pnum; ++i) {
		std::stringstream rstm;
		if (!opt::fq)
			rstm << "mreads-" << i+1 << ".fa";
		else
			rstm << "mreads-" << i+1 << ".fastq";
		rdFiles[i].open((rstm.str()).c_str());
	}
	std::ofstream msFile("lreads.sam");

	if (!opt::se)
		std::cout << "Dispatching paired-end library:\n";
	else
		std::cout << "Dispatching single-end library:\n";

	std::string readHead, readSeq, readDir, readQual, rName;
	unsigned readId=0;
	int l=0;
	std::ifstream libFile(libName);
	while (getline(libFile, rName)) {
		std::ifstream readFile[2];
		readFile[0].open(rName.c_str());
		if (!opt::se) {
			getline(libFile, rName);
			readFile[1].open(rName.c_str());
		}
		int fileNo=0;
		while (getline(readFile[fileNo], readHead)) {
			getline(readFile[fileNo], readSeq);
            std::transform(readSeq.begin(), readSeq.end(), readSeq.begin(), ::toupper);
			getline(readFile[fileNo], readDir);
			getline(readFile[fileNo], readQual);
			l = readSeq.length();
			readHead[0]=':';
			bool dspRead = false;
			#pragma omp parallel for shared(myFilters,rdFiles,dspRead) private(pIndex) schedule(static,chunk)
			for (pIndex=0; pIndex<opt::pnum; ++pIndex) {
				int j=0;
				while (j < l) {
					if (j > l-opt::bmer)
						j=l-opt::bmer;
					std::string bMer = readSeq.substr(j,opt::bmer);
					//Begin
					int p=0,hLen=(opt::bmer-1)/2;
					while (bMer[p] == rc(bMer[opt::bmer-1-p])) {
						++p;
						if(p>hLen)break;
					}
					if (bMer[p] > rc(bMer[opt::bmer-1-p])) {
						for (int lIndex = p, rIndex = opt::bmer-1-p; lIndex<rIndex; ++lIndex,--rIndex) {
							char tmp = rc(bMer[rIndex]);
							bMer[rIndex] = rc(bMer[lIndex]);
							bMer[lIndex] = tmp;
						}
					}
					//END
					if (filContain(myFilters, pIndex, bMer)) {
						#pragma omp critical
							dspRead = true;
						if (!opt::fq)
							rdFiles[pIndex] << ">" << readId << readHead << "\n" << readSeq << "\n";
						else
							rdFiles[pIndex] << "@" << readId << readHead << "\n" << readSeq << "\n"<< readDir << "\n"<< readQual << "\n";
						break;
					}
					j+=opt::bmer;
				}
			}
			if (!dspRead)
				msFile << readId << readHead << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";

			++readId;
			if (!opt::se)
				fileNo = (fileNo + 1)  % 2;
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
	std::cerr << "Running time in sec: " << omp_get_wtime() - start << "\n";
#else
	std::cerr << "Running time in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif
	return 0;
}
