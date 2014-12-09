#include "config.h"
#include "Uncompress.h"
#include "FileUtil.h"
#include "MurmurHash2.h"
#include "prt.h"

#include <mpi.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <getopt.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <algorithm>
#include <queue>
#include <omp.h>
#include <signal.h>

#define PROGRAM "dida-mpi"

using namespace std;

const unsigned READ = 0;
const unsigned WRITE = 1;

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
"  -j, --threads=N         use N parallel threads [partitions]\n"
"  -l, --alen=N            the minimum alignment length [20]\n"
"  -b, --bmer=N            size of a bmer [alen/2]\n"
"  -s, --step=N            step size used when breaking a query sequence into bmers [bmer]\n"
"  -h, --hash=N            use N hash functions for Bloom filter [6]\n"
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

	/** single-end library. */
	static int se;

	/** fastq mode dispatch. */
	static int fq;
}

static const char shortopts[] = "s:l:b:j:d:h:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "threads",	required_argument, NULL, 'j' },
	{ "bmer",	required_argument, NULL, 'b' },
	{ "alen",	required_argument, NULL, 'l' },
	{ "step",	required_argument, NULL, 's' },
	{ "hash",	required_argument, NULL, 'h' },
	{ "refdir",	required_argument, NULL, 'd' },
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

size_t getInfo(const char *aName, unsigned k){
	std::string line;
	std::ifstream faFile(aName);
	size_t totItm=0;
	size_t fasta_lines = 0;
	while(getline(faFile, line)){
		getline(faFile, line);
		size_t uLen = line.length();
		totItm+=uLen-k+1;
		++fasta_lines;
	}
	if (fasta_lines == 0) {
		cerr << PROGRAM ": error: failed to read from file `" << aName << "'." << endl;
	}
	assert(fasta_lines > 0);
	assert(totItm > 0); // Require that there be at least one item in the file
	std::cerr << "|totLen|=" << totItm << std::endl;
	faFile.close();
	return totItm;
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

struct samRec {
    unsigned SamOrd;
    std::string SamQn;
    int SamFg;
    std::string SamRf;
    int SamPs;
    int SamMq;
    std::string SamCr;
    std::string SamRn;
    int SamPn;
    int SamTl;
    std::string SamSq;
    std::string SamPh;
    int SamSc;
};

samRec recLoad(std::string& line, int partNum){
    std::istringstream iss(line);
    samRec c;
    char hSep;
    iss>>c.SamOrd>>hSep>>c.SamQn>>c.SamFg>>c.SamRf>>c.SamPs>>c.SamMq>>c.SamCr;
    c.SamRn="*";c.SamPn=0;c.SamTl=0;c.SamSq="*";c.SamPh="*";
    c.SamSc=partNum;
    return c;
}

std::ostream& operator<<(std::ostream& os, const samRec& c){
    os<<c.SamQn<<"\t"<<c.SamFg<<"\t"<<c.SamRf<<"\t"<<c.SamPs<<"\t"<<c.SamMq<<"\t"<<c.SamCr<<"\t*\t0\t0\t*\t*";
    return os;
}

bool is_empty(std::ifstream& pFile) {
    return pFile.peek() == std::ifstream::traits_type::eof();
}

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

void dida_partition(const int procRank, const char *refName) {

    if (procRank == 0) {

        bool allPrtFilesExist = true;
        for (int i = 1; i <= opt::pnum; ++i) {
            if (!fileExists(getPrtFilename(refName, i))) {
                allPrtFilesExist = false;
                break;
            }
        }

        if (allPrtFilesExist) {
            std::cerr << "rank " << procRank << ": "
                "skipping partition stage, all partition files already exist: ";
            for (int i = 1; i <= opt::pnum; ++i) {
                if (i != 1)
                    std::cerr << ", ";
                std::cerr << "'" << getPrtFilename(refName, i) << "'";
            }
            std::cerr << std::endl;
        } else {
            for (int i = 1; i <= opt::pnum; ++i)
                getPrt(refName, opt::pnum, i);
        }

    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void dida_index(const int procRank, const int procSize, const char *refName) {
    if (procRank < procSize-1 && procRank > 0) {
        std::string refPartName = getPrtFilename(refName, procRank);
        std::string indexName = refPartName + ".fm";
        if (fileExists(indexName)) {
            std::cerr << "rank " << procRank << ": "
                "skipping index stage, index file '"
                << indexName << "' already exists"
                << std::endl;
        } else {
            std::ostringstream oss;
            oss << "abyss-index " << refPartName;
            assert(oss.good());
            std::string cmd = oss.str();
            std::cerr << "Rank " << procRank << ": "
                << "calling system(\"" << cmd << "\")"
                << std::endl;
            int result = system(cmd.c_str());
            if (result != 0) {
                std::cerr << "command failed: " << cmd
                    << " (exit status: " << result << ")"
                    << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}

std::vector< std::vector<bool> > loadFilter(const char *refName) {
    int pIndex,chunk=1;
    //begin create filters
    std::vector< std::vector<bool> > myFilters(opt::pnum);
    
    std::cerr << "Loading filters ...\n";
#pragma omp parallel for shared(myFilters) private(pIndex) schedule(static,chunk)
    for (pIndex=0; pIndex<opt::pnum; ++pIndex){
        std::string refPartName = getPrtFilename(refName, pIndex+1);
        size_t filterSize = opt::ibits*getInfo(refPartName.c_str(), opt::bmer);
        myFilters[pIndex].resize(filterSize);
        
        std::ifstream uFile(refPartName.c_str());
        std::string line;
        while (getline(uFile, line)){
            getline(uFile, line);
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            long uL= line.length();
            for(long j = 0; j < uL -opt::bmer+1; ++j) {
				assert((int)line.size() > j);
//				assert((int)line.size() > j + opt::bmer);
                std::string bMer = line.substr(j,opt::bmer);
                //Begin
                int p=0,hLen=(opt::bmer-1)/2;
                while (bMer[p] == b2p[(unsigned char)bMer[opt::bmer-1-p]]) {
                    ++p;
                    if(p>hLen)break;
                }
                if (bMer[p] > b2p[(unsigned char)bMer[opt::bmer-1-p]]) {
                    for (int lIndex = p, rIndex = opt::bmer-1-p; lIndex<rIndex; ++lIndex,--rIndex) {
                        char tmp = b2p[(unsigned char)bMer[rIndex]];
                        bMer[rIndex] = b2p[(unsigned char)bMer[lIndex]];
                        bMer[lIndex] = tmp;
                    }
                }
                //End
                filInsert(myFilters, pIndex, bMer);
            }
        }
        uFile.close();
    }
    std::cerr<< "Loading BF done!\n";
    return myFilters;
}

void dispatchRead(const int procSize, const vector<string>& inFiles, const std::vector< std::vector<bool> > &myFilters) {
    int pIndex,chunk=1;
    int readLen=0, recLen=0;
    std::cerr << "Number of hash functions=" << opt::nhash << "\n";
    
    unsigned tNum = omp_get_max_threads()>opt::pnum?opt::pnum:omp_get_max_threads();
    if (opt::threads < tNum && opt::threads > 0)
        tNum = opt::threads;
    std::cerr << "Number of threads=" << tNum << std::endl;
    omp_set_num_threads(tNum);
    
    
    if (!opt::se)
        std::cerr << "Dispatching paired-end library:\n";
    else
        std::cerr << "Dispatching single-end library:\n";
    
    std::string readHead, readSeq, readDir, readQual, rName;
    unsigned readId=0, notDsp=0;
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        std::ifstream readFile[2];
		assert(file_i < inFiles.size());
        readFile[0].open(inFiles[file_i].c_str());
        if (!opt::se) {
			++file_i;
			assert(file_i < inFiles.size());
            readFile[1].open(inFiles[file_i].c_str());
        }
        int fileNo=0;
        while (getline(readFile[fileNo], readHead)) {
            getline(readFile[fileNo], readSeq);
            std::transform(readSeq.begin(), readSeq.end(), readSeq.begin(), ::toupper);
            getline(readFile[fileNo], readDir);
            getline(readFile[fileNo], readQual);
            readLen = readSeq.length();
            readHead[0]=':';
            bool dspRead = false;
#pragma omp parallel for shared(myFilters, dspRead) private(pIndex) schedule(static,chunk)
            for(pIndex=1; pIndex<=opt::pnum; ++pIndex) {
                int j=0;
                while (j < readLen) {
                    if (j > readLen-opt::bmer)
                        j=readLen-opt::bmer;
					assert((int)readSeq.size() > j);
                    std::string bMer = readSeq.substr(j,opt::bmer);
                    //Begin preprocess k-mer
					int p=0,hLen=(opt::bmer-1)/2;
					while (bMer[p] == b2p[(unsigned char)bMer[opt::bmer-1-p]]) {
						++p;
						if(p>hLen)break;
					}
					if (bMer[p] > b2p[(unsigned char)bMer[opt::bmer-1-p]]) {
						for (int lIndex = p, rIndex = opt::bmer-1-p; lIndex<rIndex; ++lIndex,--rIndex) {
							char tmp = b2p[(unsigned char)bMer[rIndex]];
							bMer[rIndex] = b2p[(unsigned char)bMer[lIndex]];
							bMer[lIndex] = tmp;
						}
					}
					//END                    
                    if (filContain(myFilters, pIndex-1, bMer)) {
#pragma omp critical
                        {
                            dspRead = true;
                            std::stringstream hstm;
                            if (!opt::fq)
                                hstm << ">" << readId << readHead << "\n" << readSeq;
                            else
                                hstm << "@" << readId << readHead << "\n" << readSeq << "\n+\n" << readQual;
                            std::string readRec = hstm.str();
                            recLen = readRec.length();

                            MPI_Send(&recLen, 1, MPI_INT, pIndex, 0, MPI_COMM_WORLD);
                            MPI_Send(&readRec[0], recLen+1, MPI_CHAR, pIndex, 0, MPI_COMM_WORLD);
                            MPI_Send(&pIndex, 1, MPI_INT, procSize-1, 0, MPI_COMM_WORLD);
                        }
                        break;
                    }
                    j+=opt::bmer_step;
                }
            }
            if (!dspRead)
#pragma omp critical
            {
                std::stringstream hstm;
                if (!opt::fq)
                    hstm << ">" << readId << readHead << "\n" << readSeq;
                else
                    hstm << "@" << readId << readHead << "\n" << readSeq << "\n+\n" << readQual;
                std::string readRec = hstm.str();
                recLen = readRec.length();
                
                ++notDsp;
                int rndIndex = notDsp%opt::pnum + 1;
                MPI_Send(&recLen, 1, MPI_INT, rndIndex, 0, MPI_COMM_WORLD);
                MPI_Send(&readRec[0], recLen+1, MPI_CHAR, rndIndex, 0, MPI_COMM_WORLD);
                MPI_Send(&rndIndex, 1, MPI_INT, procSize-1, 0, MPI_COMM_WORLD);
            }
            ++readId;
            if (!opt::se)
                fileNo = (fileNo + 1)  % 2;
        }
        readFile[0].close();
        if (!opt::se)
            readFile[1].close();
    }
    int minImpi = -1;
#pragma omp critical
    {
        for(int pIndex = 1; pIndex<=opt::pnum; ++pIndex)
            MPI_Send(&minImpi, 1, MPI_INT, pIndex, 0, MPI_COMM_WORLD);
        MPI_Send(&minImpi, 1, MPI_INT, procSize-1, 0, MPI_COMM_WORLD);
    }
}

void dida_dispatch(const int procRank, const int procSize, const vector<string>& inFiles, const char *refName) {
   	if (procRank==0) {
        std::vector< std::vector<bool> > myFilters = loadFilter(refName);
        dispatchRead(procSize, inFiles, myFilters);
    }
}

/* Redirect STDOUT to input of fd. */
void stdout_to_pipe(int fd[])
{
	close(fd[READ]);
	dup2(fd[WRITE], STDOUT_FILENO);
	close(fd[WRITE]);
}

/* Redirect the output of fd to STDIN. */
void pipe_to_stdin(int fd[])
{
	close(fd[WRITE]);
	dup2(fd[READ], STDIN_FILENO);
	close(fd[READ]);
}

void dida_align(const int procRank, const int procSize, const char *refName) {
    if (procRank < procSize-1 && procRank > 0) {
    	MPI_Status status;
        int  recLen=0;
        int pid, fd[2];
        int fd2[2];
        if (pipe(fd) == -1 || pipe(fd2) == -1) {
            perror("pipe failed");
            exit(1);
        }

        if ((pid = fork()) < 0) {
            perror("fork failed");
            exit(2);
        }

        if(pid != 0) {
            stdout_to_pipe(fd);
            pipe_to_stdin(fd2);

            #pragma omp parallel sections
            {
                // receive reads from dispatch process
                #pragma omp section
                {
                	for(;;) {
                       #pragma omp critical
					   MPI_Recv(&recLen, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
					   if(recLen==-1) break;
					   char *readbuf = new char [recLen + 1];
                       #pragma omp critical
					   MPI_Recv(readbuf,recLen+1,MPI_CHAR,0,0,MPI_COMM_WORLD,&status);
                       std::string myRead(readbuf, readbuf+recLen);
					   delete [] readbuf;
					   fprintf(stdout, "%s\n",myRead.c_str());
				   }
					fclose(stdout);
                }

                // send SAM alignments to merge process
                #pragma omp section
                {
                    char *line = NULL;
                    size_t lineByte = 0;
                    ssize_t lineLen;
                    while ((lineLen = getline(&line, &lineByte, stdin)) != -1) {
                        if (line[0]!='@')
                        #pragma omp critical
                        {
							assert(lineLen > 0);
                            MPI_Send(&lineLen, 1, MPI_INT, procSize-1, 0, MPI_COMM_WORLD);
                    		MPI_Send(line, lineLen+1, MPI_CHAR, procSize-1, 0, MPI_COMM_WORLD);
                        }
                    }
                    free(line);
                    lineLen=-1;
                    #pragma omp critical
                    MPI_Send(&lineLen, 1, MPI_INT, procSize-1, 0, MPI_COMM_WORLD);
                }
            }
        }
        else {
            stdout_to_pipe(fd2);
            pipe_to_stdin(fd);

            /*execlp("malign", "malign", (char *) 0);
            perror("malign failed");
            exit(3);*/

            /*std::ifstream reffile(amapstm.str().c_str());
            if (!reffile)
                std::cerr<<"no file\n";
            else
                std::cerr<<"file exists\n";
            if (is_empty(reffile))
                std::cerr << "reffile is empty\n";
            else
                std::cerr<<"file in not empty\n";*/

            std::string refPartName = getPrtFilename(refName, procRank);
            std::ostringstream amap_j_stm, amap_l_stm;
            amap_j_stm << "-j" << opt::threads;
            amap_l_stm << "-l" << opt::alen;

            execlp("abyss-map", "abyss-map", "--order", amap_j_stm.str().c_str(), amap_l_stm.str().c_str(), "-", refPartName.c_str(), (char *)0);

            /*std::ostringstream bow_x_stm, bow_p_stm;
            bow_x_stm << "-xmref-" << procRank;
            bow_p_stm << "-p" << omp_get_max_threads()-2;
            execlp("bowtie2-align", "bowtie2-align", "-f", "--reorder", bow_p_stm.str().c_str(), bow_x_stm.str().c_str(), "-U-", (char *)0);*/
        }
    }
}

samRec getSam(std::vector <samRec> &recBuffer){
	int mVal, s1Val, s2Val, cgV1, cgV2, cgV3;
	int rVisit_max=0, qVisit_max=0, mVisit_max=0, s1Visit_max=0, s2Visit_max=0;
	char cgC1, cgC2, cgC3;
    for (size_t i = 0; i < recBuffer.size(); ++i) {
    	std::istringstream iss(recBuffer[i].SamCr);
    	cgV1=cgV2=cgV3=0;
    	iss>>cgV1>>cgC1>>cgV2>>cgC2>>cgV3>>cgC3;
    	if (recBuffer[i].SamFg!=4) {
			if  (cgC1 == 'M') {
				s1Val=0;
				mVal=cgV1;
				s2Val = (cgC2=='S')?cgV2:0;
			}
			else {
				s1Val=cgV1;
				mVal=cgV2;
				s2Val = (cgC3=='S')?cgV3:0;
			}
			if (mVal>mVisit_max) {
				rVisit_max=i;
				qVisit_max=recBuffer[i].SamMq;
				mVisit_max=mVal;
				s1Visit_max=s1Val;
				s2Visit_max=s2Val;
			}
			else if(mVal==mVisit_max&&s1Val==s1Visit_max&&s2Val==s2Visit_max)
				qVisit_max=0;
    	}
    }
    samRec c = recBuffer[rVisit_max];
    c.SamMq = qVisit_max;
    return c;
}

void dida_merge(const int procRank, const int procSize) {
	if (procRank==procSize-1) {
		MPI_Status status;
		int isize=0, alnId=0;
		std::vector< samRec > recBuffer;
		unsigned samCounter=0;
		for(;;) {
			MPI_Recv(&alnId, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			if(alnId==-1) break;
			MPI_Recv(&isize, 1, MPI_INT, alnId, 0, MPI_COMM_WORLD, &status);
			assert(isize > 0);
			char *readbuf = new char [isize+1];
			MPI_Recv(readbuf, isize+1, MPI_CHAR, alnId, 0, MPI_COMM_WORLD, &status);
			std::string sambuf(readbuf, readbuf+isize);
			delete [] readbuf;
			samRec cRec = recLoad(sambuf,alnId);
			if (cRec.SamOrd == samCounter) {
				recBuffer.push_back(cRec);
			}
			else {
				std::cout << getSam(recBuffer) << std::endl;
				++samCounter;
				recBuffer.clear();
				recBuffer.push_back(cRec);
			}
		}
		samRec lastRec = getSam(recBuffer);
		std::cout << lastRec << std::endl;

		std::cerr << "dida finished successfully!\n";
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
		opt::bmer = opt::alen / 2;

	if (opt::bmer_step <= 0)
		opt::bmer_step = opt::bmer;

	vector<string> inFiles;
	for (int i = optind; i < argc-1; ++i) {
		string file(argv[i]);
		inFiles.push_back(file);
	}
	const char *refName(argv[argc-1]);

	int procSize, procRank, prcrNlen, provided;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

	/**
	 * The following line was needed to fix failures of the
	 * system() call under OpenMPI 1.8.3. The system() call was
	 * intermittently failing with exit status == -1, errno == 10,
	 * and perror() reporting "no child processes".
	 *
	 * It appears that MPI_Init_thread calls "signal(SIGCHLD, SIG_IGN)"
	 * which instructs the kernel to immediately cleanup completed
	 * child processes (i.e. "zombie processes") without giving the
	 * parent a chance to wait() on the child and obtain the exit
	 * status.  For further background info, see:
	 * http://stackoverflow.com/a/18438562
	 */
	signal(SIGCHLD, SIG_DFL);

	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	if (procRank == 0) {
		if (provided < MPI_THREAD_SERIALIZED) {
			std::cerr << PROGRAM ": MPI_THREAD_SERIALIZED support required.\n"
				<< "Install an MPI library which supports this level of thread safety.";
			exit(EXIT_FAILURE);
		}
		if (procSize < 3) {
			std::cerr << PROGRAM ": must specify at least 3 processes with mpirun\n";
			exit(EXIT_FAILURE);
		}

		std::cerr << "bmer=" << opt::bmer
			<< " alen=" << opt::alen
			<< " bmer_step=" << opt::bmer_step << std::endl;
	}

	if (opt::pnum < 1)
		opt::pnum = procSize - 2;

	MPI_Get_processor_name(processor_name,&prcrNlen);

	fprintf(stderr, "process %d out of %d on %s with thread level %d\n", procRank, procSize, processor_name, provided);

	dida_partition(procRank, refName);
	dida_index(procRank, procSize, refName);
	dida_dispatch(procRank, procSize, inFiles, refName);
	dida_align(procRank, procSize, refName);
	dida_merge(procRank, procSize);

    MPI_Finalize();
	return 0;
}
