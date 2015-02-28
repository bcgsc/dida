#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <list>
#include <cstdlib>
#include <stdint.h>
#include <getopt.h>

#define PROGRAM "mrg"

static const char VERSION_MESSAGE[] =
PROGRAM " Version 1.0.0 \n"
"Written by Hamid Mohamadi.\n"
"Copyright 2015 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]\n"
"Merge partial SAM files resluted from alignment step into the final compliant SAM result.\n"
"\n"
" Options:\n"
"\n"
"  -p, --partition=N       value of N partitions\n"
"  -r, --rec=N             report up to N sam records for each query[5]\n"
"  -a, --aligner=NAME      NAME aligner used for alignmet step\n"
"  -m, --mode=NAME         use NAME mode for merging. fast: fast merging, ord: ordered merging, \
mem: coupled with -r to report up to N sam records for each read, best: best quality sam record [mem]\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to hmohamadi@bcgsc.ca.\n";

namespace opt {
    /** Maximum number of sam record per query. */
    unsigned nsam = 5;
}

static const char shortopts[] = "p:a:m:r:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "partition",	required_argument, NULL, 'p' },
	{ "aligner",	required_argument, NULL, 'a' },
	{ "mode",	required_argument, NULL, 'm' },
    { "rec",	required_argument, NULL, 'r' },
	{ "help",	no_argument, NULL, OPT_HELP },
	{ "version",	no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

struct samRec {
    long SamOrd;
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
    int SamPr;
};

void getInf(unsigned &maxCont, unsigned &maxRead) {
	std::ifstream infoFile("maxinf");
	if (infoFile.good())
		infoFile >> maxCont >> maxRead;
	else{
		std::cerr << "Info file not exist!\n";
		exit(1);
	}
	if (maxCont==0 || maxRead==0) {
		std::cerr << "Error in info file!\n";
		exit(1);
	}
	infoFile.close();
}

int memory_usage() {
	int mem = 0;
    std::ifstream proc("/proc/self/status");
	std::string s;
	while (getline(proc, s), !proc.fail()) {
		if (s.substr(0, 6) == "VmSize") {
			std::stringstream convert(
				s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

bool operator>(const samRec& lhs, const samRec& rhs) {
    return lhs.SamOrd > rhs.SamOrd;
}

samRec recLoad(std::string& line, int partNum) {
    std::istringstream iss(line);
    samRec c;
    char hSep;
    iss>>c.SamOrd>>hSep>>c.SamQn>>c.SamFg>>c.SamRf>>c.SamPs>>c.SamMq>>c.SamCr;
    c.SamRn="*";c.SamPn=0;c.SamTl=0;c.SamSq="*";c.SamPh="*";
    c.SamPr=partNum;
    return c;
}

std::ostream& operator<<(std::ostream& os, const samRec& c) {
    os<<c.SamQn<<"\t"<<c.SamFg<<"\t"<<c.SamRf<<"\t"<<c.SamPs<<"\t"<<c.SamMq<<"\t"<<c.SamCr<<"\t*\t0\t0\t*\t*";
    return os;
}

void memMer(const int pNum, const std::string &alignerName) {
    std::string mapStr(alignerName);
	std::ifstream samFiles[pNum+1];

	for (int i = 0; i < pNum; ++i) {
		std::stringstream sstm;
		sstm << "aln-" << i+1 << ".sam";
		samFiles[i].open(sstm.str().c_str());
	}
	samFiles[pNum].open("lreads.sam");

    std::priority_queue< samRec, std::vector<samRec>, std::greater<samRec> > recBuffer;
    std::string line;
    std::ofstream comFile("aln.sam", std::ios_base::app);

    // Skipping @
	for (int i = 0; i < pNum; ++i) {
        while (getline(samFiles[i], line))
            if (line[0]!='@') break;
        recBuffer.push(recLoad(line,i));
    }
    
    for (int dIndex = 1; dIndex <= 3; ++dIndex)
		for (int i = 0; i < pNum+1; ++i)
			if (getline(samFiles[i], line))
				recBuffer.push(recLoad(line,i));

    long psOrd = -1, nsamCount=0;
    std::string psHead;
    bool samVal = true;
	while (!recBuffer.empty()) {
		samRec cRec=recBuffer.top();
		if (cRec.SamOrd != psOrd) {
			if (!samVal)
				comFile<<psHead<<"\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
			else {
				samVal = false;
                nsamCount=0;
            }
		}
		if (cRec.SamFg != 4) {
			samVal = true;
            if (nsamCount < opt::nsam) { //CHANGE
                comFile<<cRec<<"\n";
                ++nsamCount;
            }
		}
		recBuffer.pop();
		if (getline(samFiles[cRec.SamPr], line))
			recBuffer.push(recLoad(line,cRec.SamPr));
		psOrd = cRec.SamOrd;
		psHead = cRec.SamQn;
	}

    for (int i = 0; i < pNum+1; ++i) samFiles[i].close();
	comFile.close();
}

void fstMer(const int pNum, const std::string &alignerName) {
    std::string mapStr(alignerName);
	std::ifstream samFiles[pNum];

	for (int i = 0; i < pNum; ++i) {
		std::stringstream sstm;
		sstm << "aln-" << i+1 << ".sam";
		samFiles[i].open(sstm.str().c_str());
	}

    unsigned cntgCount=0,readCount=0;
    getInf(cntgCount, readCount);

    std::cerr<<"Maximum target ID="<<cntgCount<< "\nTotal number of queries="<<readCount<<"\n";

	std::string line, headSQ;
    // First pass, Reading
	int headEnd[pNum];
	for (int i = 0; i < pNum; ++i) {
		headEnd[i]=0;
        while (getline(samFiles[i], line)) {
            if (line[0]!='@') break;
            ++headEnd[i];
        }
	}

    // Second pass, Writing
    std::ofstream comFile("aln.sam", std::ios_base::app);

    for (int i = 0; i < pNum; ++i) {
		samFiles[i].clear();
		samFiles[i].seekg(0,samFiles[i].beg);
		for (int j = 0; j < headEnd[i]; ++j)
			getline(samFiles[i], line);
	}

    bool inIndex[pNum];
	for (int pIndex = 0; pIndex < pNum; ++pIndex)
		inIndex[pIndex] = true;
    int inCount = pNum;
    while (inCount)
		for (int pIndex = 0; pIndex < pNum; ++pIndex)
	    	if (inIndex[pIndex]) {
	    		if (getline(samFiles[pIndex],line)) {
	    			size_t pos = line.find_first_of(":");
	    			comFile<<line.substr(pos+1, std::string::npos)<<"\n";
	    		}
				else {
					inIndex[pIndex] = false;
					--inCount;
				}
	    	}

    for (int i = 0; i < pNum; ++i) samFiles[i].close();
	comFile.close();
}

void fordMer(const int pNum, const std::string &alignerName) {
    std::string mapStr(alignerName);
	std::ifstream samFiles[pNum+1];

	for (int i = 0; i < pNum; ++i) {
		std::stringstream sstm;
		sstm << "aln-" << i+1 << ".sam";
		samFiles[i].open(sstm.str().c_str());
	}
	samFiles[pNum].open("lreads.sam");

    unsigned cntgCount=0,readCount=0;
    getInf(cntgCount, readCount);

    std::cerr<<"Maximum target ID="<<cntgCount<< "\nTotal number of queries="<<readCount<<"\n";
    
	std::vector< std::vector<uint8_t> > ordList(readCount);

	unsigned readId,bitFg;
    std::string line, readHead, headSQ;

    // First pass, Reading
	int headEnd[pNum];
	for (int i = 0; i < pNum; ++i) {
		headEnd[i]=0;
        while (getline(samFiles[i], line)) {
            if (line[0]!='@') break;
            ++headEnd[i];
        }
		// inserting SAM info into ordList
        do {
            std::istringstream iss(line);
            iss>>readId;
            ordList[readId].push_back((uint8_t)(i+1));
        } while (getline(samFiles[i], line));
	}

    // Second pass, Writing
    std::ofstream comFile("aln.sam", std::ios_base::app);

    for (int i = 0; i < pNum; ++i) {
		samFiles[i].clear();
		samFiles[i].seekg(0,samFiles[i].beg);
		//Discarding @
		for (int j = 0; j < headEnd[i]; ++j)
			getline(samFiles[i], line);
	}
    samFiles[pNum].clear();
	samFiles[pNum].seekg(0,samFiles[pNum].beg);


    char colChar;
    for (unsigned i=0; i<readCount; ++i) {
		bool samVal = false;
		for (unsigned j=0; j<ordList[i].size(); ++j) {
			getline(samFiles[ordList[i][j]-1],line);
			std::istringstream iss(line);
			iss>>readId>>colChar>>readHead>>bitFg;
			if (bitFg!=4) {
				samVal=true;
				size_t pos = line.find_first_of(":");
				comFile<<line.substr(pos+1, std::string::npos)<<"\n";
				//comFile << line << "\n";
			}
		}
		if (!samVal) {
			comFile<<readHead<<"\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		}
    }

    for (int i = 0; i < pNum+1; ++i) samFiles[i].close();
	comFile.close();
}

void bestMer(const int pNum, const std::string &alignerName) {
	std::ifstream samFiles[pNum];

	for (int i = 0; i < pNum; ++i) {
		std::stringstream sstm;
		sstm << "aln-" << i+1 << ".sam";
		samFiles[i].open(sstm.str().c_str());
	}

	unsigned cntgCount=0,readCount=0;
    getInf(cntgCount, readCount);

	std::cerr<<"Maximum target ID="<<cntgCount<< "\nTotal number of queries="<<readCount<<"\n";

	uint8_t *rVisit = new uint8_t [readCount];
	for(unsigned i=0; i<readCount;++i) rVisit[i]=0;

	uint16_t *qVisit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) qVisit[i]=0;

	uint16_t *mVisit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) mVisit[i]=0;

	uint16_t *s1Visit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) s1Visit[i]=0;

	uint16_t *s2Visit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) s2Visit[i]=0;

    unsigned readId, bitFg, readPos;
    std::string line, readHead, headSQ, refId;
	uint16_t mVal, s1Val, s2Val, rQual, cgV1, cgV2, cgV3;
    char colChar, cgC1, cgC2, cgC3;

    // First pass, Reading
	int headEnd[pNum];
	for (int i = 0; i < pNum; ++i) {
		headEnd[i]=0;
        while (getline(samFiles[i], line)) {
            if (line[0]!='@') break;
            ++headEnd[i];
        }
        do {
            std::istringstream iss(line);
            cgV1=cgV2=cgV3=0;
            iss>>readId>>colChar>>readHead>>bitFg>>refId>>readPos>>rQual>>cgV1>>cgC1>>cgV2>>cgC2>>cgV3>>cgC3;
			if (bitFg != 4) {
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
				if (mVal>mVisit[readId]) {
					rVisit[readId]=(uint8_t)(i+1);
					qVisit[readId]=rQual;
					mVisit[readId]=mVal;
					s1Visit[readId]=s1Val;
					s2Visit[readId]=s2Val;
				}
				else if(mVal==mVisit[readId])
					qVisit[readId]=0;
			}
			else if (!rVisit[readId])
				rVisit[readId]=(uint8_t)(i+1);
        } while (getline(samFiles[i], line));
	}
    
    unsigned alignedCount=0;
	for (unsigned i = 0; i < readCount; ++i)
		if(qVisit[i])++alignedCount;
	std::cerr << "Number of the aligned reads: " << alignedCount << "\t" << (alignedCount*100.0/readCount) << "%\n";

    // Second pass, Writing
    std::ofstream comFile("aln.sam", std::ios_base::app);

    for (int i = 0; i < pNum; ++i) {
		samFiles[i].clear();
		samFiles[i].seekg(0,samFiles[i].beg);
		//Discarding @ 
		for (int j = 0; j < headEnd[i]; ++j)
			getline(samFiles[i], line);
	}

	delete [] mVisit;
	delete [] s1Visit;
	delete [] s2Visit;

	std::string cigStr, rNext,seqStr, qualStr;
	unsigned pNext;
	int tLen;

	std::ifstream nullSam("lreads.sam");
	bool pairSign = true;
	if (alignerName=="abyss-map") {
		getline(nullSam, line);
		std::istringstream iss(line);
		iss>> readId >> colChar >> readHead;
		if (readHead[readHead.length()-2] != '/')
			pairSign = false;
		nullSam.clear();
		nullSam.seekg(0,nullSam.beg);
	}

	for (unsigned i=0; i<readCount; ++i) {
		if (rVisit[i]) {
			while (getline(samFiles[rVisit[i]-1],line)) {
				std::istringstream iss(line);
				// initialize all fields before iss >>
				iss>>readId>>colChar>>readHead>>bitFg>>refId>>readPos>>rQual>>cigStr>>rNext>>pNext>>tLen>>seqStr>>qualStr;
				if (readId == i) {
					//new change to import Qual=0 to repetitive aln
					if (qVisit[i]!=0 || bitFg == 4) {
						size_t pos = line.find_first_of(":");
						comFile<<line.substr(pos+1, std::string::npos)<<"\n";
					}
					else
						comFile<<readHead<<"\t"<<bitFg<<"\t"<<refId<<"\t"<<readPos<<"\t"<<"0\t"<<
						cigStr<<"\t"<<rNext<<"\t"<<pNext<<"\t"<<tLen<<"\t"<<seqStr<<"\t"<<qualStr<<"\n";
					break;
				}
			}
		}
		else {
			getline(nullSam, line);
			std::istringstream iss(line);
			iss>> readId >> colChar >> readHead;
			if (pairSign)
				comFile << readHead << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
			else
				comFile << readHead << "/" << (readId % 2 + 1) << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		}
    }

    for (int i = 0; i < pNum; ++i)
    	samFiles[i].close();
    nullSam.close();
	comFile.close();
	delete [] rVisit;
	delete [] qVisit;
}

int main(int argc, char** argv) {
	clock_t sTime = clock();

	bool die = false;
    int pNum=0;
	std::string alignerName;
	std::string runMode;

	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?':
				die = true; break;
			case 'a':
				arg >> alignerName; break;
			case 'p':
				arg >> pNum; break;
			case 'm':
				arg >> runMode; break;
            case 'r':
				arg >> opt::nsam; break;
                
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

	if (pNum == 0) {
		std::cerr << PROGRAM ": missing mandatory option `-p'\n";
		die = true;
	}

	if (alignerName == "") {
		std::cerr << PROGRAM ": missing mandatory option `-a'\n";
		die = true;
	}

	if (runMode == "") {
		std::cerr << PROGRAM ": missing mandatory option `-m'\n";
		die = true;
	}

	if (die) {
		std::cerr << "Try `" << PROGRAM
		<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

    if (runMode == "fast") {
    	std::cerr << "Fast MERGE:\n";
    	fstMer(pNum, alignerName);
    }
    else if (runMode == "mem") {
    	std::cerr << "Minimum memory usage MERGE:\n";
    	memMer(pNum, alignerName);
    }
    else if (runMode == "ord") {
    	std::cerr << "Ordered and fast MERGE:\n";
    	fordMer(pNum, alignerName);
    }
    else if (runMode == "best") {
    	std::cerr << "The best quality sam record MERGE:\n";
    	bestMer(pNum, alignerName);
    }
    else
    	std::cerr << "Error, merge mode not specified! Please choose fast, mord, or ord.\n";

    std::cout << "Running time: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
    return 0;
}
