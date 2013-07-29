#include "dispatch.h"

void getFnameb(const char *filename, std::string &bName, std::string &eName){
    std::string fName(filename);
    size_t pos = fName.rfind(".");
    if(pos == std::string::npos){
        bName = fName;
		eName = "";
		return;
	}
    if(pos == 0){
        bName = "";
		eName = fName.substr(pos+1, std::string::npos);
		return;
	}
    bName = fName.substr(0, pos);
	eName = fName.substr(pos+1, std::string::npos);
}

long getInfo2(const char *aName, int k){
	std::string line;
	std::ifstream faFile(aName);
	int uLen;
	long totItm=0;
	while(getline(faFile, line)){
		getline(faFile, line);
		uLen = line.length();
		totItm+=uLen-k+1;
	}
	std::cerr << "|totLen|=" << totItm << std::endl;
	faFile.close();
	return totItm;
}

std::string getMin(const std::string& sSeq){
	std::string fsSeq(sSeq),rsSeq(sSeq);
	int sLen=sSeq.length()-1;
	for(int i=0; i<=sLen;++i)
		switch(fsSeq[i]){
		case 'A':
					rsSeq[sLen-i]='T';
					break;
		case 'a':
					rsSeq[sLen-i]='T';
					fsSeq[i]='A';
					break;
		case 'C':
					rsSeq[sLen-i]='G';
					break;
		case 'c':
					rsSeq[sLen-i]='G';
					fsSeq[i]='C';
					break;
		case 'G':
					rsSeq[sLen-i]='C';
					break;
		case 'g':
					rsSeq[sLen-i]='C';
					fsSeq[i]='G';
					break;
		case 'T':
					rsSeq[sLen-i]='A';
					break;
		case 't':
					rsSeq[sLen-i]='A';
					fsSeq[i]='T';
					break;
		default:	
					break;

		}
	int p=0,hLen=sLen/2;
	while(fsSeq[p]==rsSeq[p]){
		++p;
		if(p>hLen)break;
	}
	if(rsSeq[p]<fsSeq[p])
		return rsSeq;
	return fsSeq;
}

int getDispatch(const char *ufName, const int k, const int pNum, const char *readName){
#ifdef _OPENMP
	double start = omp_get_wtime();
#else
	clock_t sTime = clock();
#endif
	//const char *ufName = argv[1];
	//const int k = atoi(argv[2]);
	//const int pNum = atoi(argv[3]);
	//const char *readName = argv[4];//"/genesis/scratch/hmohamadi/abyss/data/mreads.fastq";
	
	std::string bName, eName;
	getFnameb(ufName, bName, eName);

	//create filters
	HashManager hashMan;
	for (unsigned int i = 0; i < HASH_NUMBER; ++i)
		hashMan.addHashFunction("CityHash64", i);

	long filterSize;
	BloomFilter* myFilters = static_cast<BloomFilter*>(::operator new (sizeof(BloomFilter)*pNum));
	for (int i=0; i<pNum; ++i){
		std::stringstream sstm;
		sstm << bName << "-" << i+1 << "." << eName;
		filterSize = BIT_PER_ITEM*getInfo2((sstm.str()).c_str(), k);
		new (&myFilters[i]) BloomFilter(filterSize, hashMan);
	}
#ifdef _OPENMP
	int tNum = omp_get_max_threads()>pNum?pNum:omp_get_max_threads();
	std::cerr << "Number of threads=" << tNum << std::endl;
	omp_set_num_threads(tNum);
#endif
	
	//load filters
	int i,chunk=1;
#pragma omp parallel for shared(myFilters) private(i) schedule(static,chunk)
	for (i=0; i<pNum; ++i){
		std::stringstream ustm;
		ustm << bName << "-" << i+1 << ".fa";
		std::ifstream uFile((ustm.str()).c_str());
		std::string line;
		int uL;
		while (getline(uFile, line)){
			getline(uFile, line);
			uL= line.length();
			for(int j=0; j< uL -k+1; ++j)
				myFilters[i].insert(getMin(line.substr(j,k)));
		}
		uFile.close();
	}
	std::cerr<< "Loading BF done!" << std::endl;
#ifdef _OPENMP
	std::cerr << "Loading in sec: " << omp_get_wtime() - start << "\n";
#else
	std::cerr << "Running time of loading in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif

	//************************************
	//whole read file
	//output partioned read files initialize
	std::ofstream rdFiles[pNum];
	std::string rbName, reName;
	getFnameb(readName, rbName, reName);
	for (int i = 0; i < pNum; ++i){
		std::stringstream rstm;
		rstm << "mreads-" << i+1 << ".fastq";
		rdFiles[i].open((rstm.str()).c_str());
	}
    std::ifstream readFile(readName);
    std::string readHead, readSeq, readDir, readQual;
    unsigned readId=0;
    int l=0;
    while(getline(readFile, readHead)){
            getline(readFile, readSeq);
            getline(readFile, readDir);
			getline(readFile, readQual);
            l = readSeq.length();
  	    #pragma omp parallel for shared(myFilters,rdFiles) private(i) schedule(static,chunk)
            for(i=0; i<pNum; ++i)
    			for(int j=0; j<l-k+1; j+=k)
    				if(myFilters[i].contains(getMin(readSeq.substr(j,k)))){
						rdFiles[i] << "@" << (readId/2) << "/" << (readId%2+1)  << "\n" << readSeq << "\n"<< readDir << "\n"<< readQual << "\n";
						break;
    				}
            ++readId;
    }
    readFile.close();
    for (int i=0; i<pNum; ++i) rdFiles[i].close();
    std::ofstream imdFile("maxinf", std::ios_base::app);
	imdFile<<readId<<"\n";
	imdFile.close();
	
	// delete filters
	for (int i=0; i<pNum; ++i)
		myFilters[i].~BloomFilter();
	::operator delete (myFilters);
#ifdef _OPENMP
	std::cerr << "Running time in sec: " << omp_get_wtime() - start << "\n";
#else
	std::cerr << "Running time in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
#endif
	return 0;
}

