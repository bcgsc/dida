#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "BloomFilter.h"
#include "HashManager.h"

#define GENOME_SIZE 20000000000
#define BIT_PER_ITEM 8
#define HASH_NUMBER 5

int memory_usage() {
	int mem = 0;
	ifstream proc("/proc/self/status");
	string s;
	while (getline(proc, s), !proc.fail()) {
		if (s.substr(0, 6) == "VmSize") {
			stringstream convert(
					s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

long getInfo(const char *libName, int bmLen){
	std::string readHead, readSeq, readDir, readQual, rName;
	int rLen;
	long totItm=0;
	std::ifstream libFile(libName);
	while(getline(libFile,rName)){
		std::ifstream rFile(rName.c_str());		
		while(getline(rFile, readHead)){
			getline(rFile, readSeq);
			getline(rFile, readDir);
			getline(rFile, readQual);
			rLen = readSeq.length();
			totItm+=rLen-bmLen+1;
		}	
		rFile.close();
	}
	libFile.close();
	std::cerr << "|totLen|=" << totItm << "\n";
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

int main(int argc, const char *argv[]){
	clock_t sTime=clock();
	int memUsage = memory_usage();

	const char *libName = argv[1];
	const int bmLen = atoi(argv[2]);

	//define filter
	HashManager hashMan;
	for (unsigned int i = 0; i < HASH_NUMBER; ++i)
		hashMan.addHashFunction("CityHash64", i);
	//long filterSize=BIT_PER_ITEM*getInfo(libName, bmLen);
    long filterSize=BIT_PER_ITEM*GENOME_SIZE;
	BloomFilter l1Filter(filterSize, hashMan);
    BloomFilter l2Filter(filterSize, hashMan);
	
	//load filter
	std::string readHead, readSeq, readDir, readQual, rName;
	int rLen;
    long c1Bmer=0,c2Bmer=0,c0Bmer=0;
	std::ifstream libFile(libName);
	while(getline(libFile,rName)){
		std::ifstream rFile(rName.c_str());
		while(getline(rFile, readHead)){
			getline(rFile, readSeq);
			getline(rFile, readDir);
			getline(rFile, readQual);
			rLen= readSeq.length();
			for(int j=0; j< rLen-bmLen+1; ++j){
                ++c0Bmer;
                std::string bMer=getMin(readSeq.substr(j,bmLen));
                if(!l1Filter.contains(bMer)){
                    ++c1Bmer;
                    l1Filter.insert(bMer);
                }
                else if(!l2Filter.contains(bMer)){
                    ++c2Bmer;
                    l2Filter.insert(bMer);
                }
            }
		}
		rFile.close();
	}
    libFile.close();
	std::cerr<< "c0Bmers="<<c0Bmer<<" c1Bmers="<<c1Bmer<<" c2Bmer="<<c2Bmer<<"\n";
	std::cerr<< "Memory usage in kbs: " << memory_usage() - memUsage << "\n";
	std::cerr<< "Running time in sec: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
	return 0;
}

