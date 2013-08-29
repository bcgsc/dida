/*
 * Tony's requested version
 * for Q=0 for same CIGAR
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

void getInf(unsigned &maxCont, unsigned &maxRead){
	std::ifstream infoFile("maxinf");
	if(infoFile.good())
		infoFile >> maxCont >> maxRead;
	else{
		std::cerr << "Info file not exist!\n";
		exit(1);
	}
	if(maxCont==0 || maxRead==0){
		std::cerr << "Error in info file!\n";
		exit(1);
	}
	infoFile.close();
}

int main(int argc, const char *argv[]){

    const int pNum = atoi(argv[1]);
    const char *inName = argv[2];

	std::ifstream samFiles[pNum];
	std::string bName("aln"), eName=("sam");

	for (int i = 0; i < pNum; ++i){
		std::stringstream sstm;
		sstm << bName << "-" << i+1 << "." << eName;
		samFiles[i].open(sstm.str().c_str());
	}

	//rVisit vector
	//old data get inf. readfile. SP1.fastq=815832422 SP2.fastq=769338690 SP3.fastq=31809228 CHR14=32621862
	unsigned cntgCount=0,readCount=0;// new_pm_40k=25175698 pm_40k=294073884 trimmed data SP1=514155660 SP2=458740024 SP3=5791684
	getInf(cntgCount, readCount);
	std::cout<<"Max_Cntg_ID="<<cntgCount<< "\nMax_Read_ID="<<readCount<<"\n";

	uint8_t *rVisit = new uint8_t [readCount];
	for(unsigned i=0; i<readCount;++i) rVisit[i]=0;

	 // unitig=160626489 rescaffold pm_40k=168182630 contig=167792650 scaffold=168182630 CHR14_Unitig=1096624 CHR14_Contig=1118721
	uint8_t *cVisit = new uint8_t [cntgCount];
	for(unsigned i=0; i<cntgCount;++i) cVisit[i]=0;

	uint16_t *qVisit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) qVisit[i]=0;


	uint16_t *mVisit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) mVisit[i]=0;

	uint16_t *s1Visit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) s1Visit[i]=0;

	uint16_t *s2Visit = new uint16_t [readCount];
	for(unsigned i=0; i<readCount;++i) s2Visit[i]=0;



	//Preprocessing BEGIN

	//Read and Write Head Line of sam files
	std::string HDline, PGline;
	for (int i = 0; i < pNum; ++i){
		getline(samFiles[i], HDline);
		getline(samFiles[i], PGline);
	}

	std::string line;
	std::string SQ1;
	char sn1,sn2,sn3;
	unsigned sqId;
	unsigned readId, pairId, bitFg, refId, readPos;
	uint16_t mVal, s1Val, s2Val, rQual, cgV1, cgV2, cgV3;
	char pSign, cgC1,cgC2, cgC3;
	for (int i = 0; i < pNum; ++i){
		while(getline(samFiles[i], line)) {
			if(line[0]!='@') break;
			std::istringstream iss(line);
			iss>>SQ1>>sn1>>sn2>>sn3>>sqId;
			cVisit[sqId]=(uint8_t)(i+1);
		}
		do{
			std::istringstream iss(line);
			cgV1=cgV2=cgV3=0;
			iss>>readId>>bitFg>>refId>>readPos>>rQual>>cgV1>>cgC1>>cgV2>>cgC2>>cgV3>>cgC3;
			//new changes for selecting best sam record according to CIGAR
			if(bitFg!=4){
				if(cgC1=='M'){
					s1Val=0;
					mVal=cgV1;
					s2Val = (cgC2=='S')?cgV2:0;
				}
				else{
					s1Val=cgV1;
					mVal=cgV2;
					s2Val = (cgC3=='S')?cgV3:0;
				}

				if(mVal>mVisit[readId]){
					rVisit[readId]=(uint8_t)(i+1);
					qVisit[readId]=rQual;
					mVisit[readId]=mVal;
					s1Visit[readId]=s1Val;
					s2Visit[readId]=s2Val;

				}
				else if(mVal==mVisit[readId]&&s1Val==s1Visit[readId]&&s2Val==s2Visit[readId])
						qVisit[readId]=0;

			}
		}while(getline(samFiles[i], line));
	}
	unsigned sum=0;
	for (unsigned i = 0; i < readCount; ++i)
		if(rVisit[i])++sum;
	std::cerr << "|reads|" << sum << std::endl;
	sum=0;
	for (unsigned i = 0; i < cntgCount; ++i)
		if(cVisit[i])++sum;
	std::cerr << "|cntig|" << sum << std::endl;
	//Preprocessing END

	//Finalizing BEGIN
	std::stringstream sstm;
	sstm<<bName<<"."<<eName;
	std::ofstream comFile(sstm.str().c_str());

	for (int i = 0; i < pNum; ++i){
		samFiles[i].clear();
		samFiles[i].seekg(0,samFiles[i].beg);
		//Skipping two header lines from each sam File
		getline(samFiles[i], HDline);
		getline(samFiles[i], PGline);
	}
	comFile << HDline << "\n" << PGline << "\n";
	//Writing head records
	for (unsigned i=0; i<cntgCount; ++i){
		if(cVisit[i]){
			getline(samFiles[cVisit[i]-1],line);
			comFile<<line<<"\n";
		}
	}
	delete [] cVisit;
	delete [] mVisit;
	delete [] s1Visit;
	delete [] s2Visit;
	//Writing read records
	std::string cigStr, rNext,seqStr, qualStr;
	unsigned pNext;
	int tLen;

	// Renaming the read ID to original one
	std::string readHead, readSeq, readDir, readQual, readHD;
    char hChar;
	std::ifstream inFile(inName);

	for (unsigned i=0; i<readCount; ++i){
		getline(inFile, readHead);
		getline(inFile, readSeq);
		getline(inFile, readDir);
		getline(inFile, readQual);
        std::istringstream issHD(readHead);
        issHD>>hChar>>readHD;
//		readHD = readHead.substr(1, std::string::npos);
		if(rVisit[i]){
			while(getline(samFiles[rVisit[i]-1],line)){
				std::istringstream iss(line);
				iss>>readId>>bitFg>>refId>>readPos>>rQual>>cigStr>>rNext>>pNext>>tLen>>seqStr>>qualStr;
				if(readId==i){
					//new change to import Qual=0 to repetitive aln
					if(qVisit[i]!=0)
						comFile<<readHD<<"\t"<<bitFg<<"\t"<<refId<<"\t"<<readPos<<"\t"
						<<rQual<<"\t"<<cigStr<<"\t"<<rNext<<"\t"<<pNext<<"\t"<<tLen<<"\t"<<seqStr<<"\t"<<qualStr<<"\n";
					else
						comFile<<readHD<<"\t"<<bitFg<<"\t"<<refId<<"\t"<<readPos<<"\t"
						<<"0\t"<<cigStr<<"\t"<<rNext<<"\t"<<pNext<<"\t"<<tLen<<"\t"<<seqStr<<"\t"<<qualStr<<"\n";
					break;
				}
			}
		}
		else{
			comFile<<readHD<<"\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		}
	}
	inFile.close();
	//Finalizing END

	for (int i = 0; i < pNum; ++i) samFiles[i].close();
	comFile.close();
	delete [] rVisit;
	delete [] qVisit;

	return 0;
}
