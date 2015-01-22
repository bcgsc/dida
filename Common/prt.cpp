#include "prt.h"
#include <cassert>

void qSort(std::vector<long> &a, std::vector<int> &d, int l, int r) {
	long i = l-1, j = r, v = a[d[r]];
	int tmpd;
	if (r <= l) return;
	for (;;) {
		while (a[d[++i]] > v);
		while (v > a[d[--j]]) if (j == l) break;
		if (i >= j) break;
		tmpd = d[i];d[i] = d[j];d[j]=tmpd;
	}
	tmpd = d[i];d[i] = d[r];d[r]=tmpd;
	qSort(a, d, l, i-1);
	qSort(a, d, i+1, r);
}

void tSort(std::vector<int> &a, std::vector<int> &d, int l, int r) {
	int i = l-1, j = r, v = a[d[r]];
	int tmpd;
	if (r <= l) return;
	for (;;) {
		while (a[d[++i]] > v);
		while (v > a[d[--j]]) if (j == l) break;
		if (i >= j) break;
		tmpd = d[i];d[i] = d[j];d[j]=tmpd;
	}
	tmpd = d[i];d[i] = d[r];d[r]=tmpd;
	tSort(a, d, l, i-1);
	tSort(a, d, i+1, r);
}

std::string proLine(std::string line) {
	unsigned i=0;
	std::string temp;
	while(i<line.size()) {
		if (line[i]=='[')
			while(line[i++]!=']');
		if (isdigit(line[i]) || line[i]=='\t' || line[i]==' ')
			temp+=line[i];
		++i;
    }
	return temp;
}

void getFname(const char *filename, std::string &bName, std::string &eName) {
    std::string fName(filename);
    size_t pos = fName.rfind(".");
    if (pos == std::string::npos) {
        bName = fName;
		eName = "";
		return;
	}
    if (pos == 0) {
        bName = "";
		assert(fName.size() > pos+1);
		eName = fName.substr(pos+1, std::string::npos);
		return;
	}
	assert(fName.size() > pos);
    bName = fName.substr(0, pos);
	eName = fName.substr(pos+1, std::string::npos);
}

std::vector< std::vector<int> > getAdj(const char *uName, std::vector<int> &lenArr, long &totLen, const int procRank) {
    std::string bName, eName;
	getFname(uName, bName, eName);
	std::stringstream ast;
	ast << bName << ".adj";

	std::string line;
	std::ifstream adjFile(ast.str().c_str());
	int maxLen;
	while(getline(adjFile, line)) {
		std::istringstream iss(line);
		iss >> maxLen;
	}
	++maxLen;
    
    //std::ofstream imdFile;
	//if (procRank==0) imdFile.open("aln.sam");

	std::vector< std::vector<int> > adjList(maxLen);
	lenArr.resize(maxLen,0);

	adjFile.clear();
	adjFile.seekg(0,adjFile.beg);

	int vertex, uLen, uSite, neighbour;
	totLen=0;
	while(getline(adjFile, line)) {
		std::istringstream iss(proLine(line));
		iss >> vertex >> uLen >> uSite;
		if (procRank==0) std::cout << "@SQ\tSN:"<< vertex << "\tLN:" << uLen << "\n";
		lenArr[vertex]=uLen;
		totLen+=uLen;
		while(iss >> neighbour) adjList[vertex].push_back(neighbour);
	}
	int compNo=0, utigNo=0;
	for (int i=0; i<maxLen; ++i) {
		if (adjList[i].size()) ++compNo;
		if (lenArr[i]) ++utigNo;
	}
	//if (procRank==0) imdFile.close();
	adjFile.close();
	return adjList;
}

std::vector< std::vector<int> > wgetAdj(const char *uName, std::vector<int> &lenArr, long &totLen) {
    std::string bName, eName;
	getFname(uName, bName, eName);
	std::stringstream ast;
	ast << bName << ".adj";
    
	std::string line;
	std::ifstream adjFile(ast.str().c_str());
	int maxLen;
	while(getline(adjFile, line)) {
		std::istringstream iss(line);
		iss >> maxLen;
	}
	++maxLen;
    
    std::ofstream infoFile("maxinf");   
    infoFile << maxLen << "\n";
    infoFile.close();
    
	std::vector< std::vector<int> > adjList(maxLen);
	lenArr.resize(maxLen,0);
    
    std::ofstream alnFile("aln.sam");
    alnFile << "@HD\tVN:0.3\n";
    alnFile << "@PG\tID:DIDA\tPN:DIDA\tVN:0.1.3\n";
    
	adjFile.clear();
	adjFile.seekg(0,adjFile.beg);
    
	int vertex, uLen, uSite, neighbour;
	totLen=0;
	while(getline(adjFile, line)) {
		std::istringstream iss(proLine(line));
		iss >> vertex >> uLen >> uSite;
        alnFile << "@SQ\tSN:"<< vertex << "\tLN:" << uLen << "\n";
		lenArr[vertex]=uLen;
		totLen+=uLen;
		while(iss >> neighbour) adjList[vertex].push_back(neighbour);
	}
	int compNo=0, utigNo=0;
	for (int i=0; i<maxLen; ++i) {
		if (adjList[i].size()) ++compNo;
		if (lenArr[i]) ++utigNo;
	}
	adjFile.close();
    alnFile.close();
	return adjList;
}

std::vector< std::vector<int> > getCom(std::vector< std::vector<int> > &adjList, std::vector<int> &lenArr) {

	int n = adjList.size();
	bool *visited = new bool [n];
	for (int i=0; i< n;++i)
		if (lenArr[i])
			visited[i] = false;
		else
			visited[i] = true;

	//for (int i=0; i< n;++i) std::cout << visited[i] << "\n";


	std::list<int> stack;
	std::vector < std::vector<int> > conComp;

	for (int i = 0; i < n; ++i) {
		if (!visited[i]) {
			std::vector<int> order;
			stack.push_back(i);
			while(!stack.empty()) {
				int top = stack.back();
				stack.pop_back();
				if (visited[top])
					continue;
				visited[top] = true;
				order.push_back(top);
				for (unsigned j=0; j<adjList[top].size(); ++j)
					if (!visited[adjList[top][j]])
						stack.push_back(adjList[top][j]);
			}
			conComp.push_back(order);
		}
	}
	//std::cout << conComp.size() << "\n";
	delete [] visited;
	return conComp;
}

int findBnode(long comLen, std::vector<long> &nodeCap, long totLen) {
	int bNode=-1;
	long minDelta=totLen;
	for (unsigned i=0; i< nodeCap.size(); ++i)
		if (nodeCap[i]-comLen>0 && nodeCap[i]-comLen<minDelta) {
			bNode = i;
			minDelta = nodeCap[i]-comLen;
		}
	if (bNode==-1) {
		long maxCap = nodeCap[0];
		bNode = 0;
		for (unsigned i=1; i< nodeCap.size(); ++i)
			if (nodeCap[i] > maxCap)
				bNode = i;
	}
	return bNode;
}

void fullAloc(std::vector<int> &comp, int bNode, long comLen,
              std::vector<int> &disVec, std::vector<long> &nodeCap) {
	for (unsigned i=0; i<comp.size(); ++i) {
		disVec[comp[i]]=bNode;
	}
	nodeCap[bNode]-=comLen;
}

void partAloc(std::vector<int> &comp, std::vector<long> &nodeCap,
              std::vector<int> &lenArr, std::vector<int> &disVec) {
	for (unsigned i=0; i< comp.size(); ++i) {
		for (unsigned j=0; j< nodeCap.size(); ++j)
			if (lenArr[comp[i]]<=nodeCap[j]) {
				disVec[comp[i]]=j;
				nodeCap[j]-=lenArr[comp[i]];
				break;
			}
	}
}

std::vector<int> compDist(std::vector< std::vector<int> > &conComp, std::vector<int> &lenArr, const long totLen, const int pNum) {
	int totUtig=0;
	std::vector<long> compLen(conComp.size());
	for (unsigned i=0; i< conComp.size();++i) {
		long comLen=0;
		for (unsigned j=0; j< conComp[i].size(); ++j) {
			comLen+=lenArr[conComp[i][j]];
		}
		compLen[i] = comLen;
		totUtig += conComp[i].size();
	}

	std::cerr << "|#utig|=" << totUtig << "\n";

	std::vector<int> compInd(compLen.size());
	for (unsigned i=0; i<compLen.size();++i) compInd[i]=i;
	qSort(compLen, compInd, 0, compLen.size()-1);

	long nodeSize = totLen/pNum + pNum+500;

	std::vector<long> nodeCap(pNum, nodeSize);
	std::vector<int> disVec(lenArr.size(),-1);
	for (unsigned i=0; i< conComp.size();++i) {
		if (compLen[compInd[i]] <= nodeSize) {
			int bNode = findBnode(compLen[compInd[i]], nodeCap, totLen);
			if (bNode!=-1)
				fullAloc(conComp[compInd[i]], bNode,compLen[compInd[i]], disVec, nodeCap);
			else
				partAloc(conComp[compInd[i]], nodeCap, lenArr, disVec);
		}
		else{

			partAloc(conComp[compInd[i]], nodeCap, lenArr, disVec);
		}
	}

    int alcUtig=0;
	for (unsigned i=0; i<disVec.size();++i) {
		if (disVec[i]!=-1)
			++alcUtig;
	}
	std::cerr << "|#aloc|=" << alcUtig << "\n";
    return disVec;
}


/**
 * Return a filename for nth part of the target
 * sequence file.
 */
std::string getPrtFilename(const char *refName, const int procRank)
{
	std::string ref(refName);
	size_t lastDotIndex = ref.find_last_of(".");
	std::ostringstream oss;
	if (lastDotIndex == std::string::npos) {
		oss << ref << "-" << procRank << ".fa";
	} else {
		assert(lastDotIndex < ref.size());
		oss << ref.substr(0, lastDotIndex)
			<< "-" << procRank
			<< ref.substr(lastDotIndex, std::string::npos);
	}
	assert(oss.good());
	return oss.str();
}

void debug_distPar(const char *uName, const int, std::vector<int> &disArr, const int procRank) {

	std::string aPath, uPath, bName, eName;
	getFname(uName, bName, eName);
	std::ifstream uFile(uName);
	std::string line;
	std::ofstream puFile;

	std::stringstream ustm;
	ustm << "debug_mref-" << procRank << ".fa";
	uPath = ustm.str();
	puFile.open(uPath.c_str());


	//begin new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|
	unsigned readId;
	char hChar;
	bool *maskVec = new bool [disArr.size()];
	for (unsigned i=0; i<disArr.size(); ++i) maskVec[i]=true;
	while(getline(uFile, line)) {
		std::istringstream iss(line);
		iss >> hChar >> readId;
		maskVec[readId]=false;
		getline(uFile, line);
	}
	for (unsigned i=0; i<disArr.size(); ++i)
		if (maskVec[i]) disArr[i]=-1;
	delete [] maskVec;
	uFile.clear();
	uFile.seekg(0,uFile.beg);
	//end new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|

	//HMXXXX make all fa rec in one ops
	for (unsigned i=0; i<disArr.size(); ++i)
		if (disArr[i]!=-1) {
			getline(uFile, line);
			if (disArr[i]==procRank-1)
				puFile << line << "\n";
			getline(uFile, line);
			if (disArr[i]==procRank-1)
				puFile << line << "\n";
		}

	puFile.close();

	uFile.close();
}

void ddistPar(const char *uName, const int pNum, std::vector<int> &disArr) {

	std::string aPath, uPath, bName, eName;
	getFname(uName, bName, eName);
	std::ifstream uFile(uName);
	std::string line;
	std::ofstream puFile[pNum];

	for (int i = 0; i < pNum; ++i) {
		uPath = getPrtFilename(uName, i+1);
		puFile[i].open(uPath.c_str());
	}

	//begin new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|
	unsigned readId;
	char hChar;
	bool *maskVec = new bool [disArr.size()];
	for (unsigned i=0; i<disArr.size(); ++i) maskVec[i]=true;
	while(getline(uFile, line)) {
		std::istringstream iss(line);
		iss >> hChar >> readId;
		maskVec[readId]=false;
		getline(uFile, line);
	}
	for (unsigned i=0; i<disArr.size(); ++i)
		if (maskVec[i]) disArr[i]=-1;
	delete [] maskVec;
	uFile.clear();
	uFile.seekg(0,uFile.beg);
	//end new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|

	for (unsigned i=0; i<disArr.size(); ++i)
		if (disArr[i]!=-1) {
			getline(uFile, line);
			puFile[disArr[i]] << line << "\n";
			getline(uFile, line);
			puFile[disArr[i]] << line << "\n";
		}

	for (int i = 0; i < pNum; ++i)
		puFile[i].close();

	uFile.close();
}

void wdistPar(const char *uName, const int pNum, std::vector<int> &disArr) {
    
	std::string aPath, uPath, bName, eName;
	getFname(uName, bName, eName);
	std::ifstream uFile(uName);
	std::string line;
	std::ofstream puFile[pNum];
    
	for (int i = 0; i < pNum; ++i) {
        uPath = getPrtFilename(uName, i+1);
		puFile[i].open(uPath.c_str());
	}
    
	//begin new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|
	unsigned readId;
	char hChar;
	bool *maskVec = new bool [disArr.size()];
	for (unsigned i=0; i<disArr.size(); ++i) maskVec[i]=true;
	while(getline(uFile, line)) {
		std::istringstream iss(line);
		iss >> hChar >> readId;
		maskVec[readId]=false;
		getline(uFile, line);
	}
	for (unsigned i=0; i<disArr.size(); ++i)
		if (maskVec[i]) disArr[i]=-1;
	delete [] maskVec;
	uFile.clear();
	uFile.seekg(0,uFile.beg);
	//end new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|
    
	for (unsigned i=0; i<disArr.size(); ++i)
		if (disArr[i]!=-1) {
			getline(uFile, line);
			puFile[disArr[i]] << line << "\n";
			getline(uFile, line);
			puFile[disArr[i]] << line << "\n";
		}
    
	for (int i = 0; i < pNum; ++i)
		puFile[i].close();
    
	uFile.close();
}

void debug_distTarget(const char *uName, const int pNum, const int procRank) {
	std::ifstream uFile(uName);
	std::vector<int> lenArr;
	std::string hLine, sLine, line;
	long totLen = 0;
	int uLen = 0, minUlen=((unsigned) 1 << 31) -1;
	getline(uFile, line);
	while (getline(uFile, line)) {
		if (line[0] != '>')
			uLen += line.length();
		else {
			if (uLen < minUlen)
				minUlen = uLen;
			totLen += uLen;
			lenArr.push_back(uLen);
			uLen = 0;
		}
	}
	totLen += uLen;
	lenArr.push_back(uLen);
    
	std::cerr << "|target|=" << totLen << "\t #seq=" << lenArr.size() << "\n";
    
	int tLen = lenArr.size();
	std::vector<int> indArr(tLen);
	for (int i=0;i<tLen;++i)
		indArr[i]=i;
	tSort(lenArr, indArr, 0, lenArr.size()-1);
    
	long nodeSize = totLen/pNum + pNum + 500;
    
	std::vector<long> nodeCap(pNum, nodeSize);
	std::vector<int> disArr(lenArr.size(),-1);
    
	for (unsigned i=0; i < lenArr.size(); ++i) {
		if (lenArr[indArr[i]] <= nodeSize) {
			int bNode = findBnode(lenArr[indArr[i]], nodeCap, totLen);
			if (bNode!=-1) {
				disArr[indArr[i]] = bNode;
				nodeCap[bNode]-=lenArr[indArr[i]];
			}
			else {
				std::cerr << "Error in distributing target! size of target is greater than available node cap.\n" << lenArr[indArr[i]] << "\n";
				exit(EXIT_FAILURE);
			}
		}
		else {
			std::cerr << "Error in distributing target! size of target is greater than whole node cap.\n";
			exit(EXIT_FAILURE);
		}
	}
    
	uFile.clear();
	uFile.seekg(0,uFile.beg);
	std::ofstream puFile;

	
    std::stringstream ustm;
	ustm << "debug_mref-" << procRank << ".fa";
    std::string uPath = ustm.str();
	puFile.open(uPath.c_str());
    
	int tIndex=0;
	getline(uFile, line);
    if (disArr[tIndex]==procRank-1)
        puFile << line << "\n";
    
	while (getline(uFile, line)) {
		if (line[0] != '>') {
            if (disArr[tIndex]==procRank-1)
                puFile << line << "\n";
        }
		else {
			++tIndex;
            if (disArr[tIndex]==procRank-1)
                puFile << line << "\n";
		}
	}
	puFile.close();
	uFile.close();
}

void ddistTarget(const char *uName, const int pNum) {
	std::ifstream uFile(uName);
	std::vector<int> lenArr;
	std::string hLine, sLine, line;
	long totLen = 0;
	int uLen = 0, minUlen=((unsigned) 1 << 31) -1;
	getline(uFile, line);
    if(line[0]=='>') {
        std::istringstream hStm(line);
        char sChar;
        hStm >> sChar >> hLine;
        //hLine = line.substr(1,std::string::npos);
    }
	while (getline(uFile, line)) {
		if (line[0] != '>') {
			uLen += line.length();
        }
		else {
			if (uLen < minUlen)
				minUlen = uLen;
			totLen += uLen;
			lenArr.push_back(uLen);
            std::cout << "@SQ\tSN:"<< hLine << "\tLN:" << uLen << "\n";
            std::istringstream hStm(line);
            char sChar;
            hStm >> sChar >> hLine;
            //hLine = line.substr(1,std::string::npos);
			uLen = 0;
		}
	}
	totLen += uLen;
	lenArr.push_back(uLen);
    
    std::cout << "@SQ\tSN:"<< hLine << "\tLN:" << uLen << "\n";

	std::cerr << "|target|=" << totLen << "\t #seq=" << lenArr.size() << "\n";

	int tLen = lenArr.size();
	std::vector<int> indArr(tLen);
	for (int i=0;i<tLen;++i)
		indArr[i]=i;
	tSort(lenArr, indArr, 0, lenArr.size()-1);

	long nodeSize = totLen/pNum + pNum + 500;

	std::vector<long> nodeCap(pNum, nodeSize);
	std::vector<int> disArr(lenArr.size(),-1);

	for (unsigned i=0; i < lenArr.size(); ++i) {
		if (lenArr[indArr[i]] <= nodeSize) {
			int bNode = findBnode(lenArr[indArr[i]], nodeCap, totLen);
			if (bNode!=-1) {
				disArr[indArr[i]] = bNode;
				nodeCap[bNode]-=lenArr[indArr[i]];
			}
			else {
				std::cerr << "Error in distributing target! size of target is greater than available node cap.\n" << lenArr[indArr[i]] << "\n";
				exit(EXIT_FAILURE);
			}
		}
		else {
			std::cerr << "Error in distributing target! size of target is greater than whole node cap.\n";
			exit(EXIT_FAILURE);
		}
	}

	uFile.clear();
	uFile.seekg(0,uFile.beg);
	std::ofstream puFile[pNum];
	for (int i = 0; i < pNum; ++i) {    
        std::string uPath = getPrtFilename(uName, i+1);
        //std::stringstream ustm;
        //ustm << "dmref-" << i+1 << ".fa";
		puFile[i].open(uPath.c_str());
	}

	int tIndex=0;
	getline(uFile, line);
	puFile[disArr[tIndex]] << line << "\n";
	while (getline(uFile, line)) {
		if (line[0] != '>')
			puFile[disArr[tIndex]] << line << "\n";
		else {
			++tIndex;
			puFile[disArr[tIndex]] << line << "\n";
		}
	}

	for (int i = 0; i < pNum; ++i)
		puFile[i].close();
	uFile.close();
}

void wdistTarget(const char *uName, const int pNum) {
	std::ifstream uFile(uName);
	std::vector<int> lenArr;
	std::string hLine, sLine, line;
	long totLen = 0;
	int uLen = 0, minUlen=((unsigned) 1 << 31) -1;
    
    std::ofstream alnFile("aln.sam");
    alnFile << "@HD\tVN:0.3\n";
    alnFile << "@PG\tID:DIDA\tPN:DIDA\tVN:0.1.3\n";
    
	getline(uFile, line);
    if(line[0]=='>') {
        std::istringstream hStm(line);
        char sChar;
        hStm >> sChar >> hLine;
        //hLine = line.substr(1,std::string::npos);
    }
	while (getline(uFile, line)) {
		if (line[0] != '>')
			uLen += line.length();
		else {
			if (uLen < minUlen)
				minUlen = uLen;
			totLen += uLen;
			lenArr.push_back(uLen);
            alnFile << "@SQ\tSN:"<< hLine << "\tLN:" << uLen << "\n";
            std::istringstream hStm(line);
            char sChar;
            hStm >> sChar >> hLine;
			uLen = 0;
		}
	}
	totLen += uLen;
	lenArr.push_back(uLen);
    
    alnFile << "@SQ\tSN:"<< hLine << "\tLN:" << uLen << "\n";
    alnFile.close();
    
	std::cerr << "|target|=" << totLen << "\t #seq=" << lenArr.size() << "\n";
    
    std::ofstream imdFile("maxinf");
	imdFile<<lenArr.size()<<"\n";
	imdFile.close();
        
	int tLen = lenArr.size();
	std::vector<int> indArr(tLen);
	for (int i=0;i<tLen;++i)
		indArr[i]=i;
	tSort(lenArr, indArr, 0, lenArr.size()-1);
    
	long nodeSize = totLen/pNum + pNum + 500;
    
	std::vector<long> nodeCap(pNum, nodeSize);
	std::vector<int> disArr(lenArr.size(),-1);
    
	for (unsigned i=0; i < lenArr.size(); ++i) {
		if (lenArr[indArr[i]] <= nodeSize) {
			int bNode = findBnode(lenArr[indArr[i]], nodeCap, totLen);
			if (bNode!=-1) {
				disArr[indArr[i]] = bNode;
				nodeCap[bNode]-=lenArr[indArr[i]];
			}
			else {
				std::cerr << "Error in distributing target! size of target is greater than available node cap.\n" << lenArr[indArr[i]] << "\n";
				exit(EXIT_FAILURE);
			}
		}
		else {
			std::cerr << "Error in distributing target! size of target is greater than whole node cap.\n";
			exit(EXIT_FAILURE);
		}
	}
    
	uFile.clear();
	uFile.seekg(0,uFile.beg);
	std::ofstream puFile[pNum];
	for (int i = 0; i < pNum; ++i) {
        std::string uPath = getPrtFilename(uName, i+1);
		puFile[i].open(uPath.c_str());
	}
    
	int tIndex=0;
	getline(uFile, line);
	puFile[disArr[tIndex]] << line << "\n";
	while (getline(uFile, line)) {
		if (line[0] != '>')
			puFile[disArr[tIndex]] << line << "\n";
		else {
			++tIndex;
			puFile[disArr[tIndex]] << line << "\n";
		}
	}
    
	for (int i = 0; i < pNum; ++i)
		puFile[i].close();
	uFile.close();
}

inline bool adjExist (const char *uName) {
    std::string bName, eName;
	getFname(uName, bName, eName);
	std::stringstream ast;
	ast << bName << ".adj";
	std::string line;
	std::ifstream adjFile(ast.str().c_str());
	if (adjFile.good()) {
		adjFile.close();
		return true;
	} else {
		adjFile.close();
		return false;
	}
}

int debug_getPrt(const char *uName, const int pNum, const int procRank) {
    clock_t sTime = clock();

	if (adjExist(uName)) {
		std::cerr << "Dynamic target partitioning.\n";
		std::vector<int> lenArr;
		long totLen;
		std::vector< std::vector<int> > adjList = getAdj(uName, lenArr, totLen, procRank);
		std::vector< std::vector<int> > conComp = getCom(adjList, lenArr);
		std::vector<int> disVec = compDist(conComp, lenArr, totLen, pNum);
		debug_distPar(uName, pNum, disVec, procRank);
	}
	else {
		std::cerr << "Static target partitioning.\n";
		debug_distTarget(uName, pNum, procRank);
	}

	std::cerr << "Running time: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
	return 0;
}

int dgetPrt(const char *uName, const int pNum, const int procRank) {
    clock_t sTime = clock();

	if (adjExist(uName)) {
		std::cerr << "Dynamic target partitioning.\n";
		std::vector<int> lenArr;
		long totLen;
		std::vector< std::vector<int> > adjList = getAdj(uName, lenArr, totLen, procRank);
		std::vector< std::vector<int> > conComp = getCom(adjList, lenArr);
		std::vector<int> disVec = compDist(conComp, lenArr, totLen, pNum);
		ddistPar(uName, pNum, disVec);
	}
	else {
		std::cerr << "Static target partitioning.\n";
		ddistTarget(uName, pNum);
	}

	std::cerr << "Running time: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
	return 0;
}

int wgetPrt(const char *uName, const int pNum) {
    clock_t sTime = clock();
    
	if (adjExist(uName)) {
		std::cerr << "Dynamic target partitioning.\n";
		std::vector<int> lenArr;
		long totLen;
		std::vector< std::vector<int> > adjList = wgetAdj(uName, lenArr, totLen);
		std::vector< std::vector<int> > conComp = getCom(adjList, lenArr);
		std::vector<int> disVec = compDist(conComp, lenArr, totLen, pNum);
		wdistPar(uName, pNum, disVec);
	}
	else {
		std::cerr << "Static target partitioning.\n";
		wdistTarget(uName, pNum);
	}
    
	std::cerr << "Running time: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
	return 0;
}