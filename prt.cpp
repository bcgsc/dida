#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <time.h>

void qSort(std::vector<long> &a, std::vector<int> &d, int l, int r){
	long i = l-1, j = r, v = a[d[r]];
	int tmpd;
	if (r <= l) return;
	for (;;){
		while (a[d[++i]] > v);
		while (v > a[d[--j]]) if (j == l) break;
		if (i >= j) break;
		tmpd = d[i];d[i] = d[j];d[j]=tmpd;
	}
	tmpd = d[i];d[i] = d[r];d[r]=tmpd;
	qSort(a, d, l, i-1);
	qSort(a, d, i+1, r);
}

std::string proLine(std::string line){
	unsigned i=0;
	std::string temp;
	while(i<line.size()){
		if(line[i]=='[')
			while(line[i++]!=']');
		if(isdigit(line[i]) || line[i]=='\t' || line[i]==' ')
			temp+=line[i];
		++i;
    }
	return temp;
}

void getInfo(const char *aName){
	std::string line;
	std::ifstream adjFile(aName);
	int vertex, uLen, uSite;
	long totLen=0;
	while(getline(adjFile, line)){
		std::istringstream iss(proLine(line));
		iss >> vertex >> uLen >> uSite;
		totLen+=uLen;
	}
    
	std::cout << "|totLen|=" << totLen << "\n";
	adjFile.close();
}

void prtVec(std::vector< std::vector<int> > &myVec){
	for(unsigned i=0; i< myVec.size();++i)
		if(myVec[i].size()){
			for (unsigned j=0; j<myVec[i].size(); ++j)
				std::cout << myVec[i][j] << "\t";
			std::cout << "\n";
		}
}

void getFname(const char *filename, std::string &bName, std::string &eName){
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

std::vector< std::vector<int> > getAdj(const char *uName, std::vector<int> &lenArr, long &totLen){
    std::string bName, eName;
	getFname(uName, bName, eName);
	std::stringstream ast;
	ast << bName << ".adj";

	std::string line;
	std::ifstream adjFile(ast.str().c_str());
	int maxLen;
	while(getline(adjFile, line)){
		std::istringstream iss(line);
		iss >> maxLen;
	}
	++maxLen;
	std::ofstream imdFile("maxinf");
	imdFile<<maxLen<<"\n";
	imdFile.close();
    
	std::vector< std::vector<int> > adjList(maxLen);
	lenArr.resize(maxLen,0);
	
	adjFile.clear();
	adjFile.seekg(0,adjFile.beg);
	
	int vertex, uLen, uSite, neighbour;
	totLen=0;
	while(getline(adjFile, line)){
		std::istringstream iss(proLine(line));
		iss >> vertex >> uLen >> uSite;
		lenArr[vertex]=uLen;
		totLen+=uLen;
		while(iss >> neighbour) adjList[vertex].push_back(neighbour);
	}
	int compNo=0, utigNo=0;
	for(int i=0; i<maxLen; ++i){
		if(adjList[i].size()) ++compNo;
		if(lenArr[i]) ++utigNo;
	}
    
	adjFile.close();
	return adjList;
}

std::vector< std::vector<int> > getCom(std::vector< std::vector<int> > &adjList, std::vector<int> &lenArr){
    
	int n = adjList.size();
	bool *visited = new bool [n];
	for(int i=0; i< n;++i)
		if(lenArr[i])
			visited[i] = false;
		else
			visited[i] = true;
	
	//for(int i=0; i< n;++i) std::cout << visited[i] << "\n";
	
	
	std::list<int> stack;
	std::vector < std::vector<int> > conComp;
	
	for(int i = 0; i < n; ++i){
		if(!visited[i]){
			std::vector<int> order;
			stack.push_back(i);
			while(!stack.empty()){
				int top = stack.back();
				stack.pop_back();
				if(visited[top])
					continue;
				visited[top] = true;
				order.push_back(top);
				for (unsigned j=0; j<adjList[top].size(); ++j)
					if(!visited[adjList[top][j]])
						stack.push_back(adjList[top][j]);
			}
			conComp.push_back(order);
		}
	}
	//std::cout << conComp.size() << "\n";
	delete [] visited;
	return conComp;
}

int findBnode(long comLen, std::vector<long> &nodeCap, long totLen){
	int bNode=-1;
	long minDelta=totLen;
	for(unsigned i=0; i< nodeCap.size(); ++i)
		if(nodeCap[i]-comLen>0 && nodeCap[i]-comLen<minDelta){
			bNode = i;
			minDelta = nodeCap[i]-comLen;
		}
	return bNode;
}

void fullAloc(std::vector<int> &comp, int bNode, long comLen,
              std::vector<int> &disVec, std::vector<long> &nodeCap){
	for(unsigned i=0; i<comp.size(); ++i){
		disVec[comp[i]]=bNode;
	}
	nodeCap[bNode]-=comLen;
}

void partAloc(std::vector<int> &comp, std::vector<long> &nodeCap,
              std::vector<int> &lenArr, std::vector<int> &disVec){
	for(unsigned i=0; i< comp.size(); ++i){
		for(unsigned j=0; j< nodeCap.size(); ++j)
			if(lenArr[comp[i]]<=nodeCap[j]){
				disVec[comp[i]]=j;
				nodeCap[j]-=lenArr[comp[i]];
				break;
			}
	}
}

std::vector<int> compDist(std::vector< std::vector<int> > &conComp, std::vector<int> &lenArr, const long totLen, const int pNum){
	int totUtig=0;
	std::vector<long> compLen(conComp.size());
	for(unsigned i=0; i< conComp.size();++i){
		long comLen=0;
		for(unsigned j=0; j< conComp[i].size(); ++j){
			comLen+=lenArr[conComp[i][j]];
		}
		compLen[i] = comLen;
		totUtig += conComp[i].size();
	}
    
	std::cout << "|#utig|=" << totUtig << "\n";
    
	std::vector<int> compInd(compLen.size());
	for(unsigned i=0; i<compLen.size();++i) compInd[i]=i;
	qSort(compLen, compInd, 0, compLen.size()-1);
    
	long nodeSize = totLen/pNum + pNum+500;
    
	std::vector<long> nodeCap(pNum, nodeSize);
	std::vector<int> disVec(lenArr.size(),-1);
	for(unsigned i=0; i< conComp.size();++i){
		if(compLen[compInd[i]] <= nodeSize){
			int bNode = findBnode(compLen[compInd[i]], nodeCap, totLen);
			if(bNode!=-1)
				fullAloc(conComp[compInd[i]], bNode,compLen[compInd[i]], disVec, nodeCap);
			else
				partAloc(conComp[compInd[i]], nodeCap, lenArr, disVec);
		}
		else{
            
			partAloc(conComp[compInd[i]], nodeCap, lenArr, disVec);
		}
	}

    int alcUtig=0;
	for(unsigned i=0; i<disVec.size();++i){
		if(disVec[i]!=-1)
			++alcUtig;
	}
	std::cout << "|#aloc|=" << alcUtig << "\n";
    return disVec;
}

void distPar(const char *uName, const int pNum, std::vector<int> &disArr){
    
	std::string aPath, uPath, bName, eName;
	getFname(uName, bName, eName);
	//std::stringstream ust;
	//ust << bName << ".fa";
	//std::string upt = ust.str();
	//const char *uName = upt.c_str();
	//std::cout << uName << "\n";
    
	//std::ifstream aFile(aName);
	std::ifstream uFile(uName);
	std::string line;
	//std::ofstream paFile[pNum];
	std::ofstream puFile[pNum];
    
	for (int i = 0; i < pNum; ++i){
		std::stringstream astm;
		std::stringstream ustm;
		//astm << bName << "-" << i+1 << ".adj";
		ustm << bName << "-" << i+1 << ".fa";
		//aPath = astm.str();
        uPath = ustm.str();
		//paFile[i].open(aPath.c_str());
		puFile[i].open(uPath.c_str());
	}
	
	//begin new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|
	unsigned readId;
	char hChar;
	bool *maskVec = new bool [disArr.size()];
	for(unsigned i=0; i<disArr.size(); ++i) maskVec[i]=true;
	while(getline(uFile, line)){
		std::istringstream iss(line);
		iss >> hChar >> readId;
		maskVec[readId]=false;
		getline(uFile, line);
	}
	for(unsigned i=0; i<disArr.size(); ++i)
		if(maskVec[i]) disArr[i]=-1;
	delete [] maskVec;
	uFile.clear();
	uFile.seekg(0,uFile.beg);
	//end new patch for tony: .fa and .adj have different sizes. |.adj| > |.fa|
    
	for(unsigned i=0; i<disArr.size(); ++i)
		if(disArr[i]!=-1){
			getline(uFile, line);
			puFile[disArr[i]] << line << "\n";
			getline(uFile, line);
			puFile[disArr[i]] << line << "\n";
		}
    
	for (int i = 0; i < pNum; ++i) puFile[i].close();
    
	//aFile.close();
	uFile.close();
}

int main(int argc, const char *argv[]){
	clock_t sTime = clock();
	const char *uName = argv[1]; // name of unitig, contig file in .adj format
	const int pNum = atoi(argv[2]); // nubmer of partitions
	
	std::vector<int> lenArr;
	long totLen;
	std::vector< std::vector<int> > adjList = getAdj(uName, lenArr, totLen);
	std::vector< std::vector<int> > conComp = getCom(adjList, lenArr);
    
	/*std::cout << "|utigs|=" << totLen << "\n";
	std::cout << "|#comp|=" << conComp.size() << "\n";
	int snglton=0;
	for(unsigned i=0; i<conComp.size();++i)
		if(conComp[i].size()==1) ++snglton;
	std::cout << "|#ccom|=" << conComp.size()-snglton << "\n";
	std::cout << "|#sgtn|=" << snglton << "\n";*/
    
	//here begin
	std::vector<int> disVec = compDist(conComp, lenArr, totLen, pNum);
	// end best fit
	
	distPar(uName, pNum, disVec);
	
	std::cout << "Running time: " << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
	return 0;
}
