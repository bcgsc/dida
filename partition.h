#ifndef PARTITION_H
#define PARTITION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <time.h>

void qSort(std::vector<long> &, std::vector<int> &, int , int );

std::string proLine(std::string);

void getInfo(const char *);

void prtVec(std::vector< std::vector<int> > &);

std::vector< std::vector<int> > getAdj(const char *, std::vector<int> &, long &);

std::vector< std::vector<int> > getCom(std::vector< std::vector<int> > &, std::vector<int> &);

int findBnode(long, std::vector<long> &, long);

void fullAloc(std::vector<int> &, int, long, std::vector<int> &, std::vector<long> &);

void partAloc(std::vector<int> &, std::vector<long> &, std::vector<int> &, std::vector<int> &);

void getFname(const char *, std::string &, std::string &);

std::vector<int> compDist(std::vector< std::vector<int> > &, std::vector<int> &, const long, const int);

void distPar(const char *, const int, std::vector<int> &);
    
int getPartition(const char *, const int);

#endif
