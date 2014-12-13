#ifndef PRT_H
#define PRT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <time.h>

void qSort(std::vector<long> &, std::vector<int> &, int, int);

void tSort(std::vector<int> &, std::vector<int> &, int, int);

std::string proLine(std::string);

void getFname(const char *, std::string &, std::string &);

std::vector< std::vector<int> > getAdj(const char *, std::vector<int> &, long &, const int);

std::vector< std::vector<int> > wgetAdj(const char *, std::vector<int> &, long &);

std::vector< std::vector<int> > getCom(std::vector< std::vector<int> > &, std::vector<int> &);

int findBnode(long, std::vector<long> &, long);

void fullAloc(std::vector<int> &, int bNode, long, std::vector<int> &, std::vector<long> &);

void partAloc(std::vector<int> &, std::vector<long> &, std::vector<int> &, std::vector<int> &);

std::vector<int> compDist(std::vector< std::vector<int> > &, std::vector<int> &, const long, const int);

std::string getPrtFilename(const char *refName, const int procRank);

void debug_distPar(const char *, const int, std::vector<int> &, const int);

void ddistPar(const char *, const int, std::vector<int> &);

void wdistPar(const char *, const int, std::vector<int> &);

void debug_distTarget(const char *, const int,const int );

void ddistTarget(const char *, const int);

void wdistTarget(const char *, const int);

inline bool adjExist (const char *);

int debug_getPrt(const char *, const int, const int);

int dgetPrt(const char *, const int, const int);

int wgetPrt(const char *, const int);

#endif
