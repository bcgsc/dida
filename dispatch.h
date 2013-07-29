#ifndef DISPATCH_H
#define DISPATCH_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "BloomFilter.h"
#include "HashManager.h"
#ifdef _OPENMP
# include <omp.h>
#endif

#define BIT_PER_ITEM 8
#define HASH_NUMBER 5

void getFnameb(const char *, std::string &, std::string &);

long getInfo2(const char *, int);

std::string getMin(const std::string&);

int getDispatch(const char *, const int, const int, const char *);

#endif
