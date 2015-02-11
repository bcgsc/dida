#ifndef MRG_H
#define MRG_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <cstdlib>
#include <stdint.h>
#include <getopt.h>


struct samHed
{
    std::string SQ1;
	char sn1;
    char sn2;
    char sn3;
	unsigned sqId;
    std::string SQ3;
    int hedPr;
};

struct samRec
{
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

std::string getSamFilename(int);

std::string getUnmappedSamFilename();

void getInf(unsigned &, unsigned &);

int memory_usage();

bool operator>(const samHed&, const samHed&);

samHed hedLoad(std::string& , int );

std::ostream& operator<<(std::ostream& , const samHed&);

bool operator>(const samRec& , const samRec&);

samRec recLoad(std::string& , int);

std::ostream& operator<<(std::ostream& , const samRec&);

void memMer(const int, const std::string &);

void fstMer(const int, const std::string &);

void fordMer(const int, const std::string &);

void bestMer(const int, const std::string &);;

int call_merger(const int, const std::string &, const std::string &, const unsigned);

#endif

