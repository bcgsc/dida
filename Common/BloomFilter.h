/*
 * BloomFilter.h
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include "HashManager.h"
#include <string>

using namespace std;

static const size_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
		0x40, 0x80 };

class BloomFilter {
public:
	//for generating a new filter
	BloomFilter(size_t filterSize, HashManager const &hashFns);
	void insert(string const &kmer);

	//for storing/restoring the filter
	void storeFilter(string const &filterFilePath) const;
	BloomFilter(size_t filterSize, string const &filterFilePath, HashManager const &hashFns);

	const bool contains(string const &kmer);
	const bool contains(vector<size_t> const &precomputed);

	virtual ~BloomFilter();
private:
	void initSize(size_t size);
	char* filter;
	size_t size;
	size_t sizeInBytes;
	HashManager multiHasher;
};

#endif /* BLOOMFILTER_H_ */
