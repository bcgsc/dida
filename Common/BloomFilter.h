/**
 * A Bloom filter
 * Copyright 2013 Shaun Jackman
 *
 * Modifications by Ben Vandervalk.
 */
#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H 1

#include "Common/IOUtil.h"
#include "Common/BitUtil.h"
#include "Common/MurmurHash2.h"
#include "Common/StringUtil.h"

#include <string>
#include <algorithm>
#include <iostream>
#include <map>

#define ARRAY_SIZE(a) (sizeof(a)/sizeof(a[0]))

/** Version number for Bloom filter file format. */
#define FILE_FORMAT_VERSION (unsigned)0

/* Pre-defined headers in Bloom filter files. */
#define VERSION_HEADER "file-format-version"
#define NUM_HASHES_HEADER "num-hashes"
#define BLOOM_SIZE_HEADER "bloom-size"

/** Marks end of header section in Bloom filter file. */
#define END_HEADER_MARKER "end-header"
	
/** A Bloom filter. */
class BloomFilter
{
  public:
	
	/** A set of headers in a Bloom filter file. */
	typedef std::map<std::string, std::string> FileHeaderMap;
	/** A header line in a Bloom filter file. */
	typedef FileHeaderMap::value_type FileHeader;
	
	/** Default constructor. */
	BloomFilter() : m_size(0), m_numHashes(0), m_array(NULL) { }

	/** Constructor. */
	BloomFilter(size_t n, unsigned numHashes=1) :
		m_size(n), m_numHashes(numHashes)
	{
		m_array = new char[(n + 7)/8]();
	}

	/** Destructor. */
	~BloomFilter()
	{
		delete[] m_array;
	}

	/** Return the size of the bit array. */
	size_t size() const { return m_size; }

	/** Return the population count, i.e. the number of set bits. */
	size_t popcount() const
	{
		return ::popcount(m_array, m_size);
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return (double)popcount() / size();
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		assert(i < m_size);
		return m_array[i / 8] & 1 << (7 - i % 8);
	}

	/** Return whether the object is present in this set. */
	bool operator[](const char* key) const
	{
		size_t len = strlen(key);
		for (unsigned i = 0; i < m_numHashes; ++i) {
			if (!(*this)[MurmurHash64A(key, len, i) % m_size])
				return false;
		}
		return true;
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t i)
	{
		assert(i < m_size);
		m_array[i / 8] |= 1 << (7 - i % 8);
	}

	/** Add the object to this set. */
	void insert(const char* key)
	{
		size_t len = strlen(key);
		for (unsigned i = 0; i < m_numHashes; ++i)
			insert(MurmurHash64A(key, len, i));
	}

	/** Operator for reading a bloom filter from a stream. */
	friend std::istream& operator>>(std::istream& in, BloomFilter& o)
	{
		o.read(in);
		return in;
	}

	/** Operator for writing the bloom filter to a stream. */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilter& o)
	{
		o.write(out);
		return out;
	}

	/** Read a bloom filter from a stream. */
	void read(std::istream& in)
	{
		FileHeaderMap headers = readFileHeader(in);
		assert(in);

		// check that all required headers are present

		const char* REQUIRED_HEADERS[] = {
			VERSION_HEADER,
			NUM_HASHES_HEADER,
			BLOOM_SIZE_HEADER
		};

		for (unsigned i = 0; i < ARRAY_SIZE(REQUIRED_HEADERS); ++i) {
			if (headers.find(REQUIRED_HEADERS[i]) == headers.end()) {
				std::cerr << "Bloom filter file is missing required "
					"header '" << REQUIRED_HEADERS[i] << "'\n";
				exit(EXIT_FAILURE);
			}
		}

		// check for compatible file format version
		
		std::stringstream ss;
		unsigned version;
		FileHeaderMap::iterator i = headers.find(VERSION_HEADER);
		assert(i != headers.end());
		stringToVal(i->second, version);

		if (version != FILE_FORMAT_VERSION) {
			std::cerr << "error: Bloom filter file format (`"
				<< version << "'), does not match version required "
				"by this program (`" << FILE_FORMAT_VERSION << "').\n";
			exit(EXIT_FAILURE);
		}

		// set number of hash functions from file

		i = headers.find(NUM_HASHES_HEADER);
		assert(i != headers.end());
		stringToVal(i->second, m_numHashes);
		
		// set Bloom size from file

		size_t size;
		i = headers.find(BLOOM_SIZE_HEADER);
		assert(i != headers.end());
		stringToVal(i->second, size);
		if (m_size != size)
			resize(size);
		
		// read in Bloom filter bits

		size_t numBytes = (size + 7) / 8;
		in.read(m_array, numBytes);
		assert(in);
	}

	/** Write a bloom filter to a stream. */
	void write(std::ostream& out) const
	{
		FileHeaderMap headers;
		headers[VERSION_HEADER] = toString(FILE_FORMAT_VERSION);
		headers[NUM_HASHES_HEADER] = toString(m_numHashes);
		headers[BLOOM_SIZE_HEADER] = toString(m_size);

		writeFileHeader(out, headers);
		assert(out);

		out.write(m_array, (m_size + 7)/8);
		assert(out);
	}

	/** Read header section of Bloom filter file. */
	static FileHeaderMap readFileHeader(std::istream& in)
	{
		FileHeaderMap headers;
		while(in) {
			std::string headerLabel;
			std::string headerValue;
			in >> headerLabel;
			assert(in);
			if (headerLabel == END_HEADER_MARKER) {
				in >> expect("\n");
				break;
			}
			in >> expect("\t") >> headerValue >> expect("\n");
			assert(in);
			headers[headerLabel] = headerValue;
		}	
		assert(in);
		return headers;
	}

	/** Resize the bloom filter (wipes the current data) */
	void resize(size_t size)
	{
		if (m_size > 0 && m_array != NULL)
			delete[] m_array;

		m_array = new char[(size + 7)/8]();
		m_size = size;
	}

  protected:

	/** Write header section of Bloom filter file. */
	static void writeFileHeader(std::ostream& out,
		 const FileHeaderMap& headers)
	{
		for (FileHeaderMap::const_iterator i = headers.begin();
			i != headers.end(); ++i) {
			out << i->first << "\t" << i->second << "\n";
		}
		out << END_HEADER_MARKER << "\n";
		assert(out);
	}

	/** Size of Bloom filter in bits. */
	size_t m_size;
	/** Number of hash functions used by Bloom filter. */
	unsigned m_numHashes;
	/** Bloom filter bit array. */
	char* m_array;
};

#endif
