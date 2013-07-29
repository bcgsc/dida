#include "Sequence.h"
#include "Options.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib> // for abort
#include <iostream>
#include <sstream>

using namespace std;

enum { A, C, G, T };
static const int cstont[4][4] = {
	{ A, C, G, T },
	{ C, A, T, G },
	{ G, T, A, C },
	{ T, G, C, A }
};

/** Return the complement of the specified nucleotide. */
char complementBaseChar(char c)
{
	char rc;
	switch (toupper(c)) {
	  case 'A': rc = 'T'; break;
	  case 'C': rc = 'G'; break;
	  case 'G': rc = 'C'; break;
	  case 'T': rc = 'A'; break;
	  case 'N': rc = 'N'; break;
	  case '.': rc = '.'; break;
	  case 'M': rc = 'K'; break; // A or C
	  case 'R': rc = 'Y'; break; // A or G
	  case 'W': rc = 'W'; break; // A or T
	  case 'S': rc = 'S'; break; // C or G
	  case 'Y': rc = 'R'; break; // C or T
	  case 'K': rc = 'M'; break; // G or T
	  case 'V': rc = 'B'; break; // A or C or G
	  case 'H': rc = 'D'; break; // A or C or T
	  case 'D': rc = 'H'; break; // A or G or T
	  case 'B': rc = 'V'; break; // C or G or T
	  default:
		cerr << "error: unexpected character: `" << c << "'\n";
		assert(false);
		abort();
	}
	return islower(c) ? tolower(rc) : rc;
}

/** Return the reverse complement of the specified sequence. */
Sequence reverseComplement(const Sequence& s)
{
	Sequence rc(s);
	reverse(rc.begin(), rc.end());
	if (!opt::colourSpace)
		transform(rc.begin(), rc.end(), rc.begin(),
				complementBaseChar);
	return rc;
}

/** Return the base enumeration for the specified character. */
uint8_t baseToCode(char base)
{
	switch (base) {
		case 'A': case '0': return 0;
		case 'C': case '1': return 1;
		case 'G': case '2': return 2;
		case 'T': case '3': return 3;
	}
	cerr << "error: unexpected character: `" << base << "'\n";
	assert(false);
	abort();
}

char codeToBase(uint8_t code)
{
	assert(code < 4);
	return (opt::colourSpace ? "0123" : "ACGT")[code];
}

char colourToNucleotideSpace(char anchor, char cs)
{
	return cs == '.' ? 'N'
		: "ACGT"[cstont[baseToCode(anchor)][baseToCode(cs)]];
}

Sequence colourToNucleotideSpace(char anchor, const Sequence& seq)
{
	int seed = baseToCode(anchor);

	ostringstream s;
	s << anchor;
	for (string::const_iterator it = seq.begin();
			it != seq.end(); ++it) {
		seed = cstont[seed][baseToCode(*it)];
		s << codeToBase(seed);
	}
	return s.str();
}

char nucleotideToColourSpace(char a, char b)
{
	if (toupper(a) == 'N' || toupper(b) == 'N')
		return islower(a) || islower(b) ? 'n' : 'N';
	return "0123"[cstont[baseToCode(a)][baseToCode(b)]];
}

/** Convert the specified ambiguity code to a bitmask. */
unsigned ambiguityToBitmask(char c)
{
	if (isdigit(c)) // colour space
		return 1 << baseToCode(c);

	static const unsigned ambiguityToBitmaskTable[26] = {
		0x1, // 'A' ---A
		0xe, // 'B' TGC-
		0x2, // 'C' --C-
		0xd, // 'D' TG-A
		0x0, // 'E'
		0x0, // 'F'
		0x4, // 'G' -G--
		0xb, // 'H' T-CA
		0x0, // 'I'
		0x0, // 'J'
		0xc, // 'K' TG--
		0x0, // 'L'
		0x3, // 'M' --CA
		0xf, // 'N' ACGT
		0x0, // 'O'
		0x0, // 'P'
		0x0, // 'Q'
		0x5, // 'R' -G-A
		0x6, // 'S' -GC-
		0x8, // 'T' T---
		0x0, // 'U'
		0x7, // 'V' -GCA
		0x9, // 'W' T--A
		0x0, // 'X'
		0xa, // 'Y' T-C-
		0x0, // 'Z'
	};
	unsigned i = toupper(c) - 'A';
	assert(i < 26);
	unsigned x = ambiguityToBitmaskTable[i];
	assert(x > 0);
	return x;
}

/** Convert the specified bitmask to an ambiguity code. */
unsigned bitmaskToAmbiguity(unsigned x)
{
	static const char bitmaskToAmbiguityTable[16] = {
		'N', //----
		'A', //---A
		'C', //--C-
		'M', //--CA
		'G', //-G--
		'R', //-G-A
		'S', //-GC-
		'V', //-GCA
		'T', //T---
		'W', //T--A
		'Y', //T-C-
		'H', //T-CA
		'K', //TG--
		'D', //TG-A
		'B', //TGC-
		'N', //TGCA
	};
	assert(x < 16);
	return bitmaskToAmbiguityTable[x];
}
