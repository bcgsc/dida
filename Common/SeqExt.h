#ifndef SEQEXT_H
#define SEQEXT_H 1

#include "Sequence.h" // for codeToBase
#include <cassert>
#include <ostream>
#include <stdint.h>

static const unsigned NUM_BASES = 4;

/** The adjacent vertices of a Kmer. */
class SeqExt
{
	public:
		SeqExt() : m_record(0) { };
		explicit SeqExt(uint8_t base) : m_record(1<<base) { };

		/** Return a SeqExt with the specified bits set. */
		static SeqExt mask(uint8_t bits)
		{
			SeqExt ext;
			ext.m_record = bits;
			return ext;
		}

		/** Return the out degree. */
		unsigned outDegree()
		{
			assert(m_record < 16);
			static const uint8_t popcount[16] = {
				0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4
			};
			return popcount[m_record];
		}

		/** Set the specified adjacency. */
		void setBase(uint8_t base)
		{
			m_record |= 1 << base;
		}

		/** Clear the specified adjacency. */
		void clearBase(uint8_t base)
		{
			m_record &= ~(1 << base);
		}

		/** Remove the specified edges. */
		void clear(SeqExt ext)
		{
			m_record &= ~ext.m_record;
		}

		/** Return wheter the specified base is adjacent. */
		bool checkBase(uint8_t base) const
		{
			return m_record & (1 << base);
		}

		/** Clear all adjacency. */
		void clear()
		{
			m_record = 0;
		}

		/** Return whether this kmer has any adjacent kmer. */
		bool hasExtension() const
		{
			return m_record > 0;
		}

		/** Return whether this kmer has more than one adjacent kmer.
		 */
		bool isAmbiguous() const
		{
			bool powerOfTwo = (m_record & (m_record - 1)) > 0;
			return m_record > 0 && powerOfTwo;
		}

		/** Return the complementary adjacency. */
		SeqExt complement() const;
		SeqExt operator ~() const { return complement(); }

		friend std::ostream& operator <<(std::ostream& out,
				const SeqExt& o)
		{
			assert(o.m_record < 1<<NUM_BASES);
			for (unsigned i = 0; i < NUM_BASES; ++i)
				if (o.checkBase(i))
					out << codeToBase(i);
			return out;
		}

	private:
		uint8_t m_record;
};

#endif
