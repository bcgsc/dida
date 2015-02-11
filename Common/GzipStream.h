#ifndef _ZIPSTREAM_H_
#define _ZIPSTREAM_H_

#include "Options.h"
#include <gzstream.h>
#include <cassert>

static inline std::istream* openInputStream(std::string path)
{
	std::istream* s;
#if HAVE_LIBZ
	if (opt::gzip)
		s = new igzstream(path.c_str());
	else
		s = new std::ifstream(path.c_str());
#else
	s = new std::ifstream(path.c_str());
#endif
	assert(*s);
	return s;
}

static inline void closeInputStream(std::istream* s)
{
	assert(s != NULL);
#if HAVE_LIBZ
	if (opt::gzip)
		((igzstream*)s)->close();
	else
		((std::ifstream*)s)->close();
#else
	((std::ifstream*)s)->close();
#endif
	delete s;
}

static inline std::ostream* openOutputStream(std::string path)
{
	std::ostream* s;
#if HAVE_LIBZ
	if (opt::gzip)
		s = new ogzstream(path.c_str());
	else
		s = new std::ofstream(path.c_str());
#else
	s = new std::ofstream(path.c_str());
#endif
	assert(*s);
	return s;
}

static inline void closeOutputStream(std::ostream* s)
{
	assert(s != NULL);
#if HAVE_LIBZ
	if (opt::gzip)
		((ogzstream*)s)->close();
	else
		((std::ofstream*)s)->close();
#else
	((std::ofstream*)s)->close();
#endif
	delete s;
}

#endif
