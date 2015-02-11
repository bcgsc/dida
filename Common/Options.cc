#include "config.h"
#include "Options.h"

namespace opt {
#if HAVE_LIBZ
	int gzip = 1;
#else
	int gzip = 0;
#endif
}
