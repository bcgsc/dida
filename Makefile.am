if HAVE_LIBZ
SUBDIRS = lib/gzstream-1.5
endif

bin_PROGRAMS = dida-mpi dida-wrapper prt dsp mrg

dida_mpi_CPPFLAGS = -I$(top_srcdir)/Common

dida_mpi_LDADD = $(MPI_LIBS)

dida_mpi_SOURCES = \
	Common/Uncompress.cpp \
	Common/Uncompress.h \
	Common/SignalHandler.cpp \
	Common/SignalHandler.h \
	Common/Fcontrol.cpp \
	Common/Fcontrol.h \
	Common/prt.cpp \
	Common/prt.h \
	Common/StringUtil.h \
	dida.cpp

dida_wrapper_CPPFLAGS = -I$(top_srcdir)/Common

if HAVE_LIBZ
dida_wrapper_CPPFLAGS += -I$(top_srcdir)/lib/gzstream-1.5
endif

dida_wrapper_LDADD = $(MPI_LIBS)

if HAVE_LIBZ
dida_wrapper_LDADD += $(top_builddir)/lib/gzstream-1.5/libgzstream.a
endif

dida_wrapper_SOURCES = \
	Common/Uncompress.cpp \
	Common/Uncompress.h \
	Common/SignalHandler.cpp \
	Common/SignalHandler.h \
	Common/Fcontrol.cpp \
	Common/Fcontrol.h \
	Common/GzipStream.h \
	Common/Options.h \
	Common/Options.cc \
	Common/prt.cpp \
	Common/prt.h \
	Common/mrg.cpp \
	Common/mrg.h \
	Common/StringUtil.h \
	dida_wrapper.cpp


prt_SOURCES = \
        prt.cpp

dsp_CPPFLAGS = -I$(top_srcdir)/Common

dsp_SOURCES = \
        Common/Uncompress.cpp \
        Common/Uncompress.h \
        Common/SignalHandler.cpp \
        Common/SignalHandler.h \
        Common/Fcontrol.cpp \
        Common/Fcontrol.h \
        Common/StringUtil.h \
        dsp.cpp

mrg_SOURCES = \
        mrg.cpp

