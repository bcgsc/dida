#!/usr/bin/make -rRf

SHELL=/bin/bash

#------------------------------------------------------------
# test input/output files
#------------------------------------------------------------

# target seq for alignments
ref_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta
ref=ref.fa
test_ref=test_ref.fa

# query seqs for alignments
reads_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_1.fastq.gz
reads=reads.fq.gz
test_reads=test_reads.fq

# output alignment files
dida_mpi_sam=dida_mpi.sam
dida_wrapper_sam=dida_wrapper.sam
abyss_map_sam=abyss_map.sam

#------------------------------------------------------------
# params
#------------------------------------------------------------

# number of MPI tasks
np?=3
# number of threads per task
j?=1
# min align length
l?=60
# num of reads to align
n?=2500
# commands to run dida
dida_mpi_run?=mpirun -np `expr $(np) + 1` dida-mpi --se -j$j -l$l $(test_reads) $(test_ref)
dida_wrapper_run?=mpirun -np $(np) dida-wrapper --se -j$j -l$l $(wrapper_opt) $(test_reads) $(test_ref)

#------------------------------------------------------------
# special targets
#------------------------------------------------------------

.PHONY: clean identity_test

default: mpi_test wrapper_test wrapper_stream_test

clean:
	rm -f $(dida_mpi_sam) $(dida_wrapper_sam) \
		$(abyss_map_sam) $(test_reads) \
		ref-* *.lines

#------------------------------------------------------------
# downloading/building test input data
#------------------------------------------------------------

# download ref
$(ref):
	curl $(ref_url) > $@

# split ref into chunks of 100,000bp or less
$(test_ref): $(ref)
	fold -w 100000 $^ | awk '{print ">"i++; print $$0}' > $@

# download some reads
$(reads):
	curl $(reads_url) > $@

# extract first $n reads
$(test_reads): $(reads)
	zcat $(reads) | paste - - - - | head -$n | \
		tr '\t' '\n' > $@

#------------------------------------------------------------
# running DIDA/abyss-map
#------------------------------------------------------------

$(dida_mpi_sam):  $(test_reads) $(test_ref)
	$(dida_mpi_run) > $@

$(dida_wrapper_sam): $(test_reads) $(test_ref)
	$(dida_wrapper_run) > $@

$(abyss_map_sam): $(test_reads) $(test_ref)
	abyss-map --order -j$j -l$l $^ > $@

#------------------------------------------------------------
# tests
#------------------------------------------------------------

mpi_test: $(abyss_map_sam) $(dida_mpi_sam)
	compare-sam $(abyss_map_sam) $(dida_mpi_sam)
	@echo $@": PASSED!"

wrapper_test: $(abyss_map_sam) $(dida_wrapper_sam)
	compare-sam $(abyss_map_sam) $(dida_wrapper_sam)
	@echo $@": PASSED!"

wrapper_stream_test: $(test_reads) $(test_ref)
	mpirun -np $(np) bash -c \
		'dida-wrapper --se -j$j -l$l $(wrapper_opt) \
		<(abyss-tofastq $(test_reads)) $(test_ref) > /dev/null'
	@echo $@": PASSED!"

simple_identity_test: $(abyss_map_sam) $(dida_mpi_sam)
	comm --nocheck-order -3 \
		<(egrep -v '^@' $(abyss_map_sam) |awk '$$5 != 0') \
		<(egrep -v '^@' $(dida_sam) |awk '$$5 != 0')
