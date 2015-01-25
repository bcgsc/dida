#!/usr/bin/make -rRf

SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# test input/output files
#------------------------------------------------------------

# target seq for alignments
ref_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta
ref=ref.fa
# query seqs for alignments
reads_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_1.fastq.gz
reads=reads.fq
reads_in=reads.in
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
n?=10000
# commands to run dida
dida_mpi_run?=mpirun -np $(np) dida-mpi --se -j$j -l$l $(reads) $(ref)
dida_wrapper_run?=mpirun -np $(np) dida-wrapper --se -j$j -l$l $(reads) $(ref)

#------------------------------------------------------------
# special targets
#------------------------------------------------------------

.PHONY: clean identity_test

default: mpi_test wrapper_test

clean:
	rm -f $(dida_mpi_sam) $(dida_wrapper_sam) $(abyss_map_sam) ref-* *.lines

#------------------------------------------------------------
# downloading/building test input data
#------------------------------------------------------------

# download ref and split into chunks of 100,000bp or less
$(ref):
	curl $(ref_url) | fold -w 100000 | awk '{print ">"i++; print $0}' > $(ref)

$(reads):
	curl $(reads_url) | gunzip -c | head -$n > $(reads)

#------------------------------------------------------------
# running DIDA/abyss-map
#------------------------------------------------------------

$(dida_mpi_sam):  $(reads) $(ref)
	$(dida_mpi_run) > $@

$(dida_wrapper_sam): $(reads_in) $(reads) $(ref)
	$(dida_wrapper_run) > $@

$(abyss_map_sam): $(reads) $(ref)
	abyss-map --order -l$l $(reads) $(ref) > $(abyss_map_sam)

#------------------------------------------------------------
# tests
#------------------------------------------------------------

mpi_test: $(abyss_map_sam) $(dida_mpi_sam)
	compare-sam $(abyss_map_sam) $(dida_mpi_sam)
	@echo $@": PASSED!"

wrapper_test: $(abyss_map_sam) $(dida_wrapper_sam)
	compare-sam $(abyss_map_sam) $(dida_wrapper_sam)
	@echo $@": PASSED!"

wrapper_stream_test: $(reads) $(ref)
	mpirun -np $(np) bash -c \
		'dida-wrapper --se -j$j -l$l \
		<(abyss-tofastq $(reads)) $(ref) > /dev/null'
	@echo $@": PASSED!"

simple_identity_test: $(abyss_map_sam) $(dida_mpi_sam)
	comm --nocheck-order -3 \
		<(egrep -v '^@' $(abyss_map_sam) |awk '$$5 != 0') \
		<(egrep -v '^@' $(dida_sam) |awk '$$5 != 0')
