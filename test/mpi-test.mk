#!/usr/bin/make -rRf

SHELL=/bin/bash -o pipefail

# target seq for alignments
ref_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta
ref=ref.fa
# query seqs for alignments
reads_url:=http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_1.fastq.gz
reads=reads.fq
# output alignment file
aligns=aln.sam

# number of MPI tasks
np:=2
# number of target partitions 
p:=2
# number of threads per task
j:=1
# kmer size for bloom filters
b:=50

default: $(aligns)

# download ref and split into chunks of 100,000bp or less
$(ref):
	curl $(ref_url) | fold -w 100000 | awk '{print ">"i++; print $0}' > $(ref)

$(reads):
	curl $(reads_url) | gunzip -c | head -10000 > $(reads)

$(aligns):  $(reads) $(ref)
	mpirun -np $(np) dida-mpi -p$p -j$j -b$b $(reads) $(ref)
