1. Compiling DIDA from source
To compile and install DIDA in /usr/local:

	./configure
	make 
	sudo make install 

To install DIDA in a specified directory:

	./configure --prefix=/opt/dida
	make 
	sudo make install 

DIDA uses OpenMP for parallelization, which requires a modern compiler such as GCC 4.2 or greater. If you have an older compiler, it is best to upgrade your compiler if possible. If you have multiple versions of GCC installed, you can specify a different compiler:

	./configure CC=gcc-4.9 CXX=g++-4.9 

If you wish to build the DIDA with MPI support, MPI should be found in /usr/include and /usr/lib or its location specified to configure:

	./configure --with-mpi=/usr/lib/openmpi 

To run DIDA, its executables should be found in your PATH. If you installed DIDA in /opt/DIDA, add /opt/DIDA/bin to your PATH:

	PATH=/opt/dida/bin:$PATH

To run DIDA, its executables should be found in your PATH. If you installed DIDA in /opt/DIDA, add /opt/DIDA/bin to your PATH:

	PATH=/opt/dida/bin:$PATH



2. Distribute Indexing and Alignment 

To align a library of paired reads named <query> against the target named <target> using the batch version of DIDA, run the set of following commands:

	prt –p <partition> <target>
	index <sub-target>
	dsp –p <partition> -b <bmer> <query>
	map <sub-query> <sub-target> 
	mrg –p <partition> -a <aligner> -m <mode>

Example of aligning a subset of chromosome 14 reads, chr14.in, against its draft assembly, CHR14.fa,CHR14.adj on 4 nodes using DIDA+ABySS-map. To install ABySS-map from ABySS Package:	http://www.bcgsc.ca/platform/bioinfo/software/abyss 

	1. prt –p4 CHR14.fa
	2. 
	for i in {1..4}
	do
			abyss-index mref-$i.fa
	done
	3. dsp -p4 -b28 chr14.in
	4. 
	for i in {1..4}
	do
			abyss-map -l28 –j12 --order mreads-$i.fa mref-$i.fa > aln-$i.sam
	done
	5. mrg -p4 -a abyss-map -m ord

If we want to use bwa-mem instead of ABySS-map steps 2, 4 and 5 will be:

	2. 
	for i in {1..4}
	do
			bwa index mref-$i.fa
	done

	4. 
	for i in {1..4}
	do
			bwa mem –t12 –k28 mref-$i.fa mreads-$i.fa > aln-$i.sam
	done

	5. mrg -p4 -a bwa-mem -m ord




3. Running DIDA on a cluster

DIDA integrates well with cluster job schedulers, such as:
	SGE (Sun Grid Engine)
	Portable Batch System (PBS)
	Load Sharing Facility (LSF)
	IBM LoadLeveler

For example, to submit the above DIDA job for aligning a subset of chromosome 14 reads, chr14.in, against its draft assembly, CHR14.fa,CHR14.adj on 4 nodes using DIDA+ABySS-map:

	qsub -N dida-chr14 -pe openmpi 5 <<< 'mpirun –np 5 dida-wrapper –b28 –l28 –a abyss-map –m ord chr14.in CHR14.fa'

DIDA has a fully streamlined version, which performs the dispatch, align, and merge steps over MPI communications. To run that version, run the command:

	qsub -N dida-chr14 -pe openmpi 6 <<< 'mpirun –np 5 dida-mpi –l28 frag_1.fastq frag_2.fastq CHR14.fa > res-chr14.sam'

