# wrapper script to compile and run bioinformatics pipeline: G-quad analysis code and subsequent statistical and summary analysis

gcc -o3 find-gquads.c -lm -o find-gquads.ce

./find-gquads.ce kang-fasta.fasta
R CMD BATCH kang-competition.R

./find-gquads.ce mtdna-sequences.fasta
./process-ncbi.sh mtdna-sequences.fasta-out.txt
gcc make-ps.c -lm -o make-ps.ce
./make-ps.ce
