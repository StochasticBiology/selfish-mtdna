# wrapper script to compile and run bioinformatics pipeline: G-quad analysis code and subsequent statistical and summary analysis

gcc -o3 find-gquads.c -lm -o find-gquads.ce
chmod +x process-ncbi.sh
chmod +x parse-kang.sh

./find-gquads.ce kang-fasta.fasta
./parse-kang.sh
R CMD BATCH kang-competition.R

./find-gquads.ce mtdna-sequences.fasta
./process-ncbi.sh mtdna-sequences.fasta-out.txt
gcc make-ps.c -lm -o make-ps.ce
./make-ps.ce 0
./make-ps.ce 1


