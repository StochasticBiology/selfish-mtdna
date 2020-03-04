# wrapper script to compile and run mtDNA simulations: RTS with mito- and cell-level selection

gcc -o3 mtdna-mito.c -lm -o mtdna-mito.ce
gcc -o3 mtdna-cell.c -lm -o mtdna-cell.ce

./mtdna-mito.ce 0 > tmp0 &
./mtdna-mito.ce 1 > tmp1 &
./mtdna-cell.ce 0 > tmp2 &
