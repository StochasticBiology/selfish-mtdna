# wrapper script to compile and run mtDNA simulations: RTS with mito- and cell-level selection

gcc -o3 mtdna-mito.c -lm -o mtdna-mito.ce
gcc -o3 mtdna-cell.c -lm -o mtdna-cell.ce

# default mito-level 
./mtdna-mito.ce 0 0 0 > tmp000 &

# default mito-level zoomed in
./mtdna-mito.ce 1 0 0 > tmp100 &

# different protein mixing
./mtdna-mito.ce 0 0 0.01 > tmp001 &
./mtdna-mito.ce 0 0 0.1 > tmp002 &

# different replication failure rates
./mtdna-mito.ce 0 0.2 0 > tmp010 &
./mtdna-mito.ce 0 0.5 0 > tmp020 &

# cell-level
./mtdna-cell.ce 0 > tmp2 &
