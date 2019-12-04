# script to process output of find-gquads.c for summary and visualisation
# takes a command-line argument: the name of the file to analyse

# separate out accession ref and G-quad statistics from output of find-gquads.c

awk 'BEGIN{FS="|";}{print $4, $5;}' $1 | awk '{print $1, $(NF-7), $(NF-3), $(NF-6), $(NF-5), $(NF-4), $(NF-2), $(NF-1), $NF;}' | sed 's/[.][0-9]//g' | sort > summary-$1

# cross-reference this output with database.txt to output haplogroups and haplotypes

awk 'NR==FNR{ a[$3]=$0; next }
     {  for(i in a) { 
            if (a[$1]) { 
                print a[$1], $0; delete a[i]; break 
            } 
        }
     }' database.txt summary-$1 > hist-scores-final.txt
awk 'BEGIN { convert="ABCDEFGHIJKLMNOPQRSTUVWXYZ" } 
       { num=index(convert,substr($0,1,1))+64; print num, $0; }' hist-scores-final.txt > tmp

# output rearrangements of these for PostScript code and reference

awk '{print $1, $2, $3, $4, $5, $8, $9, $10, $11, $12, $13;}' tmp > hist-plot-final.txt
awk '{print $4, $2, $3, $6, $7, $8, $9, $10, $11, $12, $13;}' tmp > human-set-final.txt
