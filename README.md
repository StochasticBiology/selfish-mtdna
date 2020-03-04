# selfish-mtdna

Code for simulation and analysis of selfish mtDNA dynamics, and links with G-quadruplexes in control region

C, R, shell script, and Gnuplot scripts; one preprocessing and one visualisation Mathematica notebook are included but not required to proceed.

sim/ -- RTS simulation of selfish mtDNA dynamics. these simulations should take a few hours to run on a modern machine.
  mtdna-mito.c -- simulates RTS model with mito-level selection. command-line argument determines whether parameter space is spanned (0) or a zoomed-in region is investigated (1)
  mtdna-cell.c -- simulates RTS model with cell-level selection (SI of article), command-line argument determines structure of model (0, rates constant; 1, rates scaled by copy number) 

  plot-rts.sh -- plots output of mtdna-mito.c (1) in Gnuplot (used in Fig 2 of article)
  plot-si-rts-trellis.sh -- plots output of mtdna-mito.c (0) in Gnuplot (Fig S1)
  plot-si-cell-trellis.sh -- plots output of mtdna-cell.c in Gnuplot (Fig S3)
  plot-si-rts-time-series.sh -- plots time series output of mtdna-mito.c (0) in Gnuplot (Fig S2)
  
bioinf/ -- analyses G-quad regions in sequences from NCBI and Kang et al. papers. this analysis should take a few minutes on a modern machine.
  rcrs.txt -- data file containing revised Cambridge reference sequence (RCRS)
  edits.ods -- data file containing differences between Kang et al. sequences and RCRC. Taken directly from the SI of Kang et al. 2019.
  database.txt -- data file containing haplogroups and haplotypes for a set of NCBI accessions
  (reconstruct-sequences.nb) -- Mathematica file using edits.ods and rcrs.txt to reconstruct sequences from Kang et al. paper into FASTA format
  kang-fasta.fasta -- data file, output from the above
  (mtdna-sequences.fasta) -- data file containing mtDNA sequences in FASTA format from NCBI. NOTE: this file is included to illustrate the analysis. the true data file containing all mtDNA sequences analysed in this project is too big (~0.5GB) to include in this repository, but can be downloaded from NCBI through batch downloading the accessions in human-set-article.txt
  
  find-gquads.c -- identifies and scores potential G-quadruplex sequences at particular regions in all sequences in a given FASTA file. run this on kang-fasta.fasta and mtdna-sequences.fasta to provide output for the following analyses.

  parse-kang.sh -- rearranges output of find-gquads.c applied to Kang et al. data for processing in reverting/non-reverting pairs
  kang-competition.R -- reads output of the above, compares logistic regression models, and calculates predictive power of simple model
  plot-kang-data.sh -- plots inferred relationship and mtDNA statistics (used in Fig 4 of article)

  process-ncbi.sh -- rearranges and organises output of find-gquads.c for processing
  make-ps.c -- produces PostScript figure from these data (Fig S4)

  (human-set-article.txt) -- summary output of above pipeline run on full mtDNA sequence set, containing accessions for all sequences used in full analysis
  
alignments/ -- Clustal Omega alignments of control regions from Kang et al. and mtDNA from mouse studies
  kang-alignment.clustal_num -- Kang et al. sequences
  mouse-cr-alignment.clustal_num -- Mouse sequences

summary/ -- summary and fit statistics of tissue-specific mouse segregation observations 
  mouse-segregation.ods -- spreadsheet of tissue-specific mouse segregation observations and fit statistics
  (present-mouse-tissues.nb) -- Mathematica notebook producing a graphical visualisation of these observations (used in Fig 2 of article)


