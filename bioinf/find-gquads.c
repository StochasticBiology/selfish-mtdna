// code to find and score potential G-quadruplex-forming sequences at given regions in all entries is a given FASTA file
// takes a filename from the command line: the FASTA file to analyse

#include <stdio.h>

int main(int argc, char *argv[])
{
  FILE *fp, *fpo, *fpm;
  char str[100000];
  int n;
  int i, j, k;
  char ch;
  int maxk;
  int nfound;
  int seti[1000], setk[1000];
  int c1count, gapcount, c2count, end;
  int csb2score;
  int region;

  // number of regions to analyse, and start loci for each 
  int nregions = 3;
  int regionstarts[] = {300, 16175, 16345};

  // process command-line argument
  if(argc < 2)
    {
      printf("Need file to analyse\n");
      return 0;
    }
  fp = fopen(argv[1], "r");
  if(fp == NULL)
    {
      printf("Couldn't find file %s\n", argv[1]);
      return 0;
    }

  // open files for output
  sprintf(str, "%s-out.txt", argv[1]);
  fpo = fopen(str, "w");
  sprintf(str, "%s-tas-seqs.txt", argv[1]);
  fpm = fopen(str, "w");

  // read through FASTA file. output label for each entry, and store the sequence in memory (str)
  do{ch = fgetc(fp);}while(ch != '>');
  do{
    if(ch == '>')
      {
	do{ch = fgetc(fp); if(ch != '\n') fprintf(fpo, "%c", ch);}while(ch != '\n');
      }
    fprintf(fpo, " ");
    n = 0;
    do{
      ch = fgetc(fp);
      if(ch == '>') break;
      if(ch != '\n') str[n++] = ch;
      if(feof(fp)) break;
    }while(ch != '>');
    str[n] = '\0';

    // loop through the number of specified regions to seek G-quad sequences
    for(region = 0; region < nregions; region++)
      {
	end = 0;
	// start at given region start
	for(i = regionstarts[region]; end == 0; i++)
	  {
	    // first we seek 3 Cs in a row: record the total number of Cs in this tract (c1count)
	    // then record the length of the following tract until another 3 Cs are reached (gapcount)
	    // then record the total number of Cs in the following tract (c2count)
	    c1count = gapcount = c2count = end = 0;
	    if(str[i] == 'C' && str[i+1] == 'C' && str[i+2] == 'C')
	      {
		for(j = i; end == 0; j++)
		  {
		    if(gapcount == 0 && str[j] == 'C') c1count++;
		    else if(c2count == 0 && !(str[j] == 'C' && str[j+1] == 'C' && str[j+2] == 'C')) gapcount++;
		    else if(end == 0 && str[j] == 'C') c2count++;
		    else if(c2count != 0 && str[j] != 'C') end++;
		  }

		// now examine the lengths of these tracts and output accordingly
		if(gapcount > 7)
		  {
		    // if we've got a gap of more than 7
		    if(c1count > 6)
		      {
			// the first tract is long enough to form a G-quad itself
			for(j = 0; j < c1count; j++)
			  fprintf(fpo, "%c", str[i+j]);
			gapcount = 0; c2count = 0;
		      }
		    else
		      {
			// no G-quad candidate here
			fprintf(fpo, "-");
			c1count = gapcount = c2count = 0;
		      }
		  }
		else
		  {
		    // we have a more "normal" G-quad structure
		    for(j = 0; j < c1count+gapcount+c2count; j++)
		      fprintf(fpo, "%c", str[i+j]);
		  }
	      }
	  }
	fprintf(fpo, " %i %i %i ", c1count, gapcount, c2count); 
      }
    
    fprintf(fpo, "\n");

    // for debugging and further analysis, output the 30 bases after each region start to a separate file
    for(region = 0; region < nregions; region++)
      {
	for(i = regionstarts[region]; i < regionstarts[region]+30; i++)
	  fprintf(fpm, "%c", str[i]);
	fprintf(fpm, " ");
      }
    fprintf(fpm, "\n");
  }while(!feof(fp));

  fclose(fp);
  fclose(fpo);
  fclose(fpm);
  
  return 0;
}
