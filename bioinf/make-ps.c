// code to produce a PostScript-format graphical summary (motif-plot-style) of G-quadruplex region statistics from previous analysis

// takes command-line argument: which region to examine: 0 CSB2, 1 TAS

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  FILE *fp;
  double finds[26][100];
  char ch;
  int ref;
  int r, a, b, c, tmp;
  int i, j;
  double total, sum;
  int region;
  char fstr[100];
  int skips;
  
  // process command-line argument
  region = 0;
  if(argc < 2)
    {
      printf("Region not specified, looking at CSB2\n");
    }
  else if(atoi(argv[1]) == 1)
    {
      printf("Looking at TAS\n");
      region = 1;
    }
  else printf("Looking at CSB2\n");
  
  for(i = 0; i < 26; i++)
    {
      for(j = 0; j < 100; j++)
	finds[i][j] = 0;
    }

  // very lazy code to read output of previous pipeline and extract G-quadruplex structural information
  fp = fopen("hist-plot-final.txt", "r");
  do{
    fscanf(fp, "%i", &tmp);
    if(feof(fp)) break;
    r = tmp-65;
    if(region == 1) skips = 8;
    else skips = 5;

    for(i = 0; i < skips; i++)
      {
	do{ ch = fgetc(fp); } while(ch != ' ');
      }
    fscanf(fp, "%i %i %i", &a, &b, &c);
    do{ ch = fgetc(fp); } while(ch != '\n');
    if(b > 0) ref = a*10+c;
    else ref = a;
    finds[r][ref]++;
    printf("%i %i\n", r, ref);
  }while(!feof(fp));
  fclose(fp);

  // open output file and output preamble
  sprintf(fstr, "motif-plot-%i.ps", region);
  fp = fopen(fstr, "w");
  fprintf(fp, "%%!PS-Adobe-2.0\n/Times-Roman findfont\n8 scalefont\nsetfont\nnewpath\n");

  // loop through haplogroups
  for(i = 0; i < 26; i++)
    {
      // sum up different structures
      total = 0;
      for(j = 0; j < 100; j++)
	total += finds[i][j];
      sum = 0;

      // set up output for this haplogroup
      fprintf(fp, "%i %i moveto\n", 50+(i%13)*10+3, 100-4-(i >= 13 ? 18 : 0));
      fprintf(fp, "0 0 0 setrgbcolor\n");
      fprintf(fp, "[0.4 0 0 0.4 0 0] concat\n");
      fprintf(fp, "(%c) show\n", (char)(i+65));
      fprintf(fp, "initmatrix\n");
      fprintf(fp, "%i %i moveto\n", 50+(i%13)*10+3, 100-6-(i >= 13 ? 18 : 0));
      fprintf(fp, "0 0 0 setrgbcolor\n");
      fprintf(fp, "[0.25 0 0 0.25 0 0] concat\n");
      fprintf(fp, "([%i]) show\n", (int)total);
      fprintf(fp, "initmatrix\n");


      // loop through the set of possible found structures, outputting a motif for each
      for(j = 0; j < 100; j++)
	{
	  if(finds[i][j] > 0)
	    {
	      fprintf(fp, "%i %.3f moveto\n", 50+(i%13)*10, 100+sum*8-(i >= 13 ? 18 : 0));
	      fprintf(fp, "[1 0 0 %.3f 0 0] concat\n", 1.5*finds[i][j]/total);
	      switch(j)
		{
  		case 0: fprintf(fp, "0.7 0.7 0.7 setrgbcolor\n"); break;
		case 44: fprintf(fp, "0.8 0.7 0.5 setrgbcolor\n"); break;
		case 34: fprintf(fp, "0 1 1 setrgbcolor\n"); break;
		case 86: fprintf(fp, "0 0.4 0 setrgbcolor\n"); break;
		case 76: fprintf(fp, "1 0 1 setrgbcolor\n"); break;
		case 96: fprintf(fp, "0 0 1 setrgbcolor\n"); break;
		default: fprintf(fp, "0 0 0 setrgbcolor\n"); break;
		}
	      if(j > 20) fprintf(fp, "(%i/%i) show\n", j/10, j-10*(j/10));
	      else fprintf(fp, "(%i/) show\n", j);
	      //fprintf(fp, "newpath\n[1 0 0 1 0 0] concat\n", 10*finds[i][j]/total);
	      fprintf(fp, "initmatrix\n");
	      sum += finds[i][j]/total;
	    }
	}
    }
  fclose(fp);

  return 0;
}
