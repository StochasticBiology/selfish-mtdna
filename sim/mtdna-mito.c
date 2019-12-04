// mtDNA selfishness simulation

// simulate time evolution of systems involving mixed mtDNA types which replicate and transcribe at different rates
// mitochondria are individuals; they survive or are recycled according to the complement of protein machinery that their mtDNAs have produced

// takes a command-line argument: 0 (span parameter space) or 1 (focus on a particular region of parameter space)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

// number of (stochastic) runs for each parameterisation
#define NIT 100

// time scale of simulation
#define MAXT 1000

int main(int argc, char *argv[])
{
  int m[1000], newm[1000];
  double p[1000], newp[1000];
  int n, newn;
  int t;
  int i, j;
  double lambda1, lambda2, gamma, kappa, kappap;
  double threshold;
  double *het, *tot, meanh, sdh, meant, sdt;
  int it;
  FILE *fp, *fp1;
  char str[200];
  double nowgamma;
  double l1start, l1end, l1step, l2step;
  int zoom, traceoutput;

  // process command-line arguments: if "1", we investigate a particular region of parameter space in detail, otherwise we span space
  zoom = 0;
  if(argc < 2)
    {
      printf("Zoom not specified: running default\n");
    }
  else if(atoi(argv[1]) == 1)
    {
      printf("Running zoomed version\n");
      zoom = 1;
    }
  else printf("Running default version\n");
  
  if(zoom)
    {
      l1start = 0.25; l1end = 0.3; l1step = 0.05;
      l2step = 0.01;
    }
  else
    {
      l1start = 0.1; l1end = 0.9; l1step = 0.1;
      l2step = 0.05;
    }

  // allocate memory to store statistics
  het = (double*)malloc(sizeof(double)*NIT*MAXT);
  tot = (double*)malloc(sizeof(double)*NIT*MAXT);

  // loop through values of lambda1
  for(lambda1 = l1start; lambda1 <= l1end; lambda1+= l1step)
    {
      lambda2 = 0.05; gamma = 0.1; kappa = 0.1; kappap = 0.01;
      threshold = 5;

      // open corresponding file for output 
      if(zoom)
	{
	  sprintf(str, "mtdna-mito-zoom-overall-%.3f.txt", lambda1);
	}
      else
	{
	  sprintf(str, "mtdna-mito-overall-%.3f.txt", lambda1);
	}
      fp = fopen(str, "w");

      // loop through values of lambda2
      for(lambda2 = 0.; lambda2 <= 1; lambda2 += l2step)
	{
	  // decide whether to output explicit time series or not (only for a couple of selected parameterisations)
	  if(zoom == 0 && lambda1 >= 0.29 && lambda1 <= 0.31 && ((lambda2 >= 0.19 && lambda2 <= 0.21) || (lambda2 >= 0.39 && lambda2 <= 0.41)))
	    traceoutput = 1;
	  else
	    traceoutput = 0;

	  if(traceoutput)
	    {
	      sprintf(str, "mtdna-mito-trace-%.3f-%.3f.txt", lambda1, lambda2);
	      fp1 = fopen(str, "w");
	    }

	  // loop through different values of the selective threshold
	  for(threshold = 0; threshold < 3; threshold += 0.1)
	    {
	      // loop through stochastic runs
	      for(it = 0; it < NIT; it++)
		{
		  printf("%f %f %i\n", lambda2, threshold, it);

		  // initial conditions: n is number of mitos, m[i] is the type of mtDNA in mito i; p[i] is the number of proteins in mito i
		  n = 20;
		  for(i = 0; i < n; i++)
		    {
		      m[i] = (i < 10 ? 1 : 2);
		      p[i] = 10;
		    }

		  // run through time course
		  for(t = 0; t < MAXT; t++)
		    {
		      // update gamma
		      nowgamma = gamma*(1.-n/1000.);

		      // buffer state of system
		      for(i = 0; i < n; i++)
			{
			  newm[i] = m[i];
			  newp[i] = p[i];
			}
		      newn = n;

		      // loop through mitos
		      for(i = 0; i < n; i++)
			{
			  if(RND < nowgamma && newn < 1000)
			    {
			      // this mito is going to do something
			      if((m[i] == 1 && RND < lambda1) || (m[i] == 2 && RND < lambda2))
				{
				  // this mito will replicate, producing a new one with same mtDNA type, and halving protein content between the two
				  newm[newn] = m[i];
				  newp[i] /= 2.;
				  newp[n] = p[i];
				  newn++;
				}
			      else
				{
				  // this mito will transcribe, incrementing its protein count
				  newp[i]++;
				}
			    }
			  if(RND < kappap && newp[i] > 0)
			    {
			      // this mito will lose a protein through degradation 
			      newp[i]--;
			    }
			}

		      // now loop through mitos applying selection
		      for(i = 0; i < newn; i++)
			{
			  if(RND < kappa)
			    {
			      // this mito is vulnerable
			      if(p[i] < threshold)
				{
				  // this mito is dysfunctional and will be degraded
				  // pop from mito list
				  for(j = i; j < newn-1; j++)
				    {
				      newm[j] = newm[j+1];
				      newp[j] = newp[j+1];
				    }
				  newn--;
				}
			    }
			}

		      // record statistics of system
		      het[MAXT*it + t] = 0;
		      for(i = 0; i < newn; i++)
			{
			  m[i] = newm[i];
			  p[i] = newp[i];
			  het[MAXT*it + t] += (m[i]-1);
			}
		      n = newn;
		      tot[MAXT*it + t] = n;
		      het[MAXT*it + t] /= newn;
		    }
		}
	      
	      // simulation is done: compute time series statistics of system
	      for(t = 0; t < MAXT; t++)
		{
		  meanh = sdh = 0; meant = sdt = 0;
		  for(it = 0; it < NIT; it++)
		    {
		      meanh += het[MAXT*it + t];
		      meant += tot[MAXT*it + t];
		    }
		  meanh /= NIT; meant /= NIT;
		  for(it = 0; it < NIT; it++)
		    {
		      sdh += (meanh-het[MAXT*it + t])*(meanh-het[MAXT*it + t]);
		      sdt += (meant-tot[MAXT*it + t])*(meant-tot[MAXT*it + t]);
		    }
		  sdh = sqrt(sdh / (NIT-1));
	  	  sdt = sqrt(sdt / (NIT-1));

		  if(traceoutput)
		    {
		      // output these time series to file
		      fprintf(fp1, "%.3f %i %f %f %f %f\n", threshold, t, meanh, sdh, meant, sdt);
		    }
		}
	      // output the mean genetic summary to file
	      fprintf(fp, "%.3f %.3f %f\n", lambda2, threshold, meanh);
	    }
	  fprintf(fp, "\n");
	  if(traceoutput)
	    {
	      fclose(fp1);
	    }
	}
      fclose(fp);
    }
  
  return 0;
}
