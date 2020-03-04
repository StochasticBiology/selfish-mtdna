// mtDNA selfishness simulation

// simulate time evolution of systems involving mixed mtDNA types which replicate and transcribe at different rates
// mitochondria are individuals; they survive or are recycled according to the complement of protein machinery that their mtDNAs have produced

// takes a command-line argument: 0 (span parameter space) or 1 (focus on a particular region of parameter space)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

// number of (stochastic) runs for each parameterisation
#define NIT 20

// time scale of simulation
#define MAXT 100
#define MAXM 100.

int main(int argc, char *argv[])
{
  int m1[1000], newm1[1000];
  int m2[1000], newm2[1000];
  double p[1000], newp[1000];
  int n, newn;
  int t;
  int i, j;
  double lambda1, lambda2, gamma, kappa, kappap;
  double threshold;
  double *het, *tot, *pdist, meanh, sdh, meant, sdt, meanp, sdp;
  int it;
  FILE *fp, *fp1;
  char str[200];
  double nowgamma;
  double l1start, l1end, l1step, l2step;
  int zoom, traceoutput;
  int r, killed;
  
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
  pdist = (double*)malloc(sizeof(double)*NIT*MAXT);

  // loop through values of lambda1
  for(lambda1 = l1start; lambda1 <= l1end; lambda1+= l1step)
    {
      lambda2 = 0.05; gamma = 0.1; kappa = 0.01; kappap = 0.1;
      threshold = 5;

      // open corresponding file for output 
      if(zoom)
	{
	  sprintf(str, "mtdna-mc-zoom-overall-%.3f.txt", lambda1);
	}
      else
	{
	  sprintf(str, "mtdna-mc-overall-%.3f.txt", lambda1);
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
	      sprintf(str, "mtdna-mc-trace-%.3f-%.3f.txt", lambda1, lambda2);
	      fp1 = fopen(str, "w");
	    }

	  // loop through different values of the selective threshold
	  for(threshold = 0*MAXM; threshold < 5*MAXM; threshold += 0.2*MAXM)
	    {
	      // loop through stochastic runs
	      for(it = 0; it < NIT; it++)
		{
		  printf("%f %f %i\n", lambda2, threshold, it);

		  // initial conditions: n is number of cells, m1/2[i] is the number of type 1/2 mtDNA in cell i; p[i] is the number of proteins in cell i
		  n = 100;
		  for(i = 0; i < n; i++)
		    {
		      m1[i] = RND*100;
		      m2[i] = 100-m1[i];
		      p[i] = (100*10);
		    }

		  // run through time course
		  for(t = 0; t < MAXT; t++)
		    {

		      // buffer state of system
		      for(i = 0; i < n; i++)
			{
			  newm1[i] = m1[i];
			  newm2[i] = m2[i];
			  newp[i] = p[i];
			}
		      newn = n;

		      // loop through cells
		      for(i = 0; i < n; i++)
			{
			  // update gamma
			  nowgamma = gamma*(1.-(m1[i]+m2[i])/MAXM);

			  // loop through mitos in cell
			  for(j = 0; j < m1[i]+m2[i]; j++)
			    {
			      if(RND < nowgamma && newn < 1000)
				{
				  // this mito is going to do something
				  if(j < m1[i])
				    {
				      if(RND < lambda1) newm1[i]++;
				      else newp[i]++;
				    }
				  else
				    {
				      if(RND < lambda2) newm2[i]++;
				      else newp[i]++;
				    }
				}
			      if(RND < kappa)
				{
				  if(j < m1[i]) newm1[i]--;
				  else newm2[i]--;
				}
			     
			      if(RND < kappap && newp[i] > 0)
				{
				  // this mito will lose a protein through degradation 
				  newp[i]--;
				}
			    }
			}
		      
		      // now loop through cells applying selection
		      killed = 0;
		      for(i = 0; i < newn; i++)
			{
			  //			  if(RND < kappa)
			  {
			    // this cell is vulnerable
			    if(newp[i] < threshold)
			      {
				// this cell is dysfunctional and will be degraded
				// pop from cell list
				for(j = i; j < newn-1; j++)
				  {
				    newm1[j] = newm1[j+1];
				    newm2[j] = newm2[j+1];
				    newp[j] = newp[j+1];
				  }
				newn--;
				i--;
				killed++;
			      }
			  }
			}

		      // if any cells survive, repopulate
		      if(newn > 0)
			{
		      for(i = 0; i < killed; i++)
			{
			  r = RND*newn;
			  newm1[newn+i] = newm1[r];
			  newm2[newn+i] = newm2[r];
			  newp[newn+i] = newp[r];
			}
		      newn += killed;
			}
		      
		      // record statistics of system
		      het[MAXT*it + t] = pdist[MAXT*it + t] = 0;
		      for(i = 0; i < newn; i++)
			{
			  m1[i] = newm1[i];
			  m2[i] = newm2[i];
			  p[i] = newp[i];
			  het[MAXT*it + t] += ((double)m2[i])/(m1[i]+m2[i]);
			  if(isnan(het[MAXT*it + t]))
			    printf("wtf\n");
			  pdist[MAXT*it + t] += p[i];
			}
		      n = newn;
		      tot[MAXT*it + t] = n;
		      het[MAXT*it + t] /= newn;
		    }
		}
	      
	      // simulation is done: compute time series statistics of system
	      for(t = 0; t < MAXT; t++)
		{
		  meanh = sdh = 0; meant = sdt = 0; meanp = sdp = 0;
		  for(it = 0; it < NIT; it++)
		    {
		      meanh += het[MAXT*it + t];
		      meant += tot[MAXT*it + t];
		      meanp += pdist[MAXT*it + t];
		    }
		  meanh /= NIT; meant /= NIT; meanp /= NIT;
		  for(it = 0; it < NIT; it++)
		    {
		      sdh += (meanh-het[MAXT*it + t])*(meanh-het[MAXT*it + t]);
		      sdt += (meant-tot[MAXT*it + t])*(meant-tot[MAXT*it + t]);
		      sdp += (meanp-pdist[MAXT*it + t])*(meanp-pdist[MAXT*it + t]);
		    }
		  sdh = sqrt(sdh / (NIT-1));
	  	  sdt = sqrt(sdt / (NIT-1));
		  sdp = sqrt(sdp / (NIT-1));

		  if(traceoutput)
		    {
		      // output these time series to file
		      fprintf(fp1, "%.3f %i %f %f %f %f %f %f\n", threshold, t, meanh, sdh, meant, sdt, meanp, sdp);
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
