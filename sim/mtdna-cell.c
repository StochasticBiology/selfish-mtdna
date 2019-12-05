// mtDNA selfishness simulation: RTS model with cell-level selection

// simulate time evolution of cellular systems involving mixed mtDNA types which replicate and transcribe at different rates
// cells are individuals; after simulating a whole bunch of initial conditions we look at the population to see what proportion ends up exceeding a "bioenergetic threshold" defined by mtDNA and protein numbers
// the "tissue-wide" statistics for different rate parameters are then output

// takes a command-line argument: whether to scale rates by current state (1) or not (0)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NSAMP 100000

#define BASIC 0

// maximum number of "marked" initial conditions that we can record for successfull initial condition states
#define MAXMARK 10000

#define RND drand48()

// structure to hold parameterisation for the model
typedef struct tmpParams 
{
  // c reflects the overall "activity" of mtDNA
  // lambdaw, lambdam describe how this activity is partitioned into replication (lambda) vs transcription (1-lambda)
  // nu, nup are degradation rates of mtDNA and protein respectively
  // t, model are timescale and structure of model respectively
  
  double c, nu, nup, lambdaw, lambdam;
  double t;
  int model;
} Params;

// simulate or solve the ODEs governing the system

// model 0: simple ODEs, separate rates for each species (explicit solution possible)
// model 1: rates normalised by total mtDNA content (i.e. cellular resource required for replication is not infinite but must be shared across all mtDNAs) (numerical solution required)

// simulate or solve the ODEs governing the system
// w0, m0 are mtDNA initial conditions
// (0, 0) are protein initial conditions
// w, m are mtDNA levels with time
// pw, pm are protein levels with time
// P.t is the timescale of the simulation

void Simulate(Params P, int w0, int m0, double *w, double *m, double *pw, double *pm)
{
  double t, dt = 0.01;
  double dw, dm, dpw, dpm;

  if(P.model == 0)
    {
      // explicit solutions of model 0 ODEs
      *w = w0*exp(P.t*(P.c*P.lambdaw-P.nu));
      *m = m0*exp(P.t*(P.c*P.lambdam-P.nu));
      *pw = -(P.c*exp(-P.nup*P.t)*(-1.+exp(P.t*(P.c*P.lambdaw-P.nu+P.nup)))*w0*(-1+P.lambdaw))/(P.c*P.lambdaw-P.nu+P.nup);
      *pm = -(P.c*exp(-P.nup*P.t)*(-1.+exp(P.t*(P.c*P.lambdam-P.nu+P.nup)))*m0*(-1+P.lambdam))/(P.c*P.lambdam-P.nu+P.nup);
    }
  if(P.model == 1)
    {
      // initial conditions
      *w = w0; *m = m0; *pw = 0; *pm = 0;
      // simple euler solver
      for(t = 0; t < P.t; t += dt)
	{
	  // here the replication rates are normalised by total mtDNA content
	  dw = -P.nu*(*w) + P.c*P.lambdaw*(*w)/((*w)+(*m));
	  dm = -P.nu*(*m) + P.c*P.lambdam*(*m)/((*w)+(*m));
	  dpw = -P.nup*(*pw) + P.c*(1-P.lambdaw)*(*w)/((*w)+(*m));
	  dpm = -P.nup*(*pm) + P.c*(1-P.lambdam)*(*m)/((*w)+(*m));
	  (*w) += dw*dt;
	  (*m) += dm*dt;
	  (*pw) += dpw*dt;
	  (*pm) += dpm*dt;
	}
    }
}

int main(int argc, char *argv[])
{
  double epsilon;
  Params P;
  int w0, m0;
  double w, m, pw, pm;
  double wgrid[100][100], mgrid[100][100], pwgrid[100][100], pmgrid[100][100];
  int mark0[MAXMARK], mark1[MAXMARK];
  int nmark;
  double thresh;
  FILE *fp;
  double b;
  int r;
  double meanhet;
  int i;
  int expt;
  char str[100];
  int mym0;
  double meanh, normh;

  // general approach is to loop through initial conditions m0 = {0...10} and w0 = {0...100}
  // we compute the final state of the cell for each
  // imposing a given threshold, we identify those initial conditions that give a "power" exceeding this threshold
  // we then sample the final states within this set to compute the final heteroplasmy

  // process command-line argument
  // expt determines the model structure (whether or not replication rates are scaled by total mtDNA content)
  expt = 0;
  if(argc < 2)
    {
      printf("Model structure not specified; using gamma = 0\n");
    }
  else if(atoi(argv[1]) == 1)
    {
      printf("Using gamma = 1\n");
      expt = 1;
    }
  else printf("Using gamma = 0\n");
  
  // lambdaw is wild-type replication rate, against which "foreign selfishness" is measured
  for(P.lambdaw = 0.2; P.lambdaw <= 0.8; P.lambdaw += 0.3)
    {
      // epsilon describes how negative interaction between mtDNAs changes "power" (i.e. how much mixed mtDNAs are a problem)
       for(epsilon = 0; epsilon <= 0.001; epsilon += 0.001)
	{
	  // mym0 is the initial mutant copy number
	  for(mym0 = 10; mym0 <= 90; mym0 += 10)
	    {
		{
		  sprintf(str, "mtdna-cell-%i-%i-%.3f-%.1f.txt", expt, mym0, epsilon, P.lambdaw);
		  fp = fopen(str, "w");

		  // mtDNA and protein degradation rates are set to unit rate in both models 
		  P.nu = P.nup = 1;

		  // set specific activities, timescales, and model labels for the two different model structures
		  if(expt == 0)
		    {
		      P.c = 4;
		      P.t = 10;
		      P.model = 0;
		    }
		  else
		    {
		      P.c = 400;
		      P.t = 1;
		      P.model = 1;
		    }

		  // run a dummy simulation at a fixed parameter set to check performance 
		  P.lambdam = 0.6; w0 = 50; m0 = 5;
		  Simulate(P, w0, m0, &w, &m, &pw, &pm);
		  fprintf(fp, "# %f %f %f %f\n", w, m, pw, pm);

		  // initialise the grids that will be used to store "final" states for each initial condition
		  for(w0 = 0; w0 < 100; w0++)
		    {
		      for(m0 = 0; m0 < 100; m0++)
			wgrid[w0][m0] = mgrid[w0][m0] = pwgrid[w0][m0] = pmgrid[w0][m0] = 0;
		    }

		  // loop through "foreign selfishness" (i.e. mutant replication rate relative to fixed wildtype rate)
		  for(P.lambdam = 0.1; P.lambdam <= 0.9; P.lambdam += 0.01)
		    {

		      printf("Computing %f %f %f %i %i\n", P.lambdaw, P.lambdam, epsilon, mym0, expt);
		  
		      meanh = normh = 0;
		      // loop through initial conditions, and store final state for each
		      for(w0 = 0; w0 < 100-mym0; w0++)
			{
			  for(m0 = 0; m0 < mym0; m0++)
			    {
			      Simulate(P, w0, m0, &w, &m, &pw, &pm);
			      wgrid[w0][m0] = w;
			      mgrid[w0][m0] = m;
			      pwgrid[w0][m0] = pw;
			      pmgrid[w0][m0] = pm;
		  
			      // gracefully handle mtDNA extinction and record heteroplasmy
			      if(w0 == 0)
				{
				  if(m0 != 0) { meanh += 1; normh++; }
				}
			      else
				{
				  meanh += (double)m0 / (m0+w0);
				  normh++;
				}
			    }
			}

		      // loop through "bioenergetic threshold" values
		      for(thresh = 0.; thresh < 2; thresh += 0.01)
			{
			  nmark = 0;

			  // loop through prerecorded set of initial conditions
			  for(w0 = 0; w0 < 100-mym0; w0++)
			    {
			      for(m0 = 0; m0 < mym0; m0++)
				{
				  // "power" is computed from prerecorded simulated values 
				  b = pwgrid[w0][m0] + pmgrid[w0][m0] - epsilon*pwgrid[w0][m0]*pmgrid[w0][m0];
				  if(wgrid[w0][m0] + mgrid[w0][m0] > 0)
				    {
				      // normalise power by number of mtDNAs
				      b /= (wgrid[w0][m0] + mgrid[w0][m0]);
				      if(b > thresh)
					{
					  // if our final power exceeds this threshold, store the initial conditions that gave rise to this
					  mark0[nmark] = w0;
					  mark1[nmark] = m0;
					  nmark++;
					  if(nmark > MAXMARK-1)
					    {
					      printf("too many marks\n");
					      exit(0);
					    }
					}
				    }
				}
			    }

			  // if no initial conditions exceed the threshold, store this fact
			  if(nmark == 0)
			    fprintf(fp, "%f %f -1\n", P.lambdam, thresh);
			  else
			    {
			      // otherwise, sample some of those initial conditions that did exceed the threshold, and compute the final heteroplasmy for each
			      meanhet = 0;
			      for(i = 0; i < NSAMP; i++)
				{
				  r = RND*nmark;
				  meanhet += (mgrid[mark0[r]][mark1[r]])/(wgrid[mark0[r]][mark1[r]] + mgrid[mark0[r]][mark1[r]]);

				  // we should have caught extinction cases, but tell us if something somehow slipped through
				  if(isnan(meanhet))
				    printf("undefined heteroplasmy found\n");
				}
			      // output the final heteroplasmy here
			      fprintf(fp, "%f %f %f %f %i\n", P.lambdam, thresh, meanh/normh, meanhet/NSAMP, nmark);
			    }
			}
		      fprintf(fp, "\n");
		    }
		  fclose(fp);
		}
	    }
	}
    }
  
  return 0;
}
 
