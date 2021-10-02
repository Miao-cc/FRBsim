// gcc -O3 -lm -o simulateSimplePsr simulateSimplePsr.c
// This routine simulates a simple pulsar that is defined by a single periodicity
// (i.e., this routine can not be used when modelling Doppler effects from the Earth's motion or
// from orbital effects).
// This routine also assumes that the dispersion smearing across the band is less than half the pulse period


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include <string.h>
#include "simulate.h"
#include <omp.h>
#include "ran1.c"
#include "gasdev.c"
#include "time.h"
#include "bessi0.c"

int main(int argc,char *argv[])
{
  header *head;
  int i,j,k;
  float amp = 3;
  double *sum;
  float si = -2;
  float e = 2.718281828 ;
  
  float width; // smear width 
  float *profile;
  int npsamp;
  int nsamp;
  float tval,phase;
  float fval;
  int ival;
  float t;
  float *dt_dm;
  float fref;
  float chanbw;
  float *chanfref;
  float minfref = 244.140625;      // the profile will be 0 in those channels whose fref < minfref (by mcc)
  /*float minfref = 0.;*/
  int ncap;
  float kappa;
  
  long pn=0; // Pulse number
  FILE *fout;
  char outName[128] = "pulsar.dat";
  char description[128];
  char format[128];
  int useParamFile=0;
  char paramFile[MAX_PARAM_FILES][1024];
  int nperiods;
  int delay;


  head = (header *)malloc(sizeof(header));
  simulateSetHeaderDefaults(head);

  // Read input parameters
  for (i=0;i<argc;i++)
{
      if (strcmp(argv[i],"-o")==0)
strcpy(outName,argv[++i]);
      else if (strcmp(argv[i],"-p")==0)
	{
	  strcpy(paramFile[useParamFile],argv[++i]);
	  useParamFile++;
	}
    }
  if (useParamFile>0)
    {
      for (i=0;i<useParamFile;i++)
	{
	  printf("Reading parameter file: %s\n",paramFile[i]);
	  simulateReadParamFile(head,paramFile[i]);
	}
    }
  //
  npsamp = (int)(head->p0/head->tsamp);     //samples in a period
  nsamp = (int)ceil((head->t1 - head->t0)/head->tsamp);  //samples in all file length
  //nperiods = (int)ceil((head->t1 - head->t0)/head->p0);
  nperiods = 1; // change to create one preiod simulation file
  dt_dm = (float *)malloc(sizeof(float)*head->nchan);
  chanfref = (float *)malloc(sizeof(float)*head->nchan);
  sum = (double *)malloc(sizeof(double)*head->nchan);
  chanbw = (head->f2-head->f1)/head->nchan;
  profile = (float *)malloc(sizeof(float)*nsamp*head->nchan);

  printf(" sample per period: %d\n All samples: %d\n Number of periods: %d\n  chanbw:%f\n", npsamp, nsamp, nperiods, fabs(chanbw));

  fref = 0.5*(head->f2+head->f2);

  //open file && write header

  //change
  //----------------------------------------------------------------------------------------------------
  // calculate the time smear in each channel 
  for (i=0;i<head->nchan;i++)
  {
    chanfref[i] = (head->f1+i*chanbw);
    if (chanfref[i] <= minfref)
      {
        dt_dm[i]=0.;
      }
    else
      {
        dt_dm[i] = (4.15e-3*head->dm*(pow(fref/1000.0,si)-pow(chanfref[i]/1000.0,si)));
      }
   }
  printf("Time delay (sec): %g, dm: %g,fref: %g,f1: %g MHz\n",dt_dm[0],head->dm,fref,head->f1);


  delay = pow(2,10)*1;
  //delay = head->t1 - dt_dm
  printf("Time delay offset (sec): %g\n",head->tsamp*delay);


// Make a simple profile
/*#pragma omp parallel for default(shared) private(i,j) shared(profile, sum)*/
    for (j=0;j<head->nchan;j++)
      {
        sum[j] = 0.;
        width = (fabs(2*dt_dm[j]*chanbw/(chanfref[j])) + head->width) * 2 * 3.14159265;
        if (j % 100 == 0 ) printf("j: %d width: %f ,dt_dm: %f ,fref: %f\n",j, width,dt_dm[j],chanfref[j]);

/*#pragma omp parallel for default(shared) private(i) shared(profile, sum)*/
        for (i=0 ;i<nsamp;i++) //from t=0 ~ t=t1
          {
            //tval = i*head->tsamp + dt_dm[j] ;
            //add time delay (4 subint). make sure the pulse showed completely
            //tval = i*head->tsamp + dt_dm[j] + head->tsamp*delay;
            tval = i*head->tsamp + dt_dm[j];
            phase = fmod(tval/head->p0, 1.);
            if (phase < 0.) phase += 1.;
            phase -= fmod(head->tsamp*delay/head->p0, 1.);

            phase *= 2*3.14159265;
            /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/head->width/head->width);*/
            profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width);
          }
  }

  fout = fopen(outName,"wb");
  fwrite(profile,sizeof(float),nsamp*head->nchan,fout);
  fclose(fout);

  simulateReleaseMemory(head);
  free(profile);
  free(sum);
  free(dt_dm);
  free(chanfref);
  free(head);
}
