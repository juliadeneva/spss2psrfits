#include <stdio.h>
#include <stdlib.h>

void readinfo(FILE* f, char* source, char* ra, char* dec, int* imjd, int* smjd, float* offs, int* nchan, float* dt)
{
  char* line;
  size_t len = 0;
  
  //skip directory name
  getline(&line, &len, f);
  free(line); //discard dir name

  fscanf(f, "SOURCE = %s\n", source);
  fscanf(f, "RA = %s\n", ra);
  fscanf(f, "DEC = %s\n", dec);
  fscanf(f, "STT_IMJD = %d\n", imjd);
  fscanf(f, "STT_SMJD = %d\n", smjd);
  fscanf(f, "STT_OFFS = %f\n", offs);
  fscanf(f, "CHAN_PER_FILE = %d\n", nchan);
  fscanf(f, "TBIN = %f\n", dt);
}

void readfreq(FILE* f, int* ifreq, float** freq, int nchan)
{
  int i, ibob;
  
  while (!feof(f)) {
    fscanf(f,"IBOB%d_FREQS = ", &ibob);
    //printf("ibob: %d\n",ibob);

    for (i=0; i<nchan; i++) {
      //printf("ichan: %d\n",i);
      fscanf(f,"%f,",&(freq[ibob][i]));
    }
    fscanf(f,"\n");

    if (freq[ibob][1]-freq[ibob][0] > 0)
      ifreq[ibob] = 1;
    else
      ifreq[ibob] = -1;
  }

}

