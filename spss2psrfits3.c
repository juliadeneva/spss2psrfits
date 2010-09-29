/*
  This program converts data from DSS-28 to psrfits.

  DSS-28 data is typically written as 8 separate files containing dynamic
  spectra as 8-bit values. The files have names iBOBX.spec, where X=0-7.
  Each iBOB may drops samples (complete spectra) occasionally. This program
  assumes that dropped samples are fixed/filled in. 

  There are also iBOBX.idx files. Each of these contains 32-bit integers
  that are sample sequence numbers. If there is a jump by 2 or more in the
  integer sequence, this means that a sample was dropped. This program does
  not read the .idx files, but they are necessary for filling in dropped
  samples.

  This program depends on an info.txt file in the same directory as the .spec
  and .idx files. The info.txt file contains basic pointing information like
  the source name and coordinates, and channel center frequencies for each 
  iBOB. 

  One or more iBOBs may be turned off. If data from one or more iBOBs doesn't
  exist, the weights of the respective channels in the output psrfits file
  are set to 0.0.
  

  TO DO:

  If the user decides to ignore data from one or more iBOBs, but the data
  exists, it will be written to the output psrfits file but the channel
  weights for channels from the ignored iBOB will be set to 0.0.
  

 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chkio.h"
#include "vectors.h"
#include "spss2psrfits_cmd.h"
#include "psrfits.h"


void readinfo(FILE* f, char* source, char* ra, char* dec, int* imjd, int* smjd,float* offs, int* nchan, float* dt);
void readfreq(FILE* f, int* ifreq, float** freq, int nchan);
void fill_psrfits_struct(int numbits, struct psrfits *pf, float dt, char* source, int nchantot, float df, int imjd, int smjd, float offs, char* ra, char* dec, float lofreq);

int main(int argc, char *argv[])
{
  int numfiles, ii, rownum, ichan, datidx, breakout;
  int spec_per_row, specnum, status, bytes_read, bytes_per_subint_oneibob;
  unsigned char **rawdata;
  FILE * infofile;
  FILE **specfiles;
  char specname[128], infoname[128];
  char source[128], ra[16], dec[16], tmpchar;
  int imjd, smjd, nchan, *ifreq, nchantot, nchantot_fixed, numbits;
  float tbin, **freq, df, flow, fhigh, offs;
  struct psrfits pf;
  Cmdline *cmd;
  int* ignored, *exist; //which ibobs should be ignored and which have data
  int numfilelow, numfilehigh;

  numbits = 8;
  numfiles = 8; //iBOBs 0-7
  numfilelow = -1;
  numfilehigh = -1;

  if (argc == 1) {
    Program = argv[0];
    usage();
    exit(1);
  }

    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);
#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n          SPSS to PSRFITs Conversion Code\n");
    printf("                     by J. Deneva\n\n");
    //printf("Code has git hash:  '%s'\n\n", GITHASH);

    specfiles = (FILE **) malloc(numfiles * sizeof(FILE *));
    freq = (float **) malloc(numfiles * sizeof(float *));
    ifreq = gen_ivect(numfiles);
    ignored = gen_ivect(numfiles);
    exist = gen_ivect(numfiles);

    //Which ibobs will be ignored
    if (cmd->ignored != NULL) {
      ii = 0;
      //printf("Will ignore iBOBs ");
      while (cmd->ignored[ii] != '\0')  {
	//Typecasting a single character to int converts e.g. '9'(char)
	//and 9(int) via subtracting 48
	//printf("%d ",(int)(cmd->ignored[ii])-48);
	ignored[(int)(cmd->ignored[ii])-48] = 1;
	ii++;
      }
      //printf("\n");
      //printf("(If data exists, it will be written to output, but channel weights will be set to zero.)\n");
    }    

    // Open the input files
    printf("Reading input data from:\n");
    sprintf(infoname,"%s%s",cmd->argv[cmd->argc-1],"/info.txt");
    infofile = chkfopen(infoname,"r");
    readinfo(infofile,source,ra,dec,&imjd,&smjd,&offs,&nchan,&tbin);
    for (ii = 0; ii < numfiles; ii++)
      freq[ii] = gen_fvect((long)nchan);
    readfreq(infofile,ifreq,freq,nchan);
    printf("imjd: %d smjd: %d offs: %f\n", imjd, smjd, offs);

    flow = 10000;
    fhigh = 0;
    for (ii = 0; ii < numfiles; ii++) {
        printf("  '%s'\n", cmd->argv[cmd->argc-1]);
	sprintf(specname,"%s%s%d%s",cmd->argv[cmd->argc-1],"/iBOB",ii,"fixed.spec");
        specfiles[ii] = fopen(specname,"rb");

	if (specfiles[ii] != NULL) {
	  exist[ii] = 1;
	  
	  if (ignored[ii] == 1) {
	    printf("File %s exists, user chose to ignore.\n",specname);
	    continue;
	  }

	  if (numfilelow < 0)
	    numfilelow = ii;
	  numfilehigh = ii;

	  //If channel width is negative in psrfits headers, Presto 
	  //flips the band before dedispersing. 
	  df = fabs(freq[ii][1] - freq[ii][0]);

	  printf("ii: %d freq[0]: %f freq[nchan-1]: %f\n",ii,freq[ii][0],freq[ii][nchan-1]);

	  if (ifreq[ii] > 0) {
	    if (freq[ii][0] < flow)
	      flow = freq[ii][0];
	    if (freq[ii][nchan-1] > fhigh)
	      fhigh = freq[ii][nchan-1];
	  }
	  else if (ifreq[ii] < 0) {
	    if (freq[ii][nchan-1] < flow)
	      flow = freq[ii][nchan-1];
	    if (freq[ii][0] > fhigh)
	      fhigh = freq[ii][0];
	  }

	} else {
	  exist[ii] = 0;
	  printf("File %s doesn't exist.\n",specname);
	}
    }
    
    //Total number of channels in output data
    nchantot =(int)((fhigh-flow)/fabs(df)) + 1;

    printf("numfilelow: %d numfilehigh: %d\n", numfilelow, numfilehigh);
    fprintf(stderr, "flow: %f fhigh: %f nchantot: %d\n", flow, fhigh, nchantot);

    //Fix the total number of channels in output file to be the power of 2
    //closest to actual # of (non-overlapping) channels from all ibobs
    nchantot_fixed = 32;
    while (nchantot_fixed < nchantot)
      nchantot_fixed = nchantot_fixed * 2;

    printf("nchantot_fixed: %d\n", nchantot_fixed);

    //Initialize the psrfits struct
    fill_psrfits_struct(numbits, &pf, tbin, source, nchantot_fixed, df, imjd, smjd, offs, ra, dec, flow);
    //Channel weights are set to 1.0 by fill_psrfits_struct;
    //Reset weights of excess channels to zero.
    for (ii = nchantot; ii < nchantot_fixed; ii++)
      pf.sub.dat_weights[ii] = 0;

    //Center frequencies of each channel are set by fill_psrfits_struct.
    //Set weights to zero for iBOBs in the middle that have no data
    for (ii = numfilelow; ii <= numfilehigh; ii++) {
      for (ichan = 0; ichan < nchan; ichan++) {
	datidx = (int)(roundf(freq[ii][ichan]-flow)/2.0);
	//pf.sub.dat_freqs[datidx] = freq[ii][ichan];
	//printf("ii: %d datidx: %d freq: %f\n",ii,datidx,freq[ii][ichan]);

	//to zero as well. 
	if (exist[ii] == 0 || ignored[ii] == 1)
	  pf.sub.dat_weights[datidx] = 0.0;
      }
    }
    
    spec_per_row = pf.hdr.nsblk;
    pf.multifile = 1;
    pf.rows_per_file = (int) ((cmd->outlenGB * 1073741824.0)
			      / (double) pf.sub.bytes_per_subint);
    bytes_per_subint_oneibob = spec_per_row * nchan;

    printf("PSRFITS file(s) will have %d samples per row and %d rows per file.\n\n", spec_per_row, pf.rows_per_file);
    printf("pf.sub.bytes_per_subint: %d\n", pf.sub.bytes_per_subint);
    printf("bytes_per_subint_oneibob: %d\n", bytes_per_subint_oneibob);

    // Raw data buffer--to be filled from input file
    rawdata = calloc((numfilehigh-numfilelow+1), sizeof(unsigned char*));
    for (ii = numfilelow; ii <= numfilehigh; ii++) {
      rawdata[ii-numfilelow] = gen_bvect(spec_per_row*nchan);
    }
    
    // Create the PSRFITS file
    strcpy(pf.basefilename, cmd->outfile);
    psrfits_create(&pf);

    // Loop over the data
    rownum = 0;
    breakout = 0;
    while (1) {

      printf("\rWorking on row %d\r", ++rownum);
      fflush(stdout);

      //Read in a row's worth of data from each input file.
      //If not all files have that much data left, exit loop.
      for (ii = numfilelow; ii <= numfilehigh; ii++) {

	//Read in one row's worth of data
	if (exist[ii] == 1 && ignored[ii] == 0) {
	  bytes_read = fread(rawdata[ii-numfilelow], sizeof(unsigned char), bytes_per_subint_oneibob, specfiles[ii]);
	  //fprintf(stderr,"ibob: %d bytes_read: %d \n", ii, bytes_read);
	  //If partial last row, discard it
	  if (bytes_read < bytes_per_subint_oneibob)
	    {
	      breakout = 1;
	      break;
	    }
	}
      }
      
      if (breakout == 1)
	break; //out of while loop and close output file

      for (ii = numfilelow; ii <= numfilehigh; ii++) {
	
	//If there's no data from current iBOB or if it's ignored
	if (exist[ii] == 0 || ignored[ii] == 1) 
	  continue;
      
	// Loop over all the spectra per row
	for (specnum = 0; specnum < spec_per_row; specnum++) {
	  for (ichan = 0; ichan < nchan; ichan++) {
	    //Where in output buffer should current point go
	    datidx = specnum*nchantot_fixed + (int)(roundf(freq[ii][ichan]-flow)/2.0);
	    if (datidx > pf.sub.bytes_per_subint)
	      printf("ibob: %d specnum: %d ichan: %d datidx: %d\n", ii,specnum,ichan,datidx);
	    
	    pf.sub.data[datidx] = rawdata[ii-numfilelow][specnum*nchan + ichan];
	    
	  }
	}
      }
	
      // Now write the row...
      pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
      //printf("pf.sub.offs: %f\n",pf.sub.offs);
      status = psrfits_write_subint(&pf);
      if (status) {
	printf("\nError (%d) writing PSRFITS...\n\n", status);
	break;
      }
      
    }

    printf("\n");
    
    // Close the PSRFITS file
    psrfits_close(&pf);

    // Free the structure arrays too...
    free(pf.sub.dat_freqs);
    free(pf.sub.dat_weights);
    free(pf.sub.dat_offsets);
    free(pf.sub.dat_scales);
    free(pf.sub.data);
    //free(rawdata);
    free(ignored);
    free(exist);

    return 0;
}
