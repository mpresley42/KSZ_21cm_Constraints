#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "../Parameter_files/INIT_PARAMS.H"

/*
  Program DRIVE_ZSCROLL_NOTS.C scrolls through the redshifts defined below,
  creating halo, evolved density, velocity, 21cm fields.
  NOTE: this driver assumes that the IGM has already been heated to Ts>>Tcmb.
  If you wish to compute the spin temperature, use the other driver.
*/

/*
  redshift scrolling parameters used by the drive_zscroll drive program
*/
//#define ZSTART (15.2) //inclusive
//#define ZSTART (16.6) //inclusive
//#define ZEND (8.1999999) // inclusive
//#define ZSTEP (-0.2)
#define ZSTART (30.0) //inclusive
#define ZEND (6.0) // inclusive
#define ZSTEP (-0.5)


int main(int argc, char ** argv){
  float Z, M, M_MIN;
  char cmnd[1000];
  FILE *LOG;
  time_t start_time, curr_time;

  time(&start_time);

  // make appropriate directories
  system("mkdir ../Log_files");
  system("mkdir ../Boxes");

  // open log file
  system("rm ../Log_files/*");
  LOG = fopen("../Log_files/drive_zscroll_noTs_log_file", "w");
  if (!LOG){
    fprintf(stderr, "drive_zscroll_log_file.c: Unable to open log file\n Aborting...\n");
    return -1;
  }

  Z = ZSTART;
  while (Z > (ZEND-0.0001)){
    fprintf(stderr, "*************************************\n");

    M_MIN = get_M_min_ion(Z);

    // if USE_HALO_FIELD is turned on in ANAL_PARAMS.H, run the halo finder
    if (USE_HALO_FIELD){
      //  the following only depend on redshift, not ionization field
      // find halos
      sprintf(cmnd, "./find_halos %f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);


      // shift halos accordig to their linear velocities
      sprintf(cmnd, "./update_halo_pos %f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
    }

    // shift density field and update velocity field
    sprintf(cmnd, "./perturb_field %f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
    // end of solely redshift dependent things, now do ionization stuff


    fprintf(stderr, "*************************************\n");
    fflush(NULL);
    Z += ZSTEP;
  }

  fclose(LOG);
  return 0;
}
