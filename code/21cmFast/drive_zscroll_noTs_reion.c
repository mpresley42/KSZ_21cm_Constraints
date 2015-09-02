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
  system("mkdir ../Output_files");
  system("mkdir ../Output_files/DNDLNM_files");
  system("mkdir ../Output_files/FgtrM_files");
  system("mkdir ../Output_files/Halo_lists");
  system("mkdir ../Output_files/Size_distributions");
  system("mkdir ../Output_files/Deldel_T_power_spec");



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


    // find bubbles
    sprintf(cmnd, "./find_HII_bubbles %f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);

    /*
    // generate size distributions, first ionized bubbles
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "gen_size_distr %06.2f 0 ../Boxes/xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z,Z, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	sprintf(cmnd, "gen_size_distr %06.2f 0 ../Boxes/xH_z%06.2f_nohalos_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z,Z, HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "gen_size_distr %06.2f 0 ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	sprintf(cmnd, "gen_size_distr %06.2f 0 ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
  

    // generate size distributions, then neutral regions
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "gen_size_distr %06.2f 1 ../Boxes/xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	sprintf(cmnd, "gen_size_distr %06.2f 1 ../Boxes/xH_nohalos_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "gen_size_distr %06.2f 1 ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	sprintf(cmnd, "gen_size_distr %06.2f 1 ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
    */

    // do temperature map
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_nohalos_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);

    fprintf(stderr, "*************************************\n");
    fflush(NULL);
    Z += ZSTEP;
  }

  fclose(LOG);
  return 0;
}
