/*
 *  ******************************************************************
 *  
 *  ala_work_decomp.c
 *  (C) by Huiyu Li Jan 31 2017
 * 
 *  work decomposition for one dimential particle coupling to HO 
 *  heat bath
 *
 *  Data:
 *  1. XVF for each configuration
 *  2. Parameters for MD (stp, N_conf, energy function)
 *
 *  ******************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "inter_work_decomp.h"
#include "tool_func.h"

#define WHE_II			0

int main(int argc, char* argv[])
{
  FILE* f_II_pt_traj;
  FILE* f_II_pt_force;
  FILE* f_II_pt_Hess;
  char f_II_nm_traj[LEN_FILE_NM];
  char f_II_nm_force[LEN_FILE_NM];
  char f_II_nm_Hess[LEN_FILE_NM];  

  Sys_par sys_info;
  Reset_par reset_info; 
  char data_path[LEN_FILE_NM];

  double** traj_mat;
  double** force_mat;
  double** x_dot_mat;
  double** U_xy;

  int i, j;
  int curr_dim_id;

  FILE *f_II_pt_vel;
  char f_II_nm_vel[LEN_FILE_NM];
  int check_id;

  /* Initialize system information */
  Get_para(argc, argv, &sys_info);
  
  /* Get the responding  data file */   
  Get_data_file(f_II_nm_traj, f_II_nm_force, f_II_nm_vel, f_II_nm_Hess, sys_info);
  f_II_pt_traj = fopen(f_II_nm_traj, "r");
  f_II_pt_force = fopen(f_II_nm_force, "r");
  f_II_pt_vel = fopen(f_II_nm_vel, "r") ;
 
  if ((f_II_pt_traj == NULL) || (f_II_pt_force == NULL)){
    fprintf(stderr, "Error in open data file\n");
    exit(1);
  } else{
    fprintf(stderr, "%s\n%s\n", f_II_nm_traj, f_II_nm_force);
  }

  /* Read traj/force/Hess file, calculate x_dot, and U_xy */
  traj_mat = Malloc_mat(sys_info.sys_dim, sys_info.N_conf);
  force_mat = Malloc_mat(sys_info.sys_dim, sys_info.N_conf);
  U_xy = Malloc_mat(sys_info.N_conf * sys_info.sys_dim, sys_info.sys_dim);

  Read_II_xvf(f_II_pt_traj, traj_mat, sys_info.sys_dim, sys_info.N_conf); 
  Read_II_xvf(f_II_pt_force, force_mat, sys_info.sys_dim, sys_info.N_conf);

  Cal_Hess(U_xy, traj_mat, sys_info);

  fclose(f_II_pt_traj);
  fclose(f_II_pt_force);

  x_dot_mat = Malloc_mat(sys_info.sys_dim, sys_info.N_conf);
  Read_II_xvf(f_II_pt_vel, x_dot_mat, sys_info.sys_dim, sys_info.N_conf);
  fclose(f_II_pt_vel);

  fprintf(stderr, "begin reset\n");

  /* By force, count reset points*/
  reset_info.N_reset = malloc(sys_info.sys_dim * sizeof(int));
  reset_info.portion_reset = Malloc_mat(sys_info.sys_dim, MAX_N_RESET);
  reset_info.reset_list = Malloc_int_mat(sys_info.sys_dim,  MAX_N_RESET);

  New_Count_row_mat_cross_num(force_mat,                     
        sys_info.sys_dim, sys_info.N_conf, reset_info.N_reset, 0.0, MAX_N_RESET,   
        reset_info.portion_reset, reset_info.reset_list);

  /* Do work decomposition */
  check_id = 0;
  for (i = 0; i <  reset_info.N_reset[check_id]; i++){                     
    fprintf(stderr, "%d %d   %lf   %lf   %lf\n", i, 
            reset_info.reset_list[check_id][i], 
            sys_info.stp * reset_info.reset_list[check_id][i], 
            traj_mat[check_id][reset_info.reset_list[check_id][i]],
            traj_mat[check_id][reset_info.reset_list[check_id][i] + 1]);
  } 
//  for (j = 0; j < sys_info.sys_dim; j++){
  for (j = check_id; j < (check_id + 1); j++){
    fprintf(stderr, "df %d N_reset %d\n", j, reset_info.N_reset[j]);
    Work_decomp_dw(sys_info, traj_mat, force_mat, U_xy, 
                   sys_info.sys_dim, j, x_dot_mat, reset_info);
  }

  Free_mat(traj_mat, sys_info.sys_dim);
  Free_mat(force_mat, sys_info.sys_dim);
  Free_mat(U_xy, sys_info.N_conf * sys_info.sys_dim);
  Free_mat(x_dot_mat, sys_info.sys_dim);

  Free_int_mat(reset_info.reset_list, sys_info.sys_dim);
  free(reset_info.N_reset);
  free(reset_info.portion_reset);

  return 0;
}


void Work_decomp_dw(Sys_par sys_info, double** traj_mat, 
                    double** force_mat, double** U_xy, int sys_dimension, 
                    int curr_dim_id, double** x_dot_mat, 
                    Reset_par reset_info)
{
  FILE* f_pt_outp_force;
  FILE* f_pt_outp_work;
  char f_nm_outp_force[LEN_FILE_NM];
  char f_nm_outp_work[LEN_FILE_NM];
  char res_path[LEN_FILE_NM];

  double curr_f;

  int i, j, k;
  double** fst_integ;
  double** snd_integ;
  double* U_ab_b_dot;
  double* portion_init_value;
  double* fst_integ_a_dot;
  double* force_a_dot;
  double* tot_w_by_f;
  double* tot_w_by_f_no_reset;

  FILE *f_pt_reset;
  char f_nm_reset[LEN_FILE_NM];
  FILE *f_pt_outp_reset_force;
  char f_nm_outp_reset_force;

 
  double **before_work, **after_work;
  double **before_shift_f, **after_shift_f;
  double **before_force, **after_force;

  FILE *f_pt_integrand;                                                          
  char f_nm_integrand[LEN_FILE_NM];                                              
  int index_reset, reset_conf;
 
  /* Get result path */
  snprintf(res_path, LEN_FILE_NM, "%s/result", sys_info.data_path);

  /* Main step 1 
   * Use traj coordinate and U_xy to do integral
   * compare the result with force
   */
  fprintf(stderr, "  First integral\n");
  snprintf(f_nm_outp_force, LEN_FILE_NM, 
           "%s/force/force_%d.dat", res_path, curr_dim_id);
  f_pt_outp_force = fopen(f_nm_outp_force, "w");
  U_ab_b_dot = malloc(sys_info.N_conf*sizeof(double));
  fst_integ = Malloc_mat(sys_dimension, sys_info.N_conf);

  before_shift_f = Malloc_mat(sys_dimension, reset_info.N_reset[curr_dim_id]);
  after_shift_f = Malloc_mat(sys_dimension, reset_info.N_reset[curr_dim_id]);

  for (j = 0; j < sys_dimension; j++){
    for (i = 0; i < sys_info.N_conf; i++){
      U_ab_b_dot[i] = U_xy[curr_dim_id + sys_dimension*i][j] * 
                      x_dot_mat[j][i];  
    }

    portion_init_value = malloc(reset_info.N_reset[curr_dim_id] * sizeof(double));
/*
    Cal_portion_init_value(portion_init_value, 
             U_ab_b_dot, reset_info.N_reset[curr_dim_id],
             reset_info.portion_reset[curr_dim_id],
             reset_info.reset_list[curr_dim_id], sys_info.stp);
 */
    Get_shift_f(before_shift_f[j], after_shift_f[j], 
                U_ab_b_dot, reset_info.N_reset[curr_dim_id],
                reset_info.portion_reset[curr_dim_id],
                reset_info.reset_list[curr_dim_id], sys_info.stp);

/*
    Segment_df_integ(fst_integ[j], U_ab_b_dot, sys_info.N_conf,
                     reset_info.N_reset[curr_dim_id],   
                     reset_info.reset_list[curr_dim_id], 
                     reset_info.portion_reset[curr_dim_id],
                     sys_info.stp, portion_init_value,
                     before_shift_f[j], after_shift_f[j]);
 */
    Segment_df_integ(fst_integ[j], U_ab_b_dot, sys_info.N_conf,                  
                     1,                            
                     reset_info.reset_list[curr_dim_id],                         
                     reset_info.portion_reset[curr_dim_id],                      
                     sys_info.stp, portion_init_value,                           
                     before_shift_f[j], after_shift_f[j]);


    free(portion_init_value);
  }

  Outp_force(f_pt_outp_force, force_mat[curr_dim_id], curr_dim_id,
             fst_integ, sys_info);

  /* Main step 2
   * Use first integral result and traj coordinate to do integral
   * compare the result with the work by integral force
   */
  fprintf(stderr, "  Second integral\n");
  snprintf(f_nm_outp_work, LEN_FILE_NM,
           "%s/work/work_%d.dat", res_path, curr_dim_id);
  f_pt_outp_work = fopen(f_nm_outp_work, "w");
  fst_integ_a_dot = malloc(sys_info.N_conf*sizeof(double));
  snd_integ = Malloc_mat(sys_dimension, sys_info.N_conf);

  before_work = Malloc_mat(sys_dimension, reset_info.N_reset[curr_dim_id]);
  after_work = Malloc_mat(sys_dimension, reset_info.N_reset[curr_dim_id]);

  before_force = Malloc_mat(sys_dimension, reset_info.N_reset[curr_dim_id]);
  after_force = Malloc_mat(sys_dimension, reset_info.N_reset[curr_dim_id]);

//  for (j = 0; j < 6; j++){
  for (j = 0; j < sys_dimension; j++){
    for (i = 0; i < sys_info.N_conf; i++){
      fst_integ_a_dot[i] = fst_integ[j][i] * x_dot_mat[curr_dim_id][i]; 
    }

    for (i = 0; i < sys_info.N_conf; i++){                                       
      U_ab_b_dot[i] = U_xy[curr_dim_id + sys_dimension*i][j] *                   
                      x_dot_mat[j][i];                                           
    } 

    Do_simpson_piecewise(snd_integ[j], fst_integ_a_dot, 0, sys_info.N_conf-1,
        sys_info.stp, 0.0, reset_info.reset_list[curr_dim_id], 
        reset_info.N_reset[curr_dim_id], reset_info.portion_reset[curr_dim_id],
        before_work[j], after_work[j], x_dot_mat[curr_dim_id],
        fst_integ[j], curr_dim_id, j, U_ab_b_dot, 
        before_shift_f[j], after_shift_f[j], before_force[j], after_force[j]);

    // Output check information
    snprintf(f_nm_integrand, LEN_FILE_NM, "./result/integrand_%d.dat", j);       
    f_pt_integrand = fopen(f_nm_integrand, "w");

    for (k = 1; k < reset_info.N_reset[curr_dim_id] - 1; k++){    
      reset_conf = reset_info.reset_list[curr_dim_id][k];                        
      fprintf(f_pt_integrand, 
              "%d %12.8e  %12.8e  %12.8e %12.8e %12.8e %12.8e\n",  
              k, before_work[j][k], after_work[j][k],
              snd_integ[j][reset_conf+1] - snd_integ[j][reset_conf],
              reset_info.portion_reset[curr_dim_id][k],
              before_force[j][k], after_force[j][k]);           
    }

    fclose(f_pt_integrand);
  } 
 
  fprintf(stderr, "Done\n") ;
  /* Do integral for force to calculate total work */
  force_a_dot = malloc(sys_info.N_conf*sizeof(double));
  for (i = 0; i < sys_info.N_conf; i++){
    force_a_dot[i] = force_mat[curr_dim_id][i] * x_dot_mat[curr_dim_id][i];
  }

  /* Total work without reset */
  tot_w_by_f_no_reset = malloc(sys_info.N_conf * sizeof(double));
  Do_integ_simpson(tot_w_by_f_no_reset, force_a_dot, 0, sys_info.N_conf-1,
                   sys_info.stp, 0.0);
  Outp_tot_w(f_pt_outp_work, tot_w_by_f_no_reset, curr_dim_id, 
             snd_integ, sys_info);

  // Output reset information
  snprintf(f_nm_reset, LEN_FILE_NM, "%s/work/reset_%d.dat", res_path, curr_dim_id);
  f_pt_reset = fopen(f_nm_reset, "w");                                           
  fwrite(&reset_info.N_reset[curr_dim_id], sizeof(int), 1, f_pt_reset);
  for (i = 0; i < reset_info.N_reset[curr_dim_id]; i++){                         
    fwrite(&i, sizeof(int), 1, f_pt_reset);
    fwrite(&reset_info.reset_list[curr_dim_id][i], sizeof(int), 1, f_pt_reset);
    fwrite(&reset_info.portion_reset[curr_dim_id][i], sizeof(double), 
           1, f_pt_reset);
    for (j = 0; j < sys_dimension; j++){
      fwrite(&before_work[j][i], sizeof(double), 1, f_pt_reset);
      fwrite(&after_work[j][i], sizeof(double), 1, f_pt_reset);
    }
  }                                                                              
  fclose(f_pt_reset);

  // Free space
  Free_mat(before_work, sys_dimension);
  Free_mat(after_work, sys_dimension);
  Free_mat(before_shift_f, sys_dimension);
  Free_mat(after_shift_f, sys_dimension);
  Free_mat(before_force, sys_dimension);
  Free_mat(after_force, sys_dimension);

  Free_mat(fst_integ, sys_dimension);
  free(U_ab_b_dot);
  free(force_a_dot);
  free(tot_w_by_f_no_reset);
}

void Outp_force(FILE* f_pt_outp, double* force_vec, int curr_id, 
                double** integ_result, Sys_par sys_info)
{
  int i, j;
  double sum_fw;

  double curr_f;
  double curr_time;  

  for (i = 0; i < sys_info.N_conf; i++){
    sum_fw = 0;
    for (j = 0; j < sys_info.sys_dim; j++){
      sum_fw = sum_fw + integ_result[j][i];
    }

    curr_time = sys_info.stp*i;
    curr_f = (-1)*force_vec[i];

    if (WHE_II == 0){
      fprintf(f_pt_outp, "%lf %lf %lf ", curr_time, curr_f, sum_fw);
      for (j = 0; j < sys_info.sys_dim; j++){
        fprintf(f_pt_outp, "%lf ", integ_result[j][i]);
      }
      fprintf(f_pt_outp, "\n");
    } else if (WHE_II == 1){
      fwrite(&curr_time, sizeof(double), 1, f_pt_outp);
      fwrite(&curr_f, sizeof(double), 1, f_pt_outp);
      fwrite(&sum_fw, sizeof(double), 1, f_pt_outp);

      for (j = 0; j < sys_info.sys_dim; j++){
        fwrite(&integ_result[j][i], sizeof(double), 1, f_pt_outp);
      }
    } else{                                                                      
      fprintf(stderr, "WHE_II should be 0 or 1 for text file or binary file\n"); 
      exit(1);                                                                   
    }
  }
} 

void Outp_tot_w(FILE* f_pt_outp, double* tot_w_without_reset, int curr_id, 
                double** integ_result, Sys_par sys_info)
{
  int i, j;
  double sum_w;

  double curr_tot, curr_reset, curr_time;

  for (i = 0; i < sys_info.N_conf; i++){
    sum_w = 0.0;
    for (j = 0; j < sys_info.sys_dim; j++){
      sum_w = sum_w + integ_result[j][i];
    }

    curr_time = sys_info.stp*i;
    curr_reset = (-1)*tot_w_without_reset[i];

    if (WHE_II == 0){
      fprintf(f_pt_outp, "%lf %lf %lf ", curr_time, curr_reset, sum_w); 
      for (j = 0; j < sys_info.sys_dim; j++){                                      
        fprintf(f_pt_outp, "%lf ", integ_result[j][i]);                            
      }                                                                            
      fprintf(f_pt_outp, "\n"); 
    } else if (WHE_II == 1){
      fwrite(&curr_time, sizeof(double), 1, f_pt_outp);
      fwrite(&curr_reset, sizeof(double), 1, f_pt_outp);
      fwrite(&sum_w, sizeof(double), 1, f_pt_outp);

      for (j = 0; j < sys_info.sys_dim; j++){
        fwrite(&integ_result[j][i], sizeof(double), 1, f_pt_outp);
      } 
    } else{
      fprintf(stderr, "WHE_II should be 0 or 1 for text file or binary file\n");
      exit(1);
    }

  }
}


