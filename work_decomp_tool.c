#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "inter_work_decomp.h"
#include "tool_func.h"

void Get_para(int argc, char *argv[], Sys_par *pt_sys_info)
{
  int num_input = 9;

  if (argc != num_input){
    fprintf(stderr, "Erro: Need %d input, not %d, check the list\n", 
            num_input - 1, argc - 1);
    exit(1);
  }

  pt_sys_info->stp = atof(argv[1]);
  pt_sys_info->N_conf = atoi(argv[2]);
  pt_sys_info->traj_id = atoi(argv[3]);
  pt_sys_info->num_ho = atoi(argv[4]);
  pt_sys_info->wsq_init = atof(argv[5]);
  pt_sys_info->wsq_final = atof(argv[6]);
  pt_sys_info->coup_const = atof(argv[7]);
  strncpy(pt_sys_info->data_path, argv[8], LEN_FILE_NM);

  pt_sys_info->sys_dim = pt_sys_info->num_ho + 1;
  fprintf(stderr, "******\nSystem information\n"
          "stp	%lf\nN_conf	%d\ntraj_id	%d\n"
          "num_ho	%d\nwsq_init	%lf\nwsq_final	%lf\n"
          "coup_const	%.1f\nsys_dim	%d\n",
          pt_sys_info->stp, pt_sys_info->N_conf, pt_sys_info->traj_id,
          pt_sys_info->num_ho, pt_sys_info->wsq_init, pt_sys_info->wsq_final,
          pt_sys_info->coup_const, pt_sys_info->sys_dim);
  fprintf(stderr, "Current_Dir %s\n******\n", pt_sys_info->data_path);
}

void Get_data_file(char *f_nm_traj, char *f_nm_force, char *f_nm_vel,
                   char *f_nm_Hess,
                   Sys_par sys_info)
{
  char file_dir[LEN_FILE_NM];

  snprintf(file_dir, LEN_FILE_NM, "%s/../result", sys_info.data_path);
  snprintf(f_nm_traj, LEN_FILE_NM, "%s/"
           "tps_ho_%d_coor_%d_wsq_%.0f_%.0f_coup_%.1f.dat",
           file_dir, sys_info.num_ho, sys_info.traj_id, sys_info.wsq_init,
           sys_info.wsq_final, sys_info.coup_const);
  snprintf(f_nm_force, LEN_FILE_NM, "%s/"
           "tps_ho_%d_force_%d_wsq_%.0f_%.0f_coup_%.1f.dat",
           file_dir, sys_info.num_ho, sys_info.traj_id, sys_info.wsq_init,
           sys_info.wsq_final, sys_info.coup_const); 
  snprintf(f_nm_vel, LEN_FILE_NM, "%s/"                                        
           "tps_ho_%d_vel_%d_wsq_%.0f_%.0f_coup_%.1f.dat",                     
           file_dir, sys_info.num_ho, sys_info.traj_id, sys_info.wsq_init,       
           sys_info.wsq_final, sys_info.coup_const);
}

void Read_II_xvf(FILE *f_pt, double **mat, int sys_dim, int n_conf)
{
  int i, j;
  double curr_time;

/*
  // Decimal file
  for (i = 0; i < n_conf; i++){
    fscanf(f_pt, "%lf", &curr_time);

    for (j = 0; j < sys_dim; j++){
      fscanf(f_pt, "%lf", &mat[j][i]);
    }
  }
 */
  // Binary file
  for (i = 0; i < n_conf; i++){
    fread(&curr_time, sizeof(double), 1, f_pt);

    for (j = 0; j < sys_dim; j++){
      fread(&mat[j][i], sizeof(double), 1, f_pt);
    } 
  }
}

void Cal_Hess(double **hess, double **traj_mat, Sys_par sys_info)
{
  int i, j, k; 
  int curr_line;
  double *wsq;
  double sum_rsq_over_wsq;
  FILE *f_pt_wsq;  
  char f_nm_wsq[LEN_FILE_NM];

  wsq = malloc(sys_info.num_ho * sizeof(double));
  snprintf(f_nm_wsq, LEN_FILE_NM, "%s/../ho_info/ho_%d.info", 
           sys_info.data_path, sys_info.num_ho);
  f_pt_wsq = fopen(f_nm_wsq, "r");
  Read_wsq(f_pt_wsq, wsq, sys_info.num_ho);
  fclose(f_pt_wsq);

  sum_rsq_over_wsq = 0.0;
  for (i = 0; i < sys_info.num_ho; i++){
    sum_rsq_over_wsq += pow(sys_info.coup_const, 2.0) / wsq[i]; 
  }

  for (i = 0; i < sys_info.N_conf; i++){
    Cal_one_Hess(hess, traj_mat, sys_info, i, wsq, sum_rsq_over_wsq);
  }

  free(wsq);
}

void Cal_one_Hess(double **hess, double **traj_mat, Sys_par sys_info, 
                  int conf_id, double *wsq, double sum_rsq_over_wsq)
{
  int j, k, curr_line;

  curr_line = conf_id * sys_info.sys_dim;   
  
  for (j = 0; j < sys_info.sys_dim; j++){
    for (k = 0; k < sys_info.sys_dim; k++){
      hess[curr_line + j][k] = 0.0;
    }
  } 

  // Uxx
  hess[curr_line + 0][0] = 
      4.0 * DEPTH_DW * (3.0 * pow(traj_mat[0][conf_id], 2.0) - BASIN) +
      sum_rsq_over_wsq;

  // Ujj (Uij == 0 if i != j)
  for (j = 0; j < sys_info.num_ho; j++){
    hess[curr_line + j + PAR_DIM][j + PAR_DIM] = wsq[j];
  }

  // Uxj
  for (j = 0; j < sys_info.num_ho; j++){
    hess[curr_line + 0][j + PAR_DIM] = -1.0 * sys_info.coup_const;
    hess[curr_line + j + PAR_DIM][0] = hess[curr_line + 0][j + PAR_DIM];
  }
}

void Read_wsq(FILE *f_pt, double *wsq, int num_ho)
{
  int i, readi;

  for (i = 0; i < num_ho; i++){
    fscanf(f_pt, "%d%lf\n", &readi, &wsq[i]);
    fprintf(stderr, "HO %2d %3.3f\n", readi, wsq[i]);
  }
}
