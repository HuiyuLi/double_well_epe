#ifndef _INTER_WORK_DECOMP_H
#define _INTER_WORK_DECOMP_H

#define PAR_DIM                         1
#define MAX_N_RESET			10000
                                                                                 
#define LEN_FILE_NM			300                                       

#define DEPTH_DW			10.0
#define BASIN				1.0 
                                                                                 
typedef struct sys_type{                                                         
  double stp;
  int N_conf;                                                                    
  int traj_id;
  int num_ho;
  double wsq_init;
  double wsq_final;
  double coup_const; 
  char data_path[LEN_FILE_NM];

  int sys_dim;
} Sys_par;                                                                       
                                                                                 
typedef struct reset_type{                                                       
  int* N_reset;                                                                  
  int** reset_list;                                                              
  double** portion_reset;                                                        
} Reset_par;                                                                     
                                                                                 
void Work_decomp_dw(Sys_par, double**, double**,                                 
                    double**, int, int, double**, Reset_par);                    
void Outp_force(FILE*, double*, int, double**, Sys_par);                        
void Outp_tot_w(FILE*, double*, int, double**, Sys_par);

void Get_para(int argc, char *argv[], Sys_par *pt_sys_info);
void Get_data_file(char *f_nm_traj, char *f_nm_force, char *f_nm_vel,
                   char *f_nm_Hess,           
                   Sys_par sys_info);
void Read_II_xvf(FILE *f_pt, double **mat, int sys_dim, int n_conf);
void Cal_Hess(double **hess, double **traj_mat, Sys_par sys_info); 
void Cal_one_Hess(double **hess, double **traj_mat, Sys_par sys_info,
                  int conf_id, double *wsq, double sum_rsq_over_wsq);
void Read_wsq(FILE *f_pt, double *wsq, int num_ho);

#endif
