#ifndef _TOOL_FUNC_H
#define _TOOL_FUNC_H

// FORWARD 0, BACKWARD 1, BOTH 2
#define FORWARD		0

void Segment_df_integ(double* res_vec_integ, double* inside_vec, int tot_len,    
                      int N_reset, int* curr_reset_list, double* portion_reset,  
                      double stp, double* portion_init_value,
                      double *be_shift_f, double *af_shift_f);

void Do_segment_integ(double *res_vec_integ, double *inside_vec,                 
                      int init_conf, int end_conf, double stp,               
                      double init_value, int direction);

void Do_integ_simpson(double *res_vec_integ, double *inside_vec,                 
                     int init_conf, int end_conf, double stp,                    
                     double init_value);

void Do_integ_reverse_simpson(double *res_vec_integ, double *inside_vec,         
                              int init_conf, int end_conf, double stp,           
                              double init_value);                                
                                                 

void Do_simpson_piecewise(double *res_vec_integ, double *inside_vec,             
                          int init_conf, int end_conf, double stp,               
                          double init_value, 
                          int *reset_list, int num_reset,                     
                          double *portion_reset, double *before_work,            
                          double *after_work, double *x_dot,                     
                          double *force_comp, int curr_dim, int df,
                          double *U_ab_b_dot,
                          double *be_shift_f, double *af_shift_f,
                          double *before_force, double *after_force);

void Do_reset_interval_integ(double *res_vec_integ,                              
                             double stp,                                          
                             double *portion_reset, double *before_work,          
                             double *after_work, double *x_dot,               
                             double *force_comp,                             
                             double *be_shift_f, double *af_shift_f,              
                             double *before_force, double *after_force,           
                             int m, int reset_id);

int Whe_in_list(int val, int *vector, int len);

void Get_shift_f(double *be_shift_f, double *af_shift_f,                         
                 double* inside_vec, int N_reset,                                
                 double* curr_portion_reset,                                     
                 int* curr_reset_list, double t_stp);

void Cal_portion_init_value(double *portion_init_value,                          
                            double* inside_vec, int N_reset,                     
                            double* curr_portion_reset,                          
                            int* curr_reset_list, double t_stp);

void New_Count_row_mat_cross_num(double** mat, int de_row, int de_col,           
                              int* N_cross, double value,                        
                              int max_num_reset, double** portion_reset,         
                              int **cross_list);

double** Malloc_mat(int de1, int de2); 
void Free_mat(double** mat, int ind);
int** Malloc_int_mat(int de1, int de2); 
void Free_int_mat(int** mat, int ind); 
void Init_mat(double** mat, int de1, int de2, double value); 
void Init_vec(double* vec, int len, double value); 
void Init_int_vec(int* vec, int len, int value); 


#endif
