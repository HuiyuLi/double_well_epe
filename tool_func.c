#include <stdio.h>
#include <stdlib.h>
#include "tool_func.h"

void Segment_df_integ(double *res_vec_integ, double *inside_vec, int tot_len,
                      int N_reset, int* curr_reset_list, double *portion_reset,
                      double stp, double* portion_init_value, 
                      double *be_shift_f, double *af_shift_f)
{
  int m;
  int init_conf, end_conf;
  double curr_init_value;

  if (N_reset <= 0){
    fprintf(stderr, "Error: Check traj and reset %d\n", N_reset);
    exit(1);
  }  

  // Segment before the reset reset point
  init_conf = 0;
  end_conf = curr_reset_list[0];
  curr_init_value = be_shift_f[0];
  Do_segment_integ(res_vec_integ, inside_vec, init_conf, end_conf,
                   stp, curr_init_value, -1);   

  // Segment between two reset points
  if (FORWARD == 1){
    // Integrate every segment between two adjacent reset points backward
    for (m = 0; m < (N_reset - 1); m = m + 1){                                     
      init_conf = curr_reset_list[m] + 1;                                          
      end_conf = curr_reset_list[m+1];                                           
      curr_init_value = be_shift_f[m+1];
      Do_segment_integ(res_vec_integ, inside_vec, init_conf, end_conf,             
                       stp, curr_init_value, -1);   
    }
  } else if (FORWARD == 0){
    // Integrate every segment between two adjacent reset points forward
    for (m = 0; m < (N_reset - 1); m = m + 1){ 
//      fprintf(stderr, "here we use reset %d / %d for force\n", m, N_reset); 
      init_conf = curr_reset_list[m] + 1;                                          
      end_conf = curr_reset_list[m + 1];                                           
      curr_init_value = af_shift_f[m];
      Do_segment_integ(res_vec_integ, inside_vec, init_conf, end_conf,             
                       stp, curr_init_value, 1);                                   
    }
  } else if (FORWARD == 2){
    // Integrate even forward and odd backward
    for (m = 0; m < (N_reset - 1); m = m + 2){
      init_conf = curr_reset_list[m] + 1;
      end_conf = curr_reset_list[m + 1];
      curr_init_value = af_shift_f[m];
      Do_segment_integ(res_vec_integ, inside_vec, init_conf, end_conf,
                       stp, curr_init_value, 1);     
    }
    for (m = 1; m < (N_reset - 1); m = m + 2){                                     
      init_conf = curr_reset_list[m] + 1;                                          
      end_conf = curr_reset_list[m + 1];                                           
      curr_init_value = be_shift_f[m+1];
      Do_segment_integ(res_vec_integ, inside_vec, init_conf, end_conf,             
                       stp, curr_init_value, -1); 
    } 
  } else{
    fprintf(stderr, "Error: Forward/backward/both\n");
    exit(1);
  }

  // Last segment
//  fprintf(stderr, "here we use %d / %d reset for force\n", N_reset - 1, N_reset);
  init_conf = curr_reset_list[N_reset - 1] + 1;
  end_conf = tot_len - 1;
  curr_init_value = af_shift_f[N_reset - 1];
  Do_segment_integ(res_vec_integ, inside_vec, init_conf, end_conf,
                     stp, curr_init_value, 1);    
}

void Do_segment_integ(double *res_vec_integ, double *inside_vec,
                      int init_conf, int end_conf, double stp,
                      double init_value, int direction)
{
  int i;

  if (end_conf == init_conf){
    fprintf(stderr, "special case, only one conf in one segment %d\n", 
            init_conf);
    res_vec_integ[init_conf] = init_value; 
  } else if (end_conf > init_conf){
    if (direction == 1){
      Do_integ_simpson(res_vec_integ, inside_vec, init_conf, end_conf,
                       stp, init_value); 
    } else if (direction == -1){
      Do_integ_reverse_simpson(res_vec_integ, inside_vec, init_conf, end_conf,
                       stp, init_value);
    }     
  } else{
    fprintf(stderr, "Error: integ sgement from %d to %d\n", init_conf, end_conf);
  }
}

void Do_integ_simpson(double *res_vec_integ, double *inside_vec,
                     int init_conf, int end_conf, double stp,
                     double init_value)
{
  int m;

  res_vec_integ[init_conf] = init_value;
  res_vec_integ[init_conf + 1] = res_vec_integ[init_conf] +
      0.5 * stp * (inside_vec[init_conf] + inside_vec[init_conf + 1]);

  for (m = (init_conf + 2); m <= end_conf; m++){
    res_vec_integ[m] = res_vec_integ[m - 2] + (stp / 3.0) *
        (inside_vec[m - 2] + 4.0 * inside_vec[m - 1] + inside_vec[m]);
  }
}

void Do_integ_reverse_simpson(double *res_vec_integ, double *inside_vec,
                              int init_conf, int end_conf, double stp,
                              double init_value)
{
  int m;

  res_vec_integ[end_conf] = init_value;
  res_vec_integ[end_conf - 1] = res_vec_integ[end_conf] - 
          0.5 * stp * (inside_vec[end_conf] + inside_vec[end_conf - 1]);

  for (m = (end_conf - 2); m >= init_conf; m--){
    res_vec_integ[m] = res_vec_integ[m + 2] - (stp / 3.0) *
         (inside_vec[m + 2] + 4.0 * inside_vec[m + 1] + inside_vec[m]);
  }
}

void Do_simpson_piecewise(double *res_vec_integ, double *inside_vec,
                          int init_conf, int end_conf, double stp,
                          double init_value,
                          int *reset_list, int num_reset,
                          double *portion_reset, double *before_work,
                          double *after_work, double *x_dot, 
                          double *force_comp, int curr_dim, int df,
                          double *U_ab_b_dot,
                          double *be_shift_f, double *af_shift_f,
                          double *before_force, double *after_force)
{
  int m, j;
  int add_conf;
  int reset_id;
  double reset_x_dot, after_x_dot;

  res_vec_integ[init_conf] = init_value;

  reset_id = Whe_in_list(init_conf, reset_list, num_reset);

  if (reset_id == -1){
    res_vec_integ[init_conf + 1] = res_vec_integ[init_conf] +
        0.5 * stp * (inside_vec[init_conf] + inside_vec[init_conf + 1]);
    add_conf = 2;
  } else{
    Do_reset_interval_integ(res_vec_integ, stp, portion_reset,
        before_work, after_work, x_dot, force_comp, be_shift_f,
        af_shift_f, before_force, after_force, init_conf + 1, reset_id); 

    if (((init_conf+2) <= end_conf)){
      add_conf = 3;
      fprintf(stderr, "reset %d start %d\n", reset_id, init_conf + add_conf);
    }
  }

  // Need to be super careful, whether it is reset point
  // whether it is one of two points after reset point 
  for (m = (init_conf + add_conf); m <= end_conf; m++){
    reset_id = Whe_in_list(m-1, reset_list, num_reset);
    
    if (reset_id == -1){
      res_vec_integ[m] = res_vec_integ[m-2] + (stp / 3.0) *
          (inside_vec[m-2] + 4.0 * inside_vec[m-1] + inside_vec[m]);
    } else{
//      fprintf(stderr, "here we use reset %d / %d for work\n", reset_id, num_reset);   
      Do_reset_interval_integ(res_vec_integ, stp, portion_reset, 
          before_work, after_work, x_dot, force_comp, be_shift_f,
          af_shift_f, before_force, after_force, m, reset_id);

      // Use trapezoid rule for next step, remember to update the
      // conf_id + 1
      if ((m+1) <= end_conf){
        res_vec_integ[m+1] = res_vec_integ[m] + 0.5 * stp *
            (inside_vec[m] + inside_vec[m+1]);
        m++;
      }
    }
  }
}

void Do_reset_interval_integ(double *res_vec_integ,                              
                             double stp,                                         
                             double *portion_reset, double *before_work,         
                             double *after_work, double *x_dot,                  
                             double *force_comp,                                 
                             double *be_shift_f, double *af_shift_f,             
                             double *before_force, double *after_force,          
                             int m, int reset_id)
{
  double reset_x_dot, after_x_dot;

  reset_x_dot = x_dot[m-1] +                                                 
      (x_dot[m] - x_dot[m-1]) * portion_reset[reset_id];                     
  after_x_dot = reset_x_dot;                                                 

  before_force[reset_id] = force_comp[m-1] - be_shift_f[reset_id];           
  after_force[reset_id] = force_comp[m] - af_shift_f[reset_id];              
                                                                                 
  before_work[reset_id] = 0.5 * stp * portion_reset[reset_id] *              
      (force_comp[m-1] * x_dot[m-1] + before_force[reset_id] * reset_x_dot); 
  after_work[reset_id] = 0.5 * stp * (1.0 - portion_reset[reset_id]) *       
      (force_comp[m] * x_dot[m] + after_force[reset_id] * after_x_dot);      
                                                                                 
  res_vec_integ[m] = res_vec_integ[m-1] +                                    
      before_work[reset_id] + after_work[reset_id];
}

int Whe_in_list(int val, int *vector, int len)
{
  int i;

  for (i = 0; i < len; i++){
    if (val == vector[i]){
      return i;
    }
  }

  return -1;
}

void Get_shift_f(double *be_shift_f, double *af_shift_f,
                 double* inside_vec, int N_reset,
                 double* curr_portion_reset,
                 int* curr_reset_list, double t_stp)
{
  int m, conf_id;

  for (m = 0; m < N_reset; m++){
    conf_id = curr_reset_list[m];

    be_shift_f[m] = 0.0 - 0.5 * t_stp * curr_portion_reset[m] *
                 (inside_vec[conf_id] + inside_vec[conf_id + 1]);
    af_shift_f[m] = 0.0 + 0.5 * t_stp * (1.0 -  curr_portion_reset[m]) *
                 (inside_vec[conf_id] + inside_vec[conf_id + 1]);
  }
}


void Cal_portion_init_value(double* portion_init_value, 
                            double* inside_vec, int N_reset,                  
                            double* curr_portion_reset,                       
                            int* curr_reset_list, double t_stp)               
{                                                                                
  int m, conf_id;
                                                                                 
  for (m = 0; m < N_reset; m++){        
    conf_id = curr_reset_list[m];                                        
    portion_init_value[m] = 0.0 + (inside_vec[conf_id] +              
           inside_vec[conf_id + 1]) * 0.5 *                           
           ((1.0 - curr_portion_reset[m]) * t_stp);                              
  }                                                                              
}                                                                                
  

void New_Count_row_mat_cross_num(double** mat, int de_row, int de_col,           
                              int* N_cross, double value,                        
                              int max_num_reset, double** portion_reset,         
                              int **cross_list)                                  
{                                                                                
  int m, n;                                                                      
                                                                                 
  Init_int_vec(N_cross, de_row, 0);                                              
  Init_mat(portion_reset, de_row, max_num_reset, -100.0);                          
                                                
  for (m = 0; m < de_row; m++){
    Init_int_vec(cross_list[m], max_num_reset, -1);
                                                  
    for (n = 0; n < (de_col - 1); n++){                                          
      if (((mat[m][n] > value) && (mat[m][n+1] < value))                         
          || ((mat[m][n] < value) && (mat[m][n+1] > value))){                    
        cross_list[m][N_cross[m]] = n;                                           
        portion_reset[m][N_cross[m]] = (value - mat[m][n]) /                     
                                      (mat[m][n+1] - mat[m][n]);                 
        N_cross[m] = N_cross[m] + 1;                                             
      } else if (mat[m][n] == value){                                            
        fprintf(stderr, "right crossing at point\n");                            
        cross_list[m][N_cross[m]] = n;
        portion_reset[m][N_cross[m]] = 0.0;
        N_cross[m] = N_cross[m] + 1;
      }                                                                          
    }                                                                            
  }                                                                              
}

double** Malloc_mat(int de1, int de2)                                            
{                                                                                
  int m;                                                                         
  double** mat;                                                                  
                                                                                 
  mat = malloc(de1*sizeof(double*));                                             
                                                                                 
  if (mat == NULL){                                                              
    printf("Cannot malloc **\n");                                                
  } else{                                                                        
    for (m = 0; m < de1; m++){                                                   
      mat[m] = malloc(de2*sizeof(double));                                       
      if (mat[m] == NULL){                                                       
        printf("cannot malloc * %d\n", m);                                       
      }                                                                          
    }                                                                            
  }                                                                              
                                                                                 
  return mat;                                                                    
} 

void Free_mat(double** mat, int ind)                                             
{                                                                                
  int m;                                                                         
                                                                                 
  for (m = 0; m < ind; m++){                                                     
    free(mat[m]);                                                                
  }                                                                              
  free(mat);                                                                     
}  

int** Malloc_int_mat(int de1, int de2)                                           
{                                                                                
  int m;                                                                         
  int** mat;                                                                     
                                                                                 
  mat = malloc(de1*sizeof(int*));                                                
                                                                                 
  if (mat == NULL){                                                              
    fprintf(stderr, "Cannot malloc int **\n");                                   
  } else{                                                                        
    for (m = 0; m < de1; m++){                                                   
      mat[m] = malloc(de2*sizeof(int));                                          
      if (mat[m] == NULL){                                                       
        fprintf(stderr, "Cannot malloc int* %d\n", de2);                         
      }                                                                          
    }                                                                            
  }                                                                              
                                                                                 
  return mat;                                                                    
}

void Free_int_mat(int** mat, int de)                                             
{                                                                                
  int m;                                                                         
                                                                                 
  for (m = 0; m < de; m++){                                                      
    free(mat[m]);                                                                
  }                                                                              
                                                                                 
  free(mat);                                                                     
}  

void Init_mat(double** mat, int de1, int de2, double value)                      
{                                                                                
  int m, n;                                                                      
                                                                                 
  for (m = 0; m < de1; m++){                                                     
    for (n = 0; n < de2; n++){                                                   
      mat[m][n] = value;                                                         
    }                                                                            
  }                                                                              
} 

void Init_vec(double* vec, int len, double value)                                
{                                                                                
  int m;                                                                         
                                                                                 
  for (m = 0; m < len; m++){                                                     
    vec[m] = value;                                                              
  }                                                                              
}                                                                                
                                                                                 
void Init_int_vec(int* vec, int len, int value)                                  
{                                                                                
  int m;                                                                         
                                                                                 
  for (m = 0; m < len; m++){                                                     
    vec[m] = value;                                                              
  }                                                                              
}                                                                                
    
