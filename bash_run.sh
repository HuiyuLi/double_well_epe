#!/bin/bash

dir=`pwd`
stp=0.0025
n_conf=7999
traj_id=$1
num_ho=15
wsq_init=10
wsq_final=181
coup_const=1.8

if [ -d ${dir}/result/force ]; then
  rm -r ${dir}/result/force
fi
if [ -d ${dir}/result/work ]; then
  rm -r ${dir}/result/work
fi
mkdir ${dir}/result/force
mkdir ${dir}/result/work

list_para="${stp} ${n_conf} ${traj_id} ${num_ho} ${wsq_init} ${wsq_final} ${coup_const} ${dir}"

echo ${list_para}
${dir}/do_dw_work_decomp_ala ${list_para}

if [ -d ${dir}/result/ho_${num_ho}_force_${traj_id} ]; then
  rm -r ${dir}/result/ho_${num_ho}_force_${traj_id}
fi
if [ -d ${dir}/result/ho_${num_ho}_work_${traj_id} ]; then
  rm -r ${dir}/result/ho_${num_ho}_work_${traj_id}
fi
mv ${dir}/result/force ${dir}/result/ho_${num_ho}_force_${traj_id}
mv ${dir}/result/work ${dir}/result/ho_${num_ho}_work_${traj_id}
