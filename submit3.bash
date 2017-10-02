#!/bin/bash
#set a job name  
#SBATCH --job-name=run1
#################  
#a file for job output, you can check job progress
#SBATCH --output=run1.out
#################
# a file for errors from the job 
#SBATCH --error=run1.err
#################
#time you think you need; default is one day
#in minutes in this case, hh:mm:ss
#SBATCH --time=23:59:59
#################
#number of tasks you are requesting 
#SBATCH -n 120
#SBATCH --exclusive
#################
#partition to use
#SBATCH --partition=ser-par-10g-2
#################
#number of nodes to distribute n tasks across 
#SBATCH -N 3
#################

#SB ATCH --nodes=3
#SB ATCH --ntasks-per-node=40

work=/home/li.yi3/122916_msnet

cd $work

#mpirun -np 1 R --no-save < test_snow.R

mpirun -np 1 -oversubscribe R --no-save < msnet_simu3.R
#R CMD BATCH ./snow_test.R
