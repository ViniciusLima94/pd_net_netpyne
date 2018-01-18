#!/bin/bash 
#SBATCH --job-name=pd
#SBATCH -A csd403
#SBATCH -t 3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH -o /home/salvadord/pd_net_netpyne/pd.run
#SBATCH -e /home/salvadord/pd_net_netpyne/pd.err
#SBATCH --mail-user=salvadordura@gmail.com
#SBATCH --mail-type=end


source ~/.bashrc
cd /home/salvadord/pd_net_netpyne
ibrun -np 48 nrniv -python -mpi PDnet_netpyne_izhi.py
wait
                            