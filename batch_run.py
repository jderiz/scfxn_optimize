import os

optimizers = ['RF']  # 'GP'] #, 'RF', 'GBRT', 'ET']
loss_funcs = ['pssm']
xi = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10]
kappa = [0.001, 0.01, 0.1, 1, 10, 100, 1000]

for x, k in zip(xi, kappa):
    job_dir = '%s/.jobs' % os.getcwd()

    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    job_file = x+'_'+k+'.job'
    with open(job_file, 'w') as h:
        h.write('#!/bin/bash\n\
#SBATCH --time=48:00:00\n\
#SBATCH --mem-per-cpu=6G\n\
#SBATCH --mail-type=FAIL\n\
#SBATCH --partition=clara-job\n\
#SBATCH --mail-user=jannis.deriz@gmail.com\n\
#SBATCH -n 1\n\
#SBATCH -c 64\n\
#SBATCH --job-name=k'+k+'\n\
module load Anaconda3\n\
source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh\n\
conda activate scfxn\n\
python bench.py -e RF -l pssm -evals 150 -rpc 7 -id xi_'+x+'_kappa'+k+' -xi '+str(x)+' -kappa '+str(k)+'\n\
                    ')
    os.system("sbatch %s" % job_file)
