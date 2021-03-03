import os

optimizers = ['GP'] #, 'RF', 'GBRT', 'ET']
loss_funcs = ['pssm']

xi = [0.0001, 0.001, 0.01, 0.1, 1]
kappa = [0.01, 0.1, 1, 10, 100]

for x, k in zip(xi, kappa):
    job_dir = '%s/.jobs' % os.getcwd()

    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    job_file = str(x)+'_'+str(k)+'.job'
    xs = str(x)
    ks = str(k)

    with open(job_file, 'w') as h:
        h.write('#!/bin/bash\n\
#SBATCH --time=48:00:00\n\
#SBATCH --mem-per-cpu=6G\n\
#SBATCH --mail-type=FAIL\n\
#SBATCH --partition=clara-job\n\
#SBATCH --mail-user=jannis.deriz@gmail.com\n\
#SBATCH -n 1\n\
#SBATCH -c 64\n\
#SBATCH --job-name=k'+ks+'\n\
module load Anaconda3\n\
source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh\n\
conda activate scfxn\n\
python bench.py -e GP -l pssm -evals 150 -rpc 8 -id xi_'+xs+'_kappa'+ks+' -xi '+xs+' -kappa '+ks+'\n\
')
    os.system("sbatch %s" % job_file)
