import os
optimizers = ['RF']#'GP'] #, 'RF', 'GBRT', 'ET']
loss_funcs = ['bloss62' ,'ref15']

for opti in optimizers:
    for loss in loss_funcs:
        job_dir = '%s/.jobs' %os.getcwd()
	if not os.path.exists(job_dir):
		os.makedirs(job_dir)

        job_file = opti+'_'+loss+'.job'
        with open(job_file, 'w') as h:
            h.write('#!/bin/bash\n\
#SBATCH --time=48:00:00\n\
#SBATCH --mem-per-cpu=6G\n\
#SBATCH --mail-type=FAIL\n\
#SBATCH --partition=clara-job\n\
#SBATCH --mail-user=jannis.deriz@gmail.com\n\
#SBATCH -n 1\n\
#SBATCH -c 64\n\
#SBATCH --job-name='+opti+'_'+loss+'\n\
module load Anaconda3\n\
source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh\n\
conda activate scfxn\n\
python -u bench.py '+opti+' '+loss+' > py_out/'+opti+'_'+loss+'.out\n\
                    ')
        
        os.system("sbatch %s" %job_file)



