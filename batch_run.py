import os
optimizers = ['GP', 'RF', 'GBRT', 'ET']
loss_funcs = ['SCFXN', 'BLOSS' ,'REF15']

for opti in optimizers:
    for loss in loss_funcs:
        job_file = os.path.join(os.getcwd(), '/jobs/'+opti+' '+loss+'.job') 

        with open(job_file) as h:
            h.writelines('\
                        #!/bin/bash\n\
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
                        python -u bench.py '+opti+' '+loss+' > py_out/'+opti+' '+loss+'.out\n\
                    ')
        
        os.system("sbatch %s" %job_file)



