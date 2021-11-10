import os
for prot in ['1CKK', '1D5W', '1K9P', '1Q0N', '1QUK', '1TDE', '2J5X', '2LAO', '4AKE', '6Q21']:

    job_file = prot+".job"
    with open(job_file, 'w') as h:
        h.write('#!/bin/bash\n\
    #SBATCH --time=48:00:00\n\
    #SBATCH --mem-per-cpu=6G\n\
    #SBATCH --mail-type=FAIL\n\
    #SBATCH --partition=clara-job\n\
    #SBATCH --mail-user=jannis.deriz@gmail.com\n\
    #SBATCH -n 1\n\
    #SBATCH -c 20\n\
    #SBATCH --job-name=opt_rel\n\
    module load Anaconda3\n\
    source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh\n\
    conda activate scfxn\n\
    python bench.py  -config ../configs/relax_'+prot+'_best.pkl -pdb '+prot+'  -evals 20  -id test_config_'+prot+'\n\
    ')
    os.system("sbatch %s" % job_file)
