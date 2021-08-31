import os
for prot in ['1CKK', '1D5W', '1K9P', '1Q0N', '1QUK', '1TDE', '2J5X', '2LAO', '4AKE', '6Q21']:
    os.system('python bench.py  -config default -pdb '+prot+'  -evals 20  -id fastRelax_default'+prot)
