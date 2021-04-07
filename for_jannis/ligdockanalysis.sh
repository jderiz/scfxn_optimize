dirname=$1

for f in  ${dirname}/*/ 
do
echo $f
folder=`echo $f | awk -F"/" '{print $(NF-1)}'`
protein=`echo $folder | awk -F"_" '{print $1}'`
lig=`echo $folder | awk -F"_" '{print $2}'`
name=${protein}_${lig}
cd $f
python ../../../scripts/lig_rmsd_to_native.py -n /home/muellebk/Projects/orbitals/wannier/shannon_ligands/platinum_dataset/natives/${protein}_${lig}_* -l $lig -o ${name}.png score.sc
python ../../../scripts/lig_merge.py score.sc ${name}.txt ${name}_score_rmsd.sc
cd ../../../
done
python /home/smithst/entropic_penalties/pnear_ligands/analysis.full_docking.py -scorefiles ${dirname}/*/*_score_rmsd.sc -summary ${dirname}/all_pnear.txt
