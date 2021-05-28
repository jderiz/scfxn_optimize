for file in *[^_A].pdb ; do
	# clean
	echo process $file
	~/code/Rosetta/tools/protein_tools/scripts/clean_pdb.py $file A
	# replace orig with cleaned 
	prefix=${file%.pdb}
	mv ${prefix}_A.pdb $file
	# remove artifacts
#	rm ${prefix}_A.pdb
	
done
