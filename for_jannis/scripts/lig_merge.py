import argparse
import numbers

def main():
	parser = argparse.ArgumentParser(description='Merge ligand rmsd values into rosetta scorefile')
	parser.add_argument("rosetta_score")
	parser.add_argument("lig_rmsd")
	parser.add_argument("output")
	args = parser.parse_args()

	with open(args.rosetta_score, 'r') as rose_sc:
		rose_lines = rose_sc.readlines()

	with open(args.lig_rmsd, 'r') as lig_rmsd:
		lig_lines = lig_rmsd.readlines()

	lig_name_rmsd_dict = {}
	for l in lig_lines:
		splitlines = l.split()
		if splitlines[0] != "description":
			lig_name_rmsd_dict[splitlines[0]] = splitlines[3]
	
	with open(args.output, 'w') as outfile:
		for l in rose_lines:
			splitlines = l.split()
			if len(splitlines) < 2:
				outfile.write(l)
			else:
				try:
					#DATA LINE
					tmp = float(splitlines[1]) 
					
					if splitlines[-1] in lig_name_rmsd_dict:
						for f in splitlines[0:-1]:
							outfile.write(f+"\t")
						outfile.write(lig_name_rmsd_dict[splitlines[-1]]+"\t"+splitlines[-1]+"\n")
				except:
					#DESCRIPTION LINE
					for f in splitlines[0:-1]:
						outfile.write(f+"\t")
					outfile.write("ligand_rms_no_super_X\t"+splitlines[-1]+"\n")
						

if __name__ == '__main__':
	main()
