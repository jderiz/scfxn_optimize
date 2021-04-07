#!/usr/bin/env python
#written by: Benjamin K. Mueller

import argparse
import pandas as pd
import os
import gzip
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import math

class Residue:
        def __init__(self, resiName):
               self.resi_type = resiName
               self.atom_dictionary = {}

        def add_atom(self, atomLine):
               tmpAtom = Atom(atomLine)
               self.atom_dictionary[tmpAtom.getAtomType()] = tmpAtom

        def get_Atom(self, atom_name):
               return self.atom_dictionary[atom_name]

        def get_Resi_Type(self):
               return self.resi_type

        def number_of_atoms(self):
               return len(self.atom_dictionary)

        def get_all_atoms(self):
               return self.atom_dictionary
               
               
class Atom:
        def __init__(self, atomLine):
               self.x = float(atomLine[30:38].strip())
               self.y = float(atomLine[38:46].strip())
               self.z = float(atomLine[46:54].strip())

               self.resiNum = atomLine[22:26].strip()
               self.resiName = atomLine[17:20].strip()
               self.atomType = atomLine[12:16].strip()
               self.chain = atomLine[21:22].strip()

        
        def getX(self):
               return self.x

        def getY(self):
               return self.y

        def getZ(self):
               return self.z

        def getCoor(self):
               coordinates = [self.x,self.y,self.z]
               return coordinates

        def getResiNum(self):
               return self.resiNum

        def getResiName(self):
               return self.resiName

        def getAtomType(self):
               return self.atomType

        def getElement(self):
               return self.atomType.lstrip('0123456789')[0:1]        

        def getChain(self):
               return self.chain

def import_scorefile(scorefile):
    selection = pd.read_csv(scorefile, delim_whitespace=True, skiprows=1)[["description", "total_score", "interface_delta_X"]]
    return selection

def calcRMSD(lig1, lig2):
        distSum = 0
        ligand_1_atoms = lig1.get_all_atoms()
        ligand_2_atoms = lig2.get_all_atoms()
        for atoms in ligand_1_atoms:
               if (ligand_1_atoms[atoms].getElement() != "H"):
                     coor1 = ligand_1_atoms[atoms].getCoor()
                     coor2 = ligand_2_atoms[atoms].getCoor()
                     dist = ((coor1[0]-coor2[0])**2) + ((coor1[1]-coor2[1])**2) + ((coor1[2]-coor2[2])**2);
                     distSum += dist

        return ((1.0/len(ligand_1_atoms))*distSum)**0.5

def open_file(file_name):
        if (os.path.exists(file_name)):
               if file_name.endswith('.gz'): 
                     fileObj = gzip.GzipFile(file_name, 'rb')
                     lines = fileObj.readlines()
                     fileObj.close()
                     return lines
               else:
                     f = open(file_name, "r")
                     lines = f.readlines()
                     f.close()
                     return lines
        else:
               print "Cannot open file "+file_name+", exiting...\n"
               sys.exit(1)

parser = argparse.ArgumentParser(description='From a Rosetta scorefile determine statistics of ligand pose clusters')
parser.add_argument("scorefile", default="score.sc", help="file name of the rosetta scorefile")
parser.add_argument("-n","--native", required=True, help="file name of the native pdb")
parser.add_argument("-l","--ligand", required=True, help="3 letter code of the ligand of interest")
parser.add_argument("-o","--out_png", default="score_v_rmsd.png", help="name of the png file to be plotted")
args = parser.parse_args()

native_pdb = open_file(args.native)
native_ligand = Residue(args.ligand)
out_text_file = args.out_png.replace(".png","")+".txt"

for a in native_pdb:
	if (a[0:6] == "ATOM  " or a[0:6] == "HETATM") and a[17:20] == args.ligand:
		native_ligand.add_atom(a)

scorefile_table = import_scorefile(args.scorefile)

rmsd_list = []
for index, row in scorefile_table.iterrows():
	pdbfilename = row["description"]+".pdb"
	if os.path.exists(pdbfilename):
		pdbfile_lines = open_file(pdbfilename)
	elif os.path.exists(row["description"]+".pdb.gz"):
		pdbfile_lines = open_file(row["description"]+".pdb.gz")
	else:
		print "Cannot find file "+pdbfilename+" exiting..."
		sys.exit(0)

	model_ligand = Residue(args.ligand)
	for a in pdbfile_lines:
		if (a[0:6] == "ATOM  " or a[0:6] == "HETATM") and a[17:20] == args.ligand:
			atomname = a[12:16].strip().rstrip('0123456789')
			if atomname != "H" and atomname != "X":
				model_ligand.add_atom(a)

	rmsd_list.append(calcRMSD(native_ligand, model_ligand))
scorefile_table["rmsd"] = pd.Series(rmsd_list)

top_total = scorefile_table[ scorefile_table.total_score < scorefile_table.total_score.mean() ]
bottom_total = scorefile_table[ scorefile_table.total_score >= scorefile_table.total_score.mean() ]

x_max = 10
if (math.ceil(scorefile_table.rmsd.max()) < x_max):
	x_max = math.ceil(scorefile_table.rmsd.max())

xlim = (0, x_max)
ylim = (math.floor(scorefile_table.interface_delta_X.min()-2),0)
plt.xlim(xlim)
plt.ylim(ylim)
plt.plot(top_total["rmsd"],top_total["interface_delta_X"],".",color='blue')
plt.plot(bottom_total["rmsd"],bottom_total["interface_delta_X"],".",color='blue')
plt.grid(color='gray', linestyle=':',linewidth=0.5)
plt.ylabel("interface_delta score")
plt.xlabel("RMSD (A)")
#plt.savefig(args.out_png, dpi=100, bbox_inches='tight')
plt.savefig(args.out_png, bbox_inches='tight')

#print scorefile_table.sort_values(by="interface_delta_X").head(n=10)
scorefile_table.sort_values(by="interface_delta_X").to_csv(out_text_file, sep='\t', index=False, float_format="%.2f")
