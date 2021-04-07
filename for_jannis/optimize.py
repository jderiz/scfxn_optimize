#!/usr/bin/env python

import argparse
import sys
import random
import time
import copy
import os

def main():
	parser = argparse.ArgumentParser(description='Run Nelder-Mead Optimization')
	parser.add_argument("-c","--checkpoint", type=str, help="checkpoint file")
	parser.add_argument("-s","--structcounter", type=int, default=1, help="struct counter")
	parser.add_argument("-f","--scorefunction", type=str, default="wannier_static_weights.wts", help="score function")
	parser.add_argument("-j","--joblist", type=str, default="alljobs", help="list of rosetta jobs to run")
	parser.add_argument("-r","--randomize", type=float, help="value to set randomization size")
	parser.add_argument("-u","--usehbscore", action="store_true", help="use the hbdel score as part of the target function")
	args = parser.parse_args()
	
	# constants
	consts = {"ALPHA":1.0, \
		"BETA":0.5, \
		"GAMMA":2.0, \
		"STRUCTCOUNTER":args.structcounter, \
		"DISTR_MIN":0.002, \
		"scorefunction":args.scorefunction, \
		"joblist":args.joblist, \
		"usehbscore":args.usehbscore, \
		"wt_rr_int":4.0, \
		"wt_rr_core":4.0, \
		"wt_decoy":2.0, \
		"wt_min_decoy":2.0, \
		"wt_docking":1.0, \
		"wt_seqrecov":2.0, \
		"wt_distr":10.0, \
		"wt_liqsim":0.1, \
		"wt_hbdel":0.01, \
		"wt_lig":10.0, \
		"maxiter":1000, \
		"ftol":1e-3, \
		"sttime":time.time() \
		}

	EVALONLY = 0

	# load arrays
	guess, scale, names, links = load_values_from_file()

	if args.randomize is not None:
		for i in range(0, len(guess)):
			offset = args.randomize*(random.uniform(0,2)-1)*scale[i]
			guess[i] = guess[i] + offset

	if args.checkpoint is None:
		opt, cost = MinimizeND(names, guess, scale, links, consts)

def load_values_from_file():
	f = open("vertices","r")
	lines = f.readlines()
	f.close()

	g = []
	s = []
	n = []
	lk = [] #array for linked vertices
		
	for l in lines:
		if l[0:4] == "VERT":
			split_line = l.split()
			n.append(split_line[1])	
			g.append(float(split_line[2]))	
			s.append(float(split_line[3]))
		if l[0:4] == "VLNK": #vertices that have linked parameters
			split_line = l.split()
			n.append(split_line[1])	
			g.append(float(split_line[2]))	
			s.append(float(split_line[3]))
			lk.append(split_line[1]) # same as name
	

	return g,s,n,lk
		

def MinimizeND(names, guess, scale, links, consts):
	p = ConstructVertices(guess, scale)
	y = EvaluateVertices(names, links, p, consts)
	return Amoeba(p, y, consts, names, links)
	
def ConstructVertices(guess, scale):
	p = []
	p.append(guess)

	for i in range(0,len(guess)):
		g = copy.copy(guess)
		g[i] = g[i] + scale[i]
		p.append(g)

	return p

def PrintVertexNamesValues(names, values):
	for i in range(0, len(names)):
		print names[i]+":"+str(values[i])+", ",
	print
	

def EvaluateVertices(names, links, p, consts):
	print "EVALUATE INITIAL VERTICES"

	y = []
	for i in range(0,len(p)):
		y.append(func(names, links, p[i], consts))

	return y

def get_charges_for_group(weight, item_charge, avg_net_charge, avg_sum_charge):
	return float(weight) * (float(item_charge) - avg_sum_charge) + avg_net_charge

	

def func(names, links, values, consts):
	result = ""
	resfields = []
	run = 1

	while(run):

		if not os.path.exists("opt_"+str(consts["STRUCTCOUNTER"])):
			os.makedirs("opt_"+str(consts["STRUCTCOUNTER"]))
			os.makedirs("opt_"+str(consts["STRUCTCOUNTER"])+"/ligdock")

		# open file with the full flag name to modify parameter in Rosetta
		fn = open("flag_names", "r")
		flag_names = fn.readlines()
		fn.close()

		flag_name_dict = {}
		for l in flag_names:
			if l[0:4] == "FLAG":
				split_line = l.split()
				flag_name_dict[split_line[1]] = " ".join(split_line[2:])

		# open link file
		lk = open("links","r")
		link_lines = lk.readlines()
		lk.close()

		# must be a list, not dict
		# each line in the links file corresponds to a vertex name
		# multiple lines can share the same name, as the parameters are
		# modified by the same vertex weight
		link_params_list = []
		for l in link_lines:
			if l[0:4] == "LINK":
				link_params_list.append(l.split())
			
		# Loop over all vertices, and print to flags file
		flags = open("flags_"+str(consts["STRUCTCOUNTER"]),"w")
		for n in range(0, len(names)):
			# loop over vertices w/o links
			if names[n] not in links:
				if flag_name_dict[names[n]][-1] == ":": # if name ends in colon, don't put a space between the flag and value
					flags.write("-{0}{1:.5f}\n".format(flag_name_dict[names[n]], float(values[n])))
				else:
					flags.write("-{0} {1:.5f}\n".format(flag_name_dict[names[n]], float(values[n])))
				
			# loop over vertices with links
			elif names[n] in links:
				for l in link_params_list:
					if names[n] == l[1]:
						atoms_in_group = len(l[3:])/2

						avg_net_charge = float(l[2]) / atoms_in_group
						avg_sum_charge = sum(map(float, l[4::2])) / atoms_in_group

						for ll in range(3, len(l), 2):
							scaled_value = get_charges_for_group(values[n], l[ll+1], avg_net_charge, avg_sum_charge)
							if flag_name_dict[l[ll]][-1] == ":": # if name ends in colon, don't put a space between the flag and value
								flags.write("-{0}{1:.5f}\n".format(flag_name_dict[l[ll]], scaled_value))
							else:
								flags.write("-{0} {1:.5f}\n".format(flag_name_dict[l[ll]], scaled_value))

		flags.close()

		# Run Rosetta Jobs
		print "vertex "+str(consts["STRUCTCOUNTER"])+" - ",
		PrintVertexNamesValues(names, values)
		print "Run Rosetta Jobs"
		### SLURM ###
		run_command = "sbatch --wait --export=JOB_LIST="+consts["joblist"]+",NUM="+str(consts["STRUCTCOUNTER"])+",WEIGHTS="+consts["scorefunction"]+" run_rosetta_array.slurm"
		### LOCAL ###
		#run_command = "parallel -j 15 --workdir . :::: "+consts["joblist"]+" ::: "+str(consts["STRUCTCOUNTER"])+" ::: "+consts["scorefunction"]+" &> /dev/null"

		print run_command
		os.system(run_command)

		# Clean-up jobs
		print "Clean-up Jobs"
		run_cleanup = "./run_cleanup.sh opt_"+str(consts["STRUCTCOUNTER"])+" &> /dev/null"
		print run_cleanup
		os.system(run_cleanup)

		# hbond donor and acceptor score
		if (consts["usehbscore"]):			
			Ndon, Nacc = 0,0
			hbdon, hbdon2, hbacc, hbacc2 = 0,0,0,0

			for n in range(0, len(names)):
				if names[n][0:5] == "hbdon":
					hbdon += values[n]
					hbdon2 += values[n]*values[n]
					Ndon += 1
				if names[n][0:5] == "hbacc":
					hbacc += values[n]
					hbacc2 += values[n]*values[n]
					Nacc += 1

			hbdon_std = ( Ndon*hbdon2 - hbdon*hbdon )**0.5 / Ndon;
			hbacc_std = ( Nacc*hbacc2 - hbacc*hbacc )**0.5 / Nacc;

			score_hbdel = hbdon_std + hbacc_std

		# initialize score terms

		score_rr_int, score_rr_core, score_decoy, score_min_decoy, score_docking, score_seqrecov, score_distr, score_lig = 0,0,0,0,0,0,0,0
		validate_rr_int, validate_rr_core, validate_decoy, validate_min_decoy, validate_docking, validate_seqrecov, validate_distr, validate_lig = -1,-1,-1,-1,-1,-1,-1,-1

		# Open results
		print "Evaluate Results"
		opt_folder = "opt_"+str(consts["STRUCTCOUNTER"])+"/"

		#TURN ON LATER
		#resin = open(opt_folder+"rr_interface_result")
		#result = resin.readlines()
		#resin.close()
		#resfields = result[0].rstrip().split(",")
		#print "...rr_interface_result "+result[0].rstrip()
		#score_rr_int = -1.0*float(resfields[0])
		#validate_rr_int = float(resfields[2])

		#TURN ON LATER
		#resin = open(opt_folder+"rr_core_result")
		#result = resin.readlines()
		#resin.close()
		#resfields = result[0].rstrip().split(",")
		#print "...rr_core_result "+result[0].rstrip()
		#score_rr_core = -1.0*float(resfields[0])
		#validate_rr_core = float(resfields[2])

		#TURN ON LATER
		#resin = open(opt_folder+"score_decoy_result")
		#result = resin.readlines()
		#resin.close()
		#resfields = result[0].rstrip().split()
		#print "...score_decoy_result "+result[0].rstrip()
		#score_decoy = -1.0*float(resfields[0])
		#validate_decoy = float(resfields[2])

		#TURN ON LATER
		#resin = open(opt_folder+"score_min_result")
		#result = resin.readlines()
		#resin.close()
		#resfields = result[0].rstrip().split()
		#print "...score_min_result "+result[0].rstrip()
		#score_min_decoy = -1.0*float(resfields[0])
		#validate_min_decoy = float(resfields[2])

		#TURN ON LATER
		#resin = open(opt_folder+"score_ddg_result")
		#result = resin.readlines()
		#resin.close()
		#resfields = result[0].rstrip().split()
		#print "...score_ddg_result "+result[0].rstrip()
		#score_docking = -1.0*float(resfields[0])
		#validate_docking = float(resfields[2])

		#TURN ON LATER
		#resin = open(opt_folder+"seqrecov_result")
		#result = resin.readlines()
		#resin.close()
		#resfields = result[0].rstrip().split()
		#print "...seqrecov_result "+result[0].rstrip()
		#score_seqrecov = float(resfields[0])
		#validate_seqrecov = float(resfields[1])

		#TURN ON LATER
		#resin = open(opt_folder+"distr_result")
		#result = resin.readlines()
		#resin.close()
		#resfields = result[0].rstrip().split()
		#print "...distr_result "+result[0].rstrip()
		#score_distr = float(resfields[0])
		#validate_distr = float(resfields[1])

		resin = open(opt_folder+"ligdock_result")
		result = resin.readlines()
		resin.close()
		resfields = result[0].rstrip().split()
		print "...ligdock_result "+result[0].rstrip()
		score_lig = float(resfields[0])
		validate_lig = 1 # CHANGE THIS LATER


		#TURN ON LATER
		#if (score_distr < consts["DISTR_MIN"]):
		#	score_distr = consts["DISTR_MIN"]	

		elapsed_time = time.time() - consts["sttime"]

		run = False
		#TURN ON LATER
		#if( validate_results(validate_rr_int, validate_rr_core, validate_decoy, validate_min_decoy, validate_docking, validate_seqrecov, validate_distr) == 0 ):
		#	print "[{0:4.0f} : {1:6.0f} s]  Error! rerunning!  ({2}, {3}, {4}, {5}, {6}, {7}, {8})".format(consts["STRUCTCOUNTER"], elapsed_time, validate_rr_int, validate_rr_core, validate_decoy, validate_min_decoy, validate_docking, validate_seqrecov, validate_distr)
		#	run = True # rerun rosetta jobs

	# end while loop

	# report scores
	score_tot = ( \
		consts["wt_rr_int"] * score_rr_int + \
		consts["wt_rr_core"] * score_rr_core + \
		consts["wt_decoy"] * score_decoy + \
		consts["wt_min_decoy"] * score_min_decoy + \
		consts["wt_docking"] * score_docking + \
		consts["wt_seqrecov"] * score_seqrecov + \
		consts["wt_distr"] * score_distr + \
		consts["wt_lig"] * score_lig)

	if (consts["usehbscore"]):
		score_tot += consts["wt_hbdel"]* score_hbdel

	if (consts["usehbscore"]):
		print "[{0:4.0f} : {1:6.0f} s] TOTAL: {2:8.4f} rr_int={3:8.4f} rr_core={4:8.4f} decoy={5:8.4f} decoy_min={6:8.4f} docking={7:8.4f} seqrec={8:8.4f} distr={9:8.4f} hbdel={10:8.4f} score_lig={11:8.4f}\n".format(consts["STRUCTCOUNTER"], elapsed_time, score_tot, score_rr_int, score_rr_core, score_decoy, score_min_decoy, score_docking, score_seqrecov, score_distr, score_hbdel, score_lig)
	else:
		print "[{0:4.0f} : {1:6.0f} s] TOTAL: {2:8.4f} rr_int={3:8.4f} rr_core={4:8.4f} decoy={5:8.4f} decoy_min={6:8.4f} docking={7:8.4f} seqrec={8:8.4f} distr={9:8.4f} score_lig={10:8.4f}\n".format(consts["STRUCTCOUNTER"], elapsed_time, score_tot, score_rr_int, score_rr_core, score_decoy, score_min_decoy, score_docking, score_seqrecov, score_distr, score_lig)

	consts["STRUCTCOUNTER"] += 1

	return score_tot

def Amoeba(p, y, consts, names, links):

	print "RUN NELDER-MEAD OPTIMIZATION"

	n = len(p)
	iter = 0

	# To control the recalculation of centroid
	recalc = True
	ihi_o = ""

	while(1):
		ilo, inhi, ihi = _FindMarkers(y)

		rtol = 2*abs(y[ihi]-y[ilo]) / (abs(y[ihi])+abs(y[ilo])+1e-16)
		print "Check if gradient has converged, rtol {0:.5f}".format(rtol)

		# Stopping conditions
		iter += 1
		if (rtol < consts["ftol"]):
			"Gradient has converged... Stopping."
			break
		if (iter>consts["maxiter"]):
			"Max iteration has been reached... Stopping."
			break

		# write checkpoint to disk
		f = open("amoeba.chkpt","w")
		f.write(str(n-1)+"\n")
		for i in range(0,n):
			f.write(str(y[i])+" ")
			for j in range(0,n-1):
				f.write(str(p[i][j])+" ")
			f.write("\n")
		f.close()

		# Determine the Centroid
		if (recalc):
			pbar = _CalcCentroid(p, ihi)
		else:
			pbar = _AdjustCentroid(pbar, p, ihi_o, ihi)

		# Reset the re-calculation flag, and remember the current highest
		recalc = False
	
		# Determine the reflection point, evaluate its value
		pr = _CalcReflection(pbar, p[ihi], consts["ALPHA"])
		ypr = func(names, links, pr, consts)

		# if it gives a better value than best point, try an
		# additional extrapolation by a factor gamma, accept best
		if (ypr < y[ilo]):
			print "Reflected vertex value lower than current best, Gamma Move, calculating value of expanded vertex"
			pe = _CalcReflection(pbar, pr, -1.0*consts["GAMMA"])
			ype = func(names, links, pe, consts)

			if ype < y[ilo]:
				print "Gamma vertex value lower than current best, replacing highest point"
				p[ihi] = pe
				y[ihi] = ype
			else:
				print "Gamma vertex value worse than current best, reflected point chosen"
				p[ihi] = pr
				y[ihi] = ypr
		
		# if reflected point is worse than 2nd highest
		elif (ypr >= y[inhi]):
			print "Reflected vertex value higher than 2nd highest"
			# if it is better than the highest, replace it
			if (ypr < y[ihi]):
				print "Reflected vertex value lower than highest, replacing highest point"
				p[ihi] = pr
				y[ihi] = ypr

			# look for an intermediate lower point by performing a 
			# contraction of the simplex along one dimension
			print "Beta Move, calculating value of contracted simplex"
			pc = _CalcReflection(pbar, p[ihi], -1.0*consts["BETA"])
			ypc = func(names, links, pc, consts)

			# if contraction gives an improvement, accept it
			if (ypc < y[ihi]):
				print "Contraction gives improvement, replacing highest point"
				p[ihi] = pc
				y[ihi] = ypc

			# otherwise can't seem to remove high point
			# so contract around lowest (best) point
			else:
				print "Beta move cannot remove highest point, contracting around lowest point"
				for i in range(0, n):
					if i != ilo:
						"Recomputing vertex "+str(i)
						p[i] = _CalcReflection(p[ilo], p[i], -1.0*consts["BETA"])
						y[i] = func(names, links, p[i], consts)
				recalc = True
		# if reflected point is in-between lowest and 2nd highest
		else:
			print "Reflected vertex value between lowest and 2nd highest, replacing highest point"
			p[ihi] = pr
			y[ihi] = ypr

		# Remember the replacing position and its value
		ihi_o = ihi

	return p[ilo], y[ilo]

# Determine the reflection point
def _CalcReflection(p1, p2, scale):
	print "Determining reflection of Simplex, calculating value of reflected vertex"
	n = len(p1)

	pr = []
	for j in range(0,n):
		pr.append(p1[j] + scale*(p1[j] - p2[j]))

	return pr
		

# Determine the centroid (except the highest point)
def _CalcCentroid(p, ihi):
	n = len(p) - 1

	pbar = []
	for j in range(0,n):
		pbar.append(0)
		for i in range(0, n+1):
			if (i != ihi):
				pbar[j] += p[i][j]
		pbar[j] = pbar[j] / n

	return pbar

# Adjust the centroid only
def _AdjustCentroid(pbar, p, ihi_o, ihi):
	n = len(pbar)

	if (ihi_o != ihi):
		for j in range(0, n):
			pbar[j] = pbar[j] + ((p[ihi_o][j] - p[ihi][j]) / n)

	return pbar


def _FindMarkers(y):
	# ilo - lowest vertex value, ihi - highest vertex value, inhi - 2nd highest vertex value
	n = len(y)-1
	ilo = 0

	if (y[0]>y[1]):
		ihi=0
		inhi=1
	else:
		ihi=1
		inhi=0

	for i in range(0, len(y)):
		if (y[i] < y[ilo]):
			ilo = i
		if (y[i] > y[ihi]):
			inhi = ihi
			ihi = i
		elif (y[i] > y[inhi] and ihi != i):
			inhi = i	

	return ilo, inhi, ihi
	
	
def validate_results(validate_rr_int, validate_rr_core, validate_decoy, validate_min_decoy, validate_docking, validate_seqrecov, validate_distr):
	print "Validate Results (Did Rosetta jobs run to completion) --",

	#if (validate_rr_int!=-1 and validate_rr_int < 3390): 
	#	print "FAILED!"
	#	return 0  # allow a bit of fluff (sql sometimes flakes)
	#if (validate_rr_core!=-1 and validate_rr_core < 18820): 
	#	print "FAILED!"
	#	return 0  # allow a bit of fluff (sql sometimes flakes)
	#if (validate_decoy!=-1 and validate_decoy < 100112):
	#	print "FAILED!"
	#	return 0
	#if (validate_min_decoy!=-1 and validate_min_decoy < 5340): 
	#	print "FAILED!"
	#	return 0
	#if (validate_docking!=-1 and validate_docking < 8997):
	#	print "FAILED!"
	#	return 0
	#if (validate_seqrecov!=-1 and validate_seqrecov < 4606):
	#	print "FAILED!"
	#	return 0
	#if (validate_distr!=-1 and validate_distr < 56):
	#	print "FAILED!"
	#	return 0

	print "SUCCESS!"
	return 1

if __name__ == "__main__":
	main()
