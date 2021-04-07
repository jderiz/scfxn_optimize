#!/usr/bin/env python
#Written by: Benjamin K. Mueller

import pandas
from matplotlib import pyplot as plt
import argparse
import math

parser = argparse.ArgumentParser(description='Plots E v RMSD with and without cart-bonded term.')
parser.add_argument("-s","--scorefile", help="Rosetta scorefile")
args = parser.parse_args()

full_table = pandas.read_table(args.scorefile, delim_whitespace=True, skiprows=1)

#print full_table["total_score"]

full_table["score-cart_bond"] = full_table["score"] - full_table["cart_bonded"]
print full_table[["score-cart_bond","score","cart_bonded","rms"]]

plt.scatter(full_table["rms"], full_table["score-cart_bond"])
plt.xlabel("RMSD")
plt.ylabel("Score - cart_bonded")
plt.show()

plt.scatter(full_table["rms"], full_table["score"])
plt.xlabel("RMSD")
plt.ylabel("Score")
plt.show()
