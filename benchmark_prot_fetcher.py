import os
import random 
from pyrosetta.rosetta import pose_from_rcsb, pose_from_pdb

global relaxed, crystals
relaxed = {}
crystals = {}

for pdb in os.listdir("benchmark"):
    if pdb.endswith(".pdb"):
        #  store actual Pose() object
        name = pdb.split(".")[0]
        relaxed.update({name: pose_from_pdb("benchmark/" + pdb)})
        crystals.update({name: pose_from_rcsb(name)})


def get_crystals():
    """
    returns the crystal structure of the benchmark proteins as a dict
    """

    return crystals

def get_crystal_by_pdb(name):
    return crystals[name]

def get_all():

    """Function to return the relaxed benchmark set
    @param param:  Description
    @type  param:  Type

    @return:  dict containing all pdbs by name
    @rtype :  dict
    """

    return relaxed


def get_random():

    """Function to return a random relaxed protein from the benchmark set.
    @return:  random pose object from becnhmark
    @rtype :  Pose

    """
    return random.choice(relaxed)

