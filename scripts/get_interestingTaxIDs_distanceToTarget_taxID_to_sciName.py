#!/usr/bin/env python

# This is a script that will substitue a step of a perSVade function that needs to be run on perSVade


#### ENVIRONMENT  ####
import sys
import os
import pickle
from ete3 import NCBITaxa

# get the args
ancestor_taxID, target_taxID, outfile = sys.argv[1:]

ancestor_taxID = int(ancestor_taxID)
target_taxID = int(target_taxID)

######################

####### FUNCTIONS #######

def save_object(obj, filename):
    
    """ This is for saving python objects """

    filename_tmp = "%s.tmp"%filename
    remove_file(filename_tmp)
    
    with open(filename_tmp, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

    os.rename(filename_tmp, filename)


def remove_file(f):

    if os.path.isfile(f): 

        try: run_cmd("rm %s > /dev/null 2>&1"%f)
        except: pass

#########################

# get the ncbi
ncbi = NCBITaxa()

# get the tree
tree = ncbi.get_descendant_taxa(ancestor_taxID, collapse_subspecies=False, return_tree=True, intermediate_nodes=True)
#print("You can visualize the tree in http://etetoolkit.org/treeview/: \n\n%s\n\n"%(tree.write()))

# define the interesting taxIDs and the taxID_to_distanceToTarget
if type(tree)==list:
    interesting_taxIDs = set(tree)
    taxID_to_distanceToTarget = {tree[0] : 0}

else:

    # define interesting taxIDs (the leafs and the species names that may be intermediate)
    interesting_taxIDs = {int(x) for x in set(tree.get_leaf_names()).union({n.name for n in tree.traverse() if n.rank=="species"})}.difference({target_taxID})

    # map the distance between each leave and the target
    taxID_to_distanceToTarget = {taxID : tree.get_distance(str(target_taxID), str(taxID)) for taxID in interesting_taxIDs}

# map names
taxID_to_sciName = ncbi.get_taxid_translator(interesting_taxIDs)


# save a tuple that has all the interesting objects
final_object = (interesting_taxIDs, taxID_to_distanceToTarget, taxID_to_sciName)
save_object(final_object, outfile)
