import os
import csv
import cPickle as pickle

here = os.path.dirname(os.path.realpath(__file__))

smiles_lookup = dict( (y, x) for (x, y) in (z.split() for z in open(os.path.join(here, "chembl_20.smi")) if len(z.split())==2))
