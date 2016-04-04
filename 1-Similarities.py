import os
import sys
import cPickle as pickle
import itertools
from collections import defaultdict, Counter
import benchlib.chembl as chembl
import benchlib.fingerprint_lib as flib # From benchmarking platform
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.Fraggle import FraggleSim

def readBenchmark(fname):
    for line in open(fname):
        yield line.rstrip().split()

def evaluate_similarity_method(dataset, resultsdir):
    size = 4
    ref_corr = range(size, 0, -1)
    ref_corr_b = range(0, size)

    # Setup results dir
    if not os.path.isdir(resultsdir):
        os.mkdir(resultsdir)

    writer = Writer(resultsdir)
    for i, d in enumerate(get_rdkitmols(dataset)):
        for fpName, fpCalculator in flib.fpdict.iteritems():
            ref_mol = d[0]
            ref_fp = fpCalculator(ref_mol)
            ref_nbonds = ref_mol.GetNumBonds()
            tanimotos = []
            adjusted_tanimotos = []
            for smol in d[1:]:
                sfp = fpCalculator(smol)
                if fpName in ["ap", "tt"] or fpName.startswith("ecfc") or fpName.startswith("fcfc"):
                    tanimoto = DataStructs.DiceSimilarity(ref_fp, sfp)
                else:
                    tanimoto = DataStructs.FingerprintSimilarity(ref_fp, sfp)
                tanimotos.append(tanimoto)

            label = fpName
            writer.write_result(label, tanimotos, i==0)

class Writer(object):
    def __init__(self, resultsdir):
        self.files = {}
        self.resultsdir = resultsdir
    def write_result(self, fpname, fpresult, deletefile):
        ## Write out the results
        fname = os.path.join(self.resultsdir, fpname + ".txt")
        if deletefile:
            if os.path.isfile(fname):
                os.remove(fname)
            self.files[fname] = open(fname, "w")
        f = self.files[fname]
        f.write(" ".join(map(str, fpresult)))
        f.write("\n")

def write_results(fpname, resultsdir, fpresults):
    ## Write out the results
    fname = open(os.path.join(resultsdir, fpname + ".txt"), "w")
    for fpresult in fpresults:
        fname.write(" ".join(map(str, fpresult)))
        fname.write("\n")
    fname.close()

def get_rdkitmols(dataset):
    for d in dataset:
        tmp = []
        for smi in d:
            mol = Chem.MolFromSmiles(smi)
            if "." in smi:
                frags = list(Chem.GetMolFrags(mol, asMols=True))
                frags.sort(key=lambda x:x.GetNumHeavyAtoms(), reverse=True)
                mol = frags[0]
            tmp.append(mol)
        yield tmp

if __name__ == "__main__":
    for benchmark in ["SingleAssay", "MultiAssay"]:
        # Note that the following loop is completely parallelisable
        # e.g. you could run from 0->500 on one CPU and from 500->1000 on
        #      another to finish in half the time
        if not os.path.isdir(os.path.join(benchmark, "similarities")):
            os.mkdir(os.path.join(benchmark, "similarities"))
        for M in range(1000):
            print "\nITERATION %d\n" % M
            filename = os.path.join(benchmark, "dataset", "%d.txt" % M)
            dataset = list(readBenchmark(filename))

            d = []
            for data in dataset:
                d.append([chembl.smiles_lookup[x] for x in data])
            evaluate_similarity_method(d, os.path.join(benchmark, "similarities", str(M)))
