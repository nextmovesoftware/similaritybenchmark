import os
import sys
import glob
import scipy.stats

SIZE = 4

def getCorrelations(fname, outfile):
    ref_corr = range(SIZE, 0, -1)
    with open(outfile, "w") as f:
        for line in open(fname):
            similarities = [float(x) for x in line.split()]
            if len(set(similarities)) == 1:
                corr = (0, None)
            else:
                corr = scipy.stats.spearmanr(ref_corr, similarities)
            f.write("%f\n" % corr[0])

if __name__ == "__main__":
    for benchmark in ["SingleAssay", "MultiAssay"]:
        # Note that the following loop is completely parallelisable
        # e.g. you could run from 0->500 on one CPU and from 500->1000 on
        #      another to finish in half the time
        for i in range(1000):
            print i,
            if not os.path.isdir(os.path.join(benchmark, "correlations")):
                os.mkdir(os.path.join(benchmark, "correlations"))
            inputdir = os.path.join(benchmark, "similarities", "%d" % i)
            outputdir = os.path.join(benchmark, "correlations", "%d" % i)
            if not os.path.isdir(outputdir):
                os.mkdir(outputdir)
            for fname in glob.glob(os.path.join(inputdir, "*.txt")):
                name = os.path.basename(fname)
                getCorrelations(fname, os.path.join(outputdir, name))
