import os
import glob
from collections import Counter, defaultdict
import cPickle as pickle
import doctest
from benchlib.fpnames import fpnames
import numpy as np
import random
random.seed(1)

import scipy.stats

def assignRanks(data):
    """
    >>> ans = assignRanks([(1.0, "A"), (1.0, "B"), (0.8, "C"), (0.8, "D"), (0.7, "E"), (0.5, "F")])
    >>> ans["A"]
    1
    """
    rank = 1
    ranks = [(data[0][1], rank)]
    for i in range(1, len(data)):
        rank += 1
        d = data[i]
        if d[0]==data[i-1][0]:
            newrank = ranks[i-1][1]
        else:
            newrank = rank
        ranks.append( (d[1], newrank) )
    return dict(ranks)

def analyseBenchmarkResults(folder, minsize):

    results = {}
    for fname in glob.glob(os.path.join(folder, "*.txt")):
        name = os.path.basename(fname).split(".")[0]
        results[name] = [float(x) for x in open(fname)][:minsize]

    size = len(next(results.itervalues()))
    names = results.keys()

    # For each sequence rank the methods by which was the best
    data = [[] for x in range(size)]
    for k, v in results.iteritems():
        for i in range(len(v)):
            data[i].append( (v[i], k) )

    ranks = defaultdict(list)
    better = defaultdict(int)
    worse = defaultdict(int)
    for d in data:
        d.sort(reverse=True)
        r = assignRanks(d)
        for k, v in r.iteritems():
            ranks[k].append(v)
            for a, b in r.iteritems():
                if a == k: continue
                if v < b:
                    better[ (k, a) ] += 1
                elif v > b:
                    worse[ (k, a) ] += 1

    results = {}
    for name in names:
        for nameB in names:
            if nameB < name:
                tmp = (name, nameB)
                diff = better[tmp] - worse[tmp]
                results[tmp] = diff

    # Return a dictionary of (fpA, fpB) tuples where the value is the difference between how many times
    # fpA was better than fpB compared to how many times it was worse
    return results

def getnames(allpairs):
    names = set()
    for k in allpairs.iterkeys():
        names.add(k[0])
        names.add(k[1])

    return names

def removeInsignificant(allpairs):
    names = getnames(allpairs)

    tmp = []
    distances = {}
    for name, data in allpairs.iteritems():
        if name[0] not in names or name[1] not in names: continue
        signif = scipy.stats.ttest_1samp(data, 0)[1]
        mean = np.mean(data)
        if mean > 0: # reorder the tuples so that they indicate A>B
            tmp.append( (signif, name) )
            distances[name] = mean
        else:
            tmp.append( (signif, (name[1], name[0])) )
            distances[(name[1], name[0])] = -mean
    tmp.sort()

    # Holm-Bonferroni correction (https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method)
    results, notsignif = [], []
    m = len(tmp)
    print "Number of hypotheses", m
    alpha = 0.01
    for k in range(1, m+1):
        pvalue, name = tmp[k-1]
        if pvalue > (alpha/float(m+1-k)):
            for l in range(k, m+1):
                pvalue, name = tmp[l-1]
                print "Not significant: %s > %s (%f)" % (name[0], name[1], pvalue)
                notsignif.append(name)
            break
        else:
            # It's significant
            results.append( (name, mean) )

    return list(names), dict(results), distances, notsignif

def calcBetterOrWorse(names, origpairs):
    # Use each fp to partition each of the fps
    better = defaultdict(set) # list of fps that are better than a particular fp
    worse = defaultdict(set)  # list of fps that are worse than a particular fp
    for ref in names:
        for test in names:
            if (test, ref) in origpairs:
                better[ref].add(test)
            if (ref, test) in origpairs:
                worse[ref].add(test)
    return better, worse

def removeIncomparable(names, origpairs, better, worse):

    # Create list of pairs of X > Y
    isComparable = True
    pairs = []
    N = len(names)
    for i in range(N-1):
        name = names[i]
        for j in range(i+1, N):
            nameB = names[j]
            if (name, nameB) not in origpairs and (nameB, name) not in origpairs: continue
            # ensure that name>nameB or nameB>name
            comparisons = []
            remember = []
            for testname in names:
                scoreA = scoreB = 0
                if name in better[testname]:
                    scoreA = 1
                elif name in worse[testname]:
                    scoreA = -1
                if nameB in better[testname]:
                    scoreB = 1
                elif nameB in worse[testname]:
                    scoreB = -1
                if scoreA != scoreB:
                    comparisons.append(scoreA > scoreB)
                    remember.append(testname)
            if comparisons and (all(comparisons) or all(not x for x in comparisons)):
                # It's comparable!
                pairs.append( (name, nameB) if comparisons[0] else (nameB, name))
            else:
                isComparable = False
                print "%s and %s not comparable" % (name, nameB)
                for x, testname in zip(comparisons, remember):
                    if name in better[testname]:
                        print " %s > %s" % (name, testname)
                    elif name in worse[testname]:
                        print " %s < %s" % (name, testname)
                    elif name != testname:
                        print " %s <> %s" % (name, testname)
                    if nameB in better[testname]:
                        print " %s > %s" % (nameB, testname)
                    elif nameB in worse[testname]:
                        print " %s < %s" % (nameB, testname)
                    elif nameB != testname:
                        print " %s <> %s" % (nameB, testname)


    if not isComparable:
        print "It's not comparable!"

    return pairs

def create_graph(names, better, distances):
    """Create a directed graph with edges from A->B if there is
    evidence that A>B"""
    graph = defaultdict(list)
    for x, y in better:
        graph[x].append(y)

    # prune graph
    ngraph = dict(graph)
    for name in names:
        stack = []
        stack.append((name, []))
        while stack:
            curr, path = stack.pop()
            # prune
            for p in path[:-1]:
                if curr in ngraph[p]:
                    ngraph[p].remove(curr)
            for nbr in graph[curr]:
                stack.append( (nbr, path+[curr]) )

    dot = ["digraph {"]
    for name in names:
        label = name.split("_")[-1]
        label = fpnames.get(label, label)
        dot.append("%s [label=%s fontsize=%d]" % (name, label, 50))
    for name, nbrs in ngraph.iteritems():
        for nbr in nbrs:
            dot.append("%s -> %s [label=%.0f fontsize=%d]" % (name, nbr, distances[(name, nbr)], 40))
    dot.append("}")
    return "\n".join(dot)

def countlines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def saveDifferenceMatrix(benchmark, names, distances, notsignif):
    snames = names[:]
    snames.sort()
    fname = os.path.join(benchmark, "net_difference.tsv")
    with open(fname, "w") as f:
        f.write("\t"+"\t".join([fpnames[x] for x in snames]))
        f.write("\n")
        for a in snames:
            tmp = [fpnames[a]]
            for b in snames:
                star = ""
                if a==b:
                    distance = "-"
                else:
                    lookup = (a, b)
                    if lookup not in distances:
                        lookup = (lookup[1], lookup[0])
                        distance = "%.0f" % -distances[lookup]
                    else:
                        distance = "%.0f" % distances[lookup]
                    if lookup in notsignif:
                        star = " *"
                tmp.append(distance + star)
            f.write("\t".join(tmp))
            f.write("\n")

if __name__ == "__main__":
    for benchmark in ["SingleAssay", "MultiAssay"]:
        correlationdir = os.path.join(benchmark, "correlations")
        doctest.testmod()

        minsize = 99999999
        N = 1000
        for i in range(N):
            size = countlines(os.path.join(correlationdir, "%d" % i, "ecfp0.txt"))
            if size < minsize:
                minsize = size
        print "Using minsize of %d" % minsize

        allpairs = defaultdict(list)
        for i in range(N):
            print i,
            pairs = analyseBenchmarkResults(os.path.join(correlationdir, "%d" % i), minsize)
            for k, v in pairs.iteritems():
                allpairs[k].append(v)
        print

        names, pairs, alldistances, notsignif = removeInsignificant(allpairs)
        saveDifferenceMatrix(benchmark, names, alldistances, notsignif)

        better, worse = calcBetterOrWorse(names, pairs)
        pairs = removeIncomparable(names, pairs, better, worse)
        with open(os.path.join(benchmark, "graph.gv"), "w") as f:
            f.write(create_graph(names, pairs, alldistances))

        print """
    Create a PNG of the graph using dot (from GraphViz) as follows:
              dot graph.gv -T png > tmp.png
    """
