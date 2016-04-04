The dataset
===========

Unzip the data
--------------
The data has been compressed with 7zip. It can be unzipped with p7zip on Linux, or 7zip on Windows.

Description
-----------
The benchmark data will be placed in the dataset subdirectories of the
SingleAssay and MultiAssay directories. There are 1000 files corresponding
to the 1000 repetitions. Each file contains several thousand lines of
CHEMBL IDs, where the first ID is the reference molecule, and the other
four are molecules are increasing distance (decreasing similarity) to the
reference.

How to reproduce the results
============================

Requirements
------------
1. Python 2.7
2. NumPy
3. SciPy
4. RDKit (2015.09.2)

Optional but needed to generate the graph depictions
----------------------------------------------------
1. dot (provided by GraphViz)

Get ChEMBL
----------
1. Download ChEMBL20 as an SDF file
2. Convert it to a SMILES file where the title field is the numeric portion of
   the CHEMBLID. The details are left to the reader. Once done, the file should
   look something like this:

    Cc1cc(cn1C)c2csc(n2)N=C(N)N	153534
    COc1cc(ccc1OC(=O)C23CC4CC(C2)CC(C4)C3)CC=C	265174
    Cc1cccc(c1)N2CCN(CC2)CCCON3C(=O)c4ccccc4C3=O	264472
    c1ccc2c(c1)n(c(=N)s2)CCN3CCC(CC3)c4ccc(cc4)F	405225

3. Name this file chembl_20.smi and place it in the benchlib directory.

Run and analyse the benchmark
-----------------------------
1. python 1-Similarities.py
2. python 2-Correlations.py
3. python 3-AnalyseResults.py
4. dot SingleAssay\graph.gv -T png > singleassay.png
5. dot MultiAssay\graph.gv -T png > multiassay.png

Notes
-----
1. Running the Python scripts on one CPU may take some time. To speed things up,
   you may wish to parallelise the main loops. This is left to the reader.
