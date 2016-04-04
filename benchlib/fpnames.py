import os

here = os.path.dirname(os.path.realpath(__file__))

fpnames = {}
for line in open(os.path.join(here, "fingerprint_lib.py")):
    if line.startswith("fpdict["):
        name = line.split()[0][8:-2]
        fpnames[name] = name.upper()
fpnames["avalon"] = "Avalon"
fpnames["laval"] = "LAvalon"
fpnames["hashtt"] = "HashTT"
fpnames["hashap"] = "HashAP"

if __name__ == "__main__":
    print fpnames, len(fpnames)
