import numpy as np
import matplotlib.pyplot as plt
import os

path = "Master/Quantification/Default_longest/quants_out/"
feature = "longest isoform"
tpm = {}
samples = 0
for dirs in os.listdir(path):
    if feature == "genes":
        string = "quant.genes.sf"
    else:
        string = "quant.sf"
    file = open(os.path.join(path, dirs, string))
    lines = file.read().split('\n')
    
    del lines[0]
    del lines[-1]
    
    samples += 1
    for line in lines:
        cols = line.split("\t")
        if cols[0] not in tpm:
            tpm[cols[0]] = [float(cols[3])]
        else:
            tpm[cols[0]] += [float(cols[3])]


means = [sum(tpm[i])/samples for i in tpm]
not_zero = []
for i in means:
    if i > 0:
        not_zero.append(i)

criteria = []
for ft in tpm:
    x = 0
    for i in range(samples):
        if tpm[ft][i] > 1:
            x += 1
    if x >= 3:
        criteria.append(sum(tpm[ft])/samples)

plt.title("Histogram of " + feature + " quantification")
plt.hist(x = np.log10(not_zero), bins = 200, label = "All non zero TPMs\n(" + str(len(not_zero)) + ")")
plt.hist(x = np.log10(criteria), bins = 200, label = "TPM > 1 for at least 3 samples\n(" + str(len(criteria)) + ")")
plt.legend(loc = "best")
plt.ylabel("Number of " + feature)
plt.xlabel("$log_{10}(TPM)$")
plt.show()
