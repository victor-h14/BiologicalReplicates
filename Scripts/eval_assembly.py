import numpy as np
import matplotlib.pyplot as plt

files = ["../Assembly/Assembly_default/default.Trinity/default.Trinity.fasta", "SBC.fasta", "SBDG.fasta"]
plt.figure(1, figsize=(1, 1))


for i in range(len(files)):
    f1 = open(files[i])
    text = f1.read().split('\n')
    
    del text[-1]
    
    lens = []
    
    for line in text:
        if line[0] == '>':
            t = line.split(' ')[1]
            
            lens.append(int(t[4:]))
    
    f2 = open(files[i] + ".gene_trans_map")
    text = f2.read().split("\n")
    
    del text[-1]
    
    dict_gene = {}
    
    for line in text:
        gene = line.split("\t")[0]
        if gene not in dict_gene:
            dict_gene[gene] = 1
        else:
            dict_gene[gene] += 1
    
    dict_isoforms = {"1": 0, "2":0, "3":0, "4":0, "5":0, "6 a 10":0, "> 10":0}
    for n in dict_gene.values():
        if n == 1:
            dict_isoforms["1"] += 1
        elif n == 2:
            dict_isoforms["2"] += 1
        elif n == 3:
            dict_isoforms["3"] += 1
        elif n == 4:
            dict_isoforms["4"] += 1
        elif n == 5:
            dict_isoforms["5"] += 1
        elif n > 10:
            dict_isoforms["> 10"] += 1
        elif n > 5:
            dict_isoforms["6 a 10"] += 1
    
    
    
    slices = len(dict_isoforms)
    labels = dict_isoforms.keys()
    sizes = dict_isoforms.values()
    
    plt.subplot(len(files),1, i+1)
    plt.title("Assembly " + str(i + 1))
    plt.bar(labels, sizes)
    
    print("number of genes:", len(dict_gene))
    print("number of transcripts:", len(lens))
    print("mean length of transcripts:", np.round(sum(lens)/len(lens), 2))
    print("median length of transcripts:", np.median(lens))
    
    print("genes w/ 1 transcript:", np.round(sum([i == 1 for i in dict_gene.values()])/len(dict_gene.values()), 3))
    print("genes w/ 2 transcript:", np.round(sum([i == 2 for i in dict_gene.values()])/len(dict_gene.values()), 3))
    print("genes w/ 3 transcript:", np.round(sum([i == 3 for i in dict_gene.values()])/len(dict_gene.values()), 3))
    print("genes w/ > 4 transcript:", np.round(sum([i > 3 for i in dict_gene.values()])/len(dict_gene.values()), 3))
    print("largest number of transcripts for a gene:", max(dict_gene.values()))
    
    max_value = max(dict_gene.values())
    for key in dict_gene.keys():
        if  dict_gene[key] == max_value:
            print(key)

plt.show()

