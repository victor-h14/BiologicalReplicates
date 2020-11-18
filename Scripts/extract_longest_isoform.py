f1 = open("Master/Assembly/Assembly_default/default.Trinity.fasta")
text = f1.read().split('\n')

del text[-1]

dict_gene = {}
for i in range(len(text)):
    if text[i][0] == '>':
        line = text[i]
        head = line.split(' ')
        gene = head[0][:line.index('i')-1]
        if gene not in dict_gene:
            dict_gene[gene] = text[i+1]
        elif len(dict_gene[gene]) < len(text[i+1]):
            dict_gene[gene] = text[i+1]
f2 = open("Master/Assembly/Assembly_default/longest-isoform.fasta", "w")

for gene in dict_gene.keys():
    f2.write(gene + "\n" + dict_gene[gene] + "\n")
    
f2.flush()
f2.close()
