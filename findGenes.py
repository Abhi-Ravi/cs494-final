from Bio import Entrez, SeqIO
from Bio.SeqUtils import CodonUsage as CU
from CAI import CAI

import matplotlib.pyplot as plt
import time
import heapq

fp = open("salmonella.gtf")
tmp = fp.readlines()

max_size = float('-inf')
max_gene = ''
max_heap = []
heapq.heapify(max_heap)

min_size = float('inf')
min_gene = ''
min_heap = []
heapq.heapify(min_heap)

gene_arr = []

for line in tmp:
    line = line.replace(';', ' ')
    line_array = (line.strip()).split()
    if len(line_array) >= 3 and line_array[2] == 'gene':
        gene_start = line_array[3]
        gene_end = line_array[4]

        gene_size = abs(int(gene_end)-int(gene_start))

        gene_name = line_array[9]

        gene_id = line_array[11].split(':')
        gene_id = gene_id[1]

        gene_data = [gene_size, gene_name, gene_start, gene_end]
        print(line_array)
        heapq.heappush(max_heap, gene_data) 

        gene_arr.append(gene_data)

print(len(gene_arr))
largest_genes = heapq.nlargest(50, max_heap)
smallest_genes = heapq.nsmallest(50,max_heap)

Entrez.email = 'aravi@vols.utk.edu'

gene_seqs = []

ecoli_record = list(SeqIO.parse("e-coli.fna", "fasta"))
ecoli_seq = ecoli_record[0].seq

while (len(ecoli_seq) % 3) != 0:
    ecoli_seq = (ecoli_seq)[:-1]

x_large = []
y_large = []

x_small = []
y_small = []

fig, ax = plt.subplots(nrows=1, ncols=2)

for gene in largest_genes:
    handle = Entrez.efetch(db="nuccore", id="1109557564", seq_start=gene[2], seq_stop=gene[3],rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    gene_seq = str(record.seq)

    while (len(gene_seq) % 3) != 0:
        gene_seq = (gene_seq)[:-1]

    gene_cai = CAI(ecoli_seq, reference=[gene_seq])
    gene.append(gene_cai)
    x_large.append(gene[0])
    y_large.append(gene[4])

for gene in smallest_genes:
    handle = Entrez.efetch(db="nuccore", id="1109557564", seq_start=gene[2], seq_stop=gene[3],rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    gene_seq = str(record.seq)

    while (len(gene_seq) % 3) != 0:
        gene_seq = (gene_seq)[:-1]

    gene_cai = CAI(ecoli_seq, reference=[gene_seq])
    gene.append(gene_cai)
    
    x_small.append(gene[0])
    y_small.append(gene[4])


print("LARGEST:", largest_genes)
print("SMALLEST: ", smallest_genes)

ax[0].scatter(x_small,y_small)
ax[0].set_title("Salmonella Enterica's 50 smallest gene lengths relative to CAI in reference to E. Coli")
ax[0].set_xlabel("Number of bases")
ax[0].set_ylabel("CAI INDEX")
ax[0].set_ylim([0,1])

ax[1].scatter(x_large, y_large)
ax[1].set_title("Salmonella Enterica's 50 largest gene lengths relative to CAI in reference to E. Coli")
ax[1].set_xlabel("Number of bases")
ax[1].set_ylabel("CAI INDEX")
ax[1].set_ylim([0,1])

plt.show()
