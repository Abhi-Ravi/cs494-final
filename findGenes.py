from Bio import Entrez, SeqIO
from Bio.SeqUtils import CodonUsage as CU
from CAI import CAI

import heapq

fp = open("e-coli.gtf")
tmp = fp.readlines()

max_size = float('-inf')
max_gene = ''
max_heap = []
heapq.heapify(max_heap)

min_size = float('inf')
min_gene = ''
min_heap = []
heapq.heapify(min_heap)

for line in tmp:
    line = line.replace(';', ' ')
    line_array = (line.strip()).split()
    if len(line_array) >= 3 and line_array[2] == 'gene':
        gene_start = int(line_array[3])
        gene_end = int(line_array[4])

        gene_size = abs(gene_end-gene_start)

        gene_name = line_array[9]

        gene_data = (gene_size, gene_name, gene_start, gene_end)

        heapq.heappush(max_heap, gene_data) 




largest_genes = heapq.nlargest(5, max_heap)
smallest_genes = heapq.nsmallest(5,max_heap)


print(largest_genes)
print(smallest_genes)

# handle = Entrez.efetch(db="nuccore", id="1109557564", seq_start="4209089", seq_stop="4209144",rettype="fasta", retmode="text")
# record = SeqIO.read(handle, "fasta")
# # handle.close()
# # print(record.seq)

# # genome_record = SeqIO.read()
# fp = open("e-coli.fna")
# fp.readline()

# tmp = fp.read().splitlines()
# seq = ''.join(tmp)

# # lol = CAI(seq, reference=record.seq)
# # print(lol)
records = list(SeqIO.parse("e-coli.fna", "fasta"))
# # print(records[0].seq)

# while (len(record.seq) % 3) != 0:
#     record.seq = (record.seq)[:-1]
#     print(record.seq)

# lol = CAI(records[0].seq, reference=[record.seq])

# lol = CAI(records[0].seq, reference=["TTT"])

# # for record in SeqIO.parse("salmonella.fna", "fasta"):
# #     print(record.id)
# # lol = SeqIO.read("salmonella.fna", "fasta")
# print(records[0].seq)
# print(lol)

seq = str(records[0].seq)
myIndex = CU.CodonAdaptationIndex()
myIndex.generate_index("salmonella.fna")
myIndex.print_index()