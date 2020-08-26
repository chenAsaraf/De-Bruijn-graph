from Bio import SeqIO


def some_function(name, sequence):
    print(name, sequence)

contigs_file = "output_prefixx_contigs.unitigs.fa"
fasta_sequences = SeqIO.parse(open(contigs_file),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    new_sequence = some_function(name, sequence)


