#from Bio import AlignIO
#align = AlignIO.read("clustalw.aln", "clustal")

#print(align)
#print(len(align))

#for record in align:
   #print("%s %i" % (record.seq, len(record)))

#print(align[0].seq)

#print(align[-1].seq)

"""for row in align:
    print(row.id[1], end='')
"""

"""for row in align:
    print(row.id[ :10],row.seq[:])
"""  
#for row in align:
 #   print(row.id[:10]  + row.id[-10:])

"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

a = SeqRecord(Seq("AAAACGT"), id="Alpha", description="", name="", dbxrefs=[], features=[], annotations={"molecule_type": "DNA"})
b = SeqRecord(Seq("AAA-CGT"), id="Beta", description="", name="", dbxrefs=[], features=[], annotations={"molecule_type": "DNA"})
c = SeqRecord(Seq("AAAAGGT"), id="Gamma", description="", name="", dbxrefs=[], features=[], annotations={"molecule_type": "DNA"})

align = MultipleSeqAlignment([a, b, c],
                             annotations={"tool": "demo"},
                             column_annotations={"stats": "CCCXCCC"})
print(align)
"""
"""
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

seq_alpha = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha", description="")
seq_beta = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta", description="")
seq_gamma = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma", description="")
align = MultipleSeqAlignment([seq_alpha, seq_beta, seq_gamma])
print(align)
"""
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

seq_alpha = Seq("ACTGCTAGCTAG")
seq_beta = Seq("ACT-CTAGCTAG")
seq_gamma = Seq("ACTGCTAGATAG")

record_alpha = SeqRecord(seq_alpha, id="Alpha", description="")
record_beta = SeqRecord(seq_beta, id="Beta", description="")
record_gamma = SeqRecord(seq_gamma, id="Gamma", description="")

align = MultipleSeqAlignment([record_alpha, record_beta, record_gamma])

with open("output.fasta", "w") as output_handle:
    for record in align:
        output_handle.write(">%s\n%s\n" % (record.id, record.seq))

with open("output.fasta", "r") as output_handle:
    for line in output_handle:
        if line.startswith(">"):
            print(line.strip())
        else:
            print(line.strip().upper())
"""
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

seq_alpha = Seq("ACTGCTAGCTAG")
seq_beta = Seq("ACT-CTAGCTAG")
seq_gamma = Seq("ACTGCTAGATAG")

record_alpha = SeqRecord(seq_alpha, id="Alpha", description="")
record_beta = SeqRecord(seq_beta, id="Beta", description="")
record_gamma = SeqRecord(seq_gamma, id="Gamma", description="")

align = MultipleSeqAlignment([record_alpha, record_beta, record_gamma])

with open("output.phy", "w") as output_handle:
    AlignIO.write(align, output_handle, "phylip-relaxed")
    
with open("output.phy", "r") as output_handle:
    print(output_handle.read())
"""
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

seq_alpha = Seq("ACTGCTAGCTAG")
seq_beta = Seq("ACT-CTAGCTAG")
seq_gamma = Seq("ACTGCTAGATAG")

record_alpha = SeqRecord(seq_alpha, id="Alpha", description="")
record_beta = SeqRecord(seq_beta, id="Beta", description="")
record_gamma = SeqRecord(seq_gamma, id="Gamma", description="")

align = MultipleSeqAlignment([record_alpha, record_beta, record_gamma])
for record in align:
    print(record.id)
    print(record.seq)

"""
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

seq_alpha = Seq("ACTGCTAGCTAG")
seq_beta = Seq("ACT-CTAGCTAG")
seq_gamma = Seq("ACTGCTAGATAG")

record_alpha = SeqRecord(seq_alpha, id="Alpha")
record_beta = SeqRecord(seq_beta, id="Beta")
record_gamma = SeqRecord(seq_gamma, id="Gamma")

align = MultipleSeqAlignment([record_alpha, record_beta, record_gamma])

    
print("Alignment length:", align.get_alignment_length())
print(len(align))

"""

"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

a = SeqRecord(Seq("AAAACGT"), id="Alpha")
b = SeqRecord(Seq("AAA-CGT"), id="Beta")
c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
d = SeqRecord(Seq("AAAACGT"), id="Delta")
e = SeqRecord(Seq("AAA-GGT"), id="Epsilon")

align = MultipleSeqAlignment([a, b, c, d, e])
align.extend([d, e])
print(align)
#align = MultipleSeqAlignment([a, b, c])
#print(align)
"""
"""
from Bio import AlignIO
align = AlignIO.read("clustalw.aln", "clustal")
#print(align)
#print(len(align))

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
dummy = SeqRecord(Seq("N"*30), id="dummy")
align.append(dummy)
print(align)
print(len(align))
"""
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

a1 = SeqRecord(Seq("AAAAC"), id="Alpha")
b1 = SeqRecord(Seq("AAA-C"), id="Beta")
c1 = SeqRecord(Seq("AAAAG"), id="Gamma")
a2 = SeqRecord(Seq("GT"), id="Alpha")
b2 = SeqRecord(Seq("GT"), id="Beta")
c2 = SeqRecord(Seq("GT"), id="Gamma")

left = MultipleSeqAlignment([a1, b1, c1],
                            annotations={"tool": "demo", "name": "start"},
                            column_annotations={"stats": "CCCXC"})
right = MultipleSeqAlignment([a2, b2, c2],
                             annotations={"tool": "demo", "name": "end"},
                             column_annotations={"stats": "CC"})
print(left)
print(right)

combined = left + right
print(combined)
print(len(left))
print(len(right))
print(len(combined))
"""
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
a = SeqRecord(Seq("AAAACGT"), id="Alpha")
b = SeqRecord(Seq("AAA-CGT"), id="Beta")
c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
d = SeqRecord(Seq("AAAACGT"), id="Delta")
e = SeqRecord(Seq("AAA-GGT"), id="Epsilon")
align = MultipleSeqAlignment([a, b, c, d, e])
#first_record = align[0]
#print("%s %s" % (first_record.id, first_record.seq))

#last_record = align[-1]
#print("%s %s" % (last_record.id, last_record.seq))

#sub_alignment = align[2:5]
#print(sub_alignment)

#sub_alignment = align[::2]
#print(sub_alignment)

#rev_alignment = align[::-1]
#print(rev_alignment)

print(align[3, 4])
print(align[3][4])
print(align[3].seq[4])
print(align[:, 4])
print(align[1:3, 4])
print(align[1:5, 3:6])
"""
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils import GC

align1 = MultipleSeqAlignment([
             SeqRecord(Seq("ACGT"), id="Human"),
             SeqRecord(Seq("ACGG"), id="Mouse"),
             SeqRecord(Seq("ACGC"), id="Chicken"),
         ])

align2 = MultipleSeqAlignment([
             SeqRecord(Seq("CGGT"), id="Mouse"),
             SeqRecord(Seq("CGTT"), id="Human"),
             SeqRecord(Seq("CGCT"), id="Chicken"),
         ])
print(align1 + align2)
align1.sort()
align2.sort()
print(align1 + align2)

print(align1)

align1.sort(reverse=True)
print(align1)
"""
"""
from Bio import Align
aligner = Align.PairwiseAligner()
alignments = aligner.align("GAACT", "GAT")
alignment = alignments[0]
print(alignment)

alignment = alignments[1]
print(alignment)

aligner.mismatch_score = -10
alignments = aligner.align("AAACAAA", "AAAGAAA")
print(len(alignments))

print(alignments[0])
print(alignments[1])
"""

"""
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -0.5
aligner.extend_gap_score = -0.1
aligner.target_end_gap_score = 0.0
aligner.query_end_gap_score = 0.0
aligner.mode = 'local'
for alignment in aligner.align("TACCG", "ACG"):
    print("Score = %.1f:" % alignment.score)
    print(alignment)
"""


from Bio.Align import substitution_matrices, PairwiseAligner
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
alignments = aligner.align("KEVLA", "EVL")
alignments = list(alignments)
print("Number of alignments: %d" % len(alignments))

alignment = alignments[0]
print("Score = %.1f" % alignment.score)

print(alignment)
