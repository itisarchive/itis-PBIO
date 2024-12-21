from Bio import Align

seq_1 = "AAGAAATTCCAAGTCCAGGGATACACAAACAGGTGTACAGCAAATCATGTAGGTGGTACTTTTCCCCTAAGTTATAATTT\
AATCTATACCCTAGGAAAATGCCAAAGTCACAATTGGGTGGATTGGGTGATTTTCCAGTAGAAAGAAAATTCCATCCCAT"
seq_2 = "AAGAAATTCACAGTCCAGGGATACACAAACAGGTGTACAAAGCAAATCATGTAGGTGGTACTTTTCCCCTAAGTTATAATTT\
AATCTATACCCTAGGAAAATGCCAAAGTCGGAATTGGGTGATTGGGTGTTTCCAGTAGAAAGAAAATTCCATCCCAT"

aligner = Align.PairwiseAligner()
alignments = aligner.align(seq_1, seq_2)

max_len = 0
start_position = None

for alignment in alignments:
    seqA = alignment[0]
    seqB = alignment[1]

    current_len = 0
    current_start = 0

    for i, (a, b) in enumerate(zip(seqA, seqB)):
        if a == b and a != '-':
            current_len += 1
        else:
            if current_len > max_len:
                max_len = current_len
                start_position = i - current_len
            current_len = 0

    if current_len > max_len:
        max_len = current_len
        start_position = len(seqA) - current_len

print("Najdłuższy fragment pasujących nukleotydów:", max_len)
print("Początek fragmentu (w wyrównaniu):", start_position)
