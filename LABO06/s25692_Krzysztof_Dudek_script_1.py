from Bio import Align

seq_1 = "AAAACCCCTTGGTTTT"
seq_2 = "AAACACCCTTTGTTTA"

aligner = Align.PairwiseAligner()

score = aligner.score(seq_1, seq_2)
print(f"Wartość dopasowania pomiędzy sekwencjami: {score}")

alignments = aligner.align(seq_1, seq_2)
print("Pierwsze 3 przyrównania:")
for i, alignment in enumerate(alignments):
    if i < 3:  # wyświetlmy tylko 3 pierwsze
        print(f"Alignment #{i + 1}:")
        print(alignment)
        print()
    else:
        break

for i, alignment in enumerate(alignments):
    if i < 3:
        print(f"Alignment #{i + 1} macierz koordynat:")
        print(alignment.coordinates)
        print()
    else:
        print(
            "Dzięki tym macierzom można zorientować się, między jakimi indeksami znajdują się dopasowane do siebie sekwencje\n")
        break


def count_mismatches(alignments):
    """
    Zwraca liczbę sekwencji z conajmniej jednym niedopadowaniem
    alignments: obiekt PairwiseAlignments z biblioteki Bio.Align
    """
    count_nonzero_mismatch = 0
    for alignment in alignments:
        if alignment.counts().mismatches > 0:
            count_nonzero_mismatch += 1
    return count_nonzero_mismatch


print("Liczba przyrównań z co najmniej jednym niedopasowaniem:", count_mismatches(alignments))
