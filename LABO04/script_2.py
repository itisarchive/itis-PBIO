import csv
from collections import Counter

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import nt_search


def count_nucleotides(seq):
    """
    Funkcja zwracająca liczbę nukleotydów w sekwencji DNA.
    Zamiast manualnego liczenia każdego nukleotydu, użyto modułu collections.Counter,
    aby uzyskać bardziej zwięzły i wydajny sposób na policzenie wystąpień nukleotydów.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: liczba nukleotydów A, T, C, G
    """
    counts = Counter(seq)
    return counts.get("A"), counts.get("T"), counts.get("C"), counts.get("G")


def gc_content(seq):
    """
    Funkcja zwracająca zawartość GC w sekwencji DNA.
    Korzysta bezpośrednio z funkcji gc_fraction z modułu Bio.SeqUtils,
    która oblicza stosunek G i C do długości sekwencji.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: zawartość procentowa GC w sekwencji DNA
    """
    return gc_fraction(seq)


def find_motifs(seq):
    """
    Funkcja zwracająca pozycje wybranych motywów w sekwencji DNA.
    Wykorzystuje funkcję nt_search z modułu Bio.SeqUtils, która zwraca pozycje
    wszystkich wystąpień danego motywu w sekwencji.
    Pierwszy element zwracanej listy to wartość motywu, dlatego pomijamy go.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: słownik z pozycjami motywów w sekwencji
    """
    motifs = ["ATG", "TATA", "GAATTC"]
    positions = {}
    for motif in motifs:
        positions[motif] = [pos for pos in nt_search(str(seq.seq), motif)[1:]]
    return positions


def reverse_complement(seq):
    """
    Funkcja zwracająca odwrotne dopełnienie sekwencji DNA.
    Korzysta z metody reverse_complement() z klasy Seq w Biopython,
    która automatycznie odwraca kolejność nukleotydów i zamienia je na ich komplementarne bazy.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: odwrotne dopełnienie sekwencji DNA
    """
    return seq.seq.reverse_complement()


def translate(seq):
    """
    Funkcja zwracająca długości translacji sekwencji DNA w trzech ramkach odczytu.
    Dla każdej ramki odczytu obliczana jest długość translacji do pierwszego kodonu STOP.
    W obu kierunkach (5' -> 3' i 3' -> 5') obliczane są translacje w trzech ramkach odczytu (od 1 do 3).
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: słownik z długościami translacji w trzech ramkach odczytu
    """
    translations = {}
    for frame in range(3):
        protein = seq[frame:].translate(to_stop=True)
        translations[f"Frame{frame + 1}"] = len(protein)
    for frame in range(3):
        protein = seq.reverse_complement()[frame:].translate(to_stop=True)
        translations[f"RevFrame{frame + 1}"] = len(protein)
    return translations


def save_as_csv(filename, data):
    """
    Funkcja zapisująca wyniki analizy sekwencji DNA do pliku CSV.
    Wykorzystuje moduł csv do zapisania danych w formacie CSV.
    :param filename: nazwa pliku CSV
    :param data: słownik z wynikami analizy sekwencji DNA
    """
    with open(filename, "w", newline="", encoding="UTF-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["SeqID", "Nucleotide_Counts", "GC_Content", "Motif_Positions",
                         "Reverse_Complement", "Translation_Lengths"])
        for seq_id, values in data.items():
            writer.writerow([seq_id, values["Nucleotide_Counts"], values["GC_Content"],
                             values["Motif_Positions"], values["Reverse_Complement"],
                             values["Translation_Lengths"]])


def main():
    """
    Główna funkcja programu, która wczytuje sekwencje z pliku FASTA, analizuje je i zapisuje wyniki do pliku CSV.
    """
    sequences = list(SeqIO.parse("ls_orchid.fasta", "fasta"))
    analized_data = {}
    for seq in sequences:
        a, t, c, g = count_nucleotides(seq)
        gc = gc_content(seq)
        motifs = find_motifs(seq)
        reverse = reverse_complement(seq)
        translations = translate(seq)
        analized_data[seq.id] = {"Nucleotide_Counts": {'A': a, 'T': t, 'C': c, 'G': g}, "GC_Content": f"{gc:.2f}%",
                                 "Motif_Positions": motifs, "Reverse_Complement": reverse,
                                 "Translation_Lengths": translations}
    save_as_csv("sequence_analysis.csv", analized_data)


if __name__ == "__main__":
    main()
