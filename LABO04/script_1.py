from collections import Counter

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


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


def transcribe(seq):
    """
    Funkcja zwracająca sekwencję RNA transkrybowaną z sekwencji DNA.
    Wykorzystuje wbudowaną metodę transcribe() z klasy Seq w Biopython.
    Metoda automatycznie zamienia wszystkie T w DNA na U w RNA.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: sekwencja RNA
    """
    return seq.transcribe()


def translate(seq):
    """
    Funkcja zwracająca sekwencję białka translatowaną z sekwencji DNA.
    Używa metody translate() z klasy Seq w Biopython.
    Dodatkowo argument to_stop=True zapewnia, że translacja kończy się na pierwszym kodonie STOP,
    co odpowiada naturalnemu procesowi biologicznemu.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: sekwencja białka
    """
    return seq.translate(to_stop=True)


def reverse_complement(seq):
    """
    Funkcja zwracająca odwrotne dopełnienie sekwencji DNA.
    Korzysta z metody reverse_complement() z klasy Seq w Biopython,
    która automatycznie odwraca kolejność nukleotydów i zamienia je na ich komplementarne bazy.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: odwrotne dopełnienie sekwencji DNA
    """
    return seq.reverse_complement()


def analyze_sequence(seq):
    """
    Funkcja analizująca sekwencję DNA.
    Zbiera wszystkie wyniki analizy sekwencji DNA w jednej strukturze danych (słowniku).
    Pozwala to uniknąć wielokrotnego przeliczania tych samych wartości.
    :param seq: sekwencja DNA (obiektu Seq z Biopython)
    :return: słownik z wynikami analizy (np. liczba nukleotydów, zawartość GC, RNA, białko, odwrotne dopełnienie)
    """
    return {
        "original": seq,
        "nucleotides": count_nucleotides(seq),
        "gc_content": gc_content(seq),
        "transcribed_rna": transcribe(seq),
        "translated_protein": translate(seq),
        "reverse_complement": reverse_complement(seq),
    }


def save_to_file(filename, results):
    """
    Funkcja zapisująca wyniki analizy sekwencji DNA do pliku.
    Przyjmuje gotowe wyniki analizy jako argument, co pozwala uniknąć ich ponownego obliczania.
    Zapisuje dane w czytelnej formie, umożliwiając dalszą analizę wyników offline.
    :param filename: nazwa pliku, do którego zostaną zapisane wyniki
    :param results: słownik z wynikami analizy sekwencji DNA
    """
    with open(filename, "w", encoding="UTF-8") as file:
        a, t, c, g = results['nucleotides']
        file.write(f"Oryginalna sekwencja DNA: {results['original']}\n"
                   f"Liczba nukleotydów:\n"
                   f"A: {a}\nT: {t}\nC: {c}\nG: {g}\n"
                   f"Zawartość GC: {results['gc_content']:.2f}%\n"
                   f"Transkrybowany RNA: {results['transcribed_rna']}\n"
                   f"Translowane białko: {results['translated_protein']}\n"
                   f"Odwrotne dopełnienie: {results['reverse_complement']}\n")


def main():
    """
    Funkcja główna programu.
    Sekwencja DNA jest tworzona jako obiekt Seq z Biopython.
    Następnie wywoływana jest funkcja analyze_sequence, która zwraca wyniki analizy w postaci słownika.
    Wyniki są wyświetlane w konsoli i zapisywane do pliku.
    """
    dna = "AAGAAATTCCAAGTCCAGGGATACACAAACAGGTGTACAGCAAATCATGTAGGTGGTACTTTTCCCCTAAGTTATAATATT"
    seq = Seq(dna)
    results = analyze_sequence(seq)
    a, t, c, g = results['nucleotides']
    print(f"Oryginalna sekwencja DNA: {results['original']}\n"
          f"Liczba nukleotydów:\n"
          f"A: {a}\nT: {t}\nC: {c}\nG: {g}\n"
          f"Zawartość GC: {results['gc_content']:.2f}%\n"
          f"Transkrybowany RNA: {results['transcribed_rna']}\n"
          f"Translowane białko: {results['translated_protein']}\n"
          f"Odwrotne dopełnienie: {results['reverse_complement']}")
    save_to_file("sequence_analysis.txt", results)


if __name__ == "__main__":
    main()
