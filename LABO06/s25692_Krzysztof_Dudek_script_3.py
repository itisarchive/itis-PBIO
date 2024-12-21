import primer3

def analyze_primers(filename):
    results = []
    with open(filename, 'r') as f:
        for line in f:
            primer = line.strip()
            if not primer:
                continue

            tm = primer3.calc_tm(primer)
            homodimer_result = primer3.bindings.calc_homodimer(primer)
            if homodimer_result.dg < 0:
                homodimer = "Tak"
            else:
                homodimer = "Nie"
            hairpin_result = primer3.bindings.calc_hairpin(primer)
            if hairpin_result.dg < 0:
                hairpin = "Tak"
            else:
                hairpin = "Nie"
            results.append((primer, tm, homodimer, hairpin))
    return results

def main():
    filename = "primers.txt"
    results = analyze_primers(filename)

    print("Primer,Tm,Homodimer,Hairpin")
    with open("primer_analysis.csv", "w") as file:
        for (primer, tm, homodimer, hairpin) in results:
            print(f"{primer},{tm:.2f},{homodimer},{hairpin}")
            file.write(f"{primer},{tm:.2f},{homodimer},{hairpin}\n")


if __name__ == "__main__":
    main()
