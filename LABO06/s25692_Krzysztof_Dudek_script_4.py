import primer3

sequence = "GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG"

primer_design_settings = {
    'SEQUENCE_TEMPLATE': sequence,
    'SEQUENCE_TARGET': [100, 250],
    'PRIMER_NUM_RETURN': 2
}

result = primer3.bindings.design_primers(primer_design_settings, {})

with open("PCR_primers.txt", "w") as file:
    for i in range(primer_design_settings['PRIMER_NUM_RETURN']):
        left_primer_key = f'PRIMER_LEFT_{i}_SEQUENCE'
        right_primer_key = f'PRIMER_RIGHT_{i}_SEQUENCE'

        left_primer_seq = result[left_primer_key]
        right_primer_seq = result[right_primer_key]

        left_pos, left_len = result[f'PRIMER_LEFT_{i}']
        right_pos, right_len = result[f'PRIMER_RIGHT_{i}']

        primer_start = left_pos
        primer_end = right_pos + right_len - 1

        print(f"{left_primer_seq} -- {right_primer_seq} {primer_start} {primer_end}")
        file.write(f"{left_primer_seq} -- {right_primer_seq} {primer_start} {primer_end}\n")
