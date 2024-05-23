import csv

amino_acids_codons = {
    'Ala': ['CGA', 'CGG', 'CGT', 'CGC'], 
    'Arg': ['GCA', 'GCG', 'GCT', 'GCC', 'TCT', 'TCC'], 
    'Asn': ['TTA', 'TTG'], 
    'Asp': ['CTA', 'CTG'], 
    'Cys': ['ACA', 'ACG'], 
    'Gln': ['GTT', 'GTC'], 
    'Glu': ['CTT', 'CTC'], 
    'Gly': ['CCA', 'CCG', 'CCT', 'CCC'], 
    'His': ['GTA', 'GTG'], 
    'Ile': ['TAA', 'TAG', 'TAT'], 
    'Leu': ['AAT', 'AAC', 'GAA', 'GAG', 'GAT', 'GAC'], 
    'Lys': ['TTT', 'TTC'], 
    'Met': ['TAC'], 
    'Phe': ['AAA', 'AAG'], 
    'Pro': ['GGA', 'GGG', 'GGT', 'GGC'], 
    'Ser': ['AGA', 'AGG', 'AGT', 'AGC', 'TCA', 'TCG'], 
    'Thr': ['TGA', 'TGG', 'TGT', 'TGC'], 
    'Trp': ['ACC'], 
    'Tyr': ['ATA', 'ATG'], 
    'Val': ['CAA', 'CAG', 'CAT', 'CAC'], 
    'Stop': ['ATT', 'ATC', 'ACT']
    }


dna_to_rna_complements = {
    'A': 'U',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

# Time Complexity is O(nlogn)
def partition_sort(data, start, end):
    pivot = data[end]
    i = start - 1
    for j in range(start, end):
        if data[j] <= pivot:
            i += 1
            data[i], data[j] = data[j], data[i]
    data[i + 1], data[end] = data[end], data[i + 1]
    return i + 1

def quickSort(data, start, end):
    if start < end:
        pivot_index = partition_sort(data, start, end)
        quickSort(data, start, pivot_index - 1)
        quickSort(data, pivot_index + 1, end)



'''def suffix_search(gene, pattern):
    
    suffix_array = sorted((i, gene[i:]) for i in range(len(gene)))
    suffix_array = {k: v for k,v in suffix_array}
   
    indices = []
    for i in suffix_array:
        
        suffix = gene[i:]
        
        if suffix[:len(pattern)] == pattern:
            indices.append((i, i + len(pattern)))
    
    quickSort(indices, 0, len(indices)-1)
    return indices'''

# The suffix_search function was used during testing to ensure fm-index was outputing correctly

        
def bwt_and_suffix_array(gene): # Time Complexity is O(n)
    
    gene = gene + '$'
    rotations = [gene[i:] + gene[:i] for i in range(len(gene))]
    sorted_rotations = []
    counter = 0
    for rotation in rotations:
        sorted_rotations.append((rotation, counter))
        counter += 1
    sorted_rotations = sorted(sorted_rotations)
    bwt_string = ''.join(rotation[len(rotation)-1] for rotation, _ in sorted_rotations)
    suffix_array_map = [index for _, index in sorted_rotations]
    return bwt_string, suffix_array_map

def rank_bwt(bw):   # Time Complexity is O(n**2)
    totals = {}
    occ_ranks = []
    for char in bw:
        if char not in totals:
            totals[char] = 0
        totals[char] += 1
        # Create a new snapshot dictionary for this rank
        snapshot = {char: count for char, count in totals.items()}
        occ_ranks.append(snapshot)
    return occ_ranks


def create_c_table(bw): # Time Complexity is O(nlogn)
    total = 0
    c_table = {}
    unique_nucleotides = []
    counts = {'A': 0, 'C': 0, 'G':0, 'T': 0, 'U': 0, '$': 0}
    for base in bw:
        if base not in unique_nucleotides:
            unique_nucleotides.append(base)
        counts[base] += 1

    base_occ = []
    for n_base in unique_nucleotides:
        base_occ.append((n_base, counts[n_base]))
    del counts
    quickSort(base_occ, 0, len(base_occ)-1)

    for char, count in base_occ:
        c_table[char] = total
        total += count
    return c_table

def search_with_fm_index(pattern, bwt_string, c_table, occ_ranks, suffix_array):    # Time Complexity is O(nlogn)

    l, r = 0, len(bwt_string) - 1
    for i in range(len(pattern) - 1, -1, -1):
        char = pattern[i]
        if char not in c_table:  # Character not in BWT string, pattern cannot exist
            return -1
        if l > 0 and char in occ_ranks[l-1]:
            l = c_table[char] + occ_ranks[l-1][char]
        else:
            l = c_table[char]
        if char in occ_ranks[r]:
            r = c_table[char] + occ_ranks[r][char] - 1
        else:
            r = c_table[char] - 1
        
        if r < l:  # Pattern not found
            return -1

    indices = []
    # Reconstruct indices using suffix array mapping
    for i in range(l, r + 1):
        indices.append((suffix_array[i],suffix_array[i] + len(pattern)))
    
    quickSort(indices, 0, len(indices) - 1)
    return indices
    
def csv_reader(protein, sequence = True, name = False): # Time Complexity is O(n+m)
    myfile = open("proteinsCSV.csv", "r")
    reader = csv.reader(myfile)
    matrix = list(reader)

    index = matrix[1].index(protein)
    if index < 0:
        return 'Protein Not Found'
    name1 = matrix[index+2][0]
    sequence1 = matrix[index+2][1]
    del matrix
    if name and not sequence:
        return name1
    elif sequence and not name:
        return sequence1
    elif sequence and name:
        return name1, sequence1

def mutations(gene, mut_type, coordinates, subSequence = None): # Time Complexity is O(n**2)

    start, end = coordinates
    if mut_type == 'Insertion':

        gene = gene[:start] + subSequence + gene[start:]

    elif mut_type == 'Substitution':

        gene = gene[:start] + subSequence + gene[end:]

    elif mut_type == 'Deletion':

        gene = gene[:start] + gene[end:]
    
    bwt_string, suffix_array = bwt_and_suffix_array(gene)
    c_table = create_c_table(bwt_string)
    occ_ranks = rank_bwt(bwt_string)
    
    return gene, bwt_string, suffix_array, c_table, occ_ranks

def gene_translation(gene):  #TC is O(n**2)
    polypeptide = []
    for i in range(0, len(gene), 3):
        codon =gene[i:i+3]
        for amino_acid, codons in amino_acids_codons.items():
            if codon in codons:
                if amino_acid == 'Stop':
                    if i < len(gene) - 3:
                        position = i // 3 + 1
                        if position == 1:
                            ordinal_suffix = "st"
                        elif position == 2:
                            ordinal_suffix = "nd"
                        elif position == 3:
                            ordinal_suffix = "rd"
                        else:
                            ordinal_suffix = "th"

                        error = f'Non-sense mutation encountered at {position}{ordinal_suffix} triplet codon'

                        return '-'.join(polypeptide), polypeptide, error
                    break
                polypeptide.append(amino_acid)
    return '-'.join(polypeptide), polypeptide

def dna_to_mrna(gene):  # Time Complexity is O(n)
    mrna = ''
    for base in gene:
        mrna += dna_to_rna_complements[base]
    return mrna

def coding_to_template(coding_sequence): #TC O(n)
    template_sequence = ""
    for base in coding_sequence:
        if base == "A":
            template_sequence += "T"
        elif base == "T":
            template_sequence += "A"
        elif base == "C":
            template_sequence += "G"
        elif base == "G":
            template_sequence += "C"
    return template_sequence[::-1]

def validation_check(sequence): #Validating User Entry (Template Sequence); TC is O(n)
    if sequence == "":
        return "Please Enter A Sequence To Begin..."
    lst_codons = "ATCG"
    sequence = sequence.upper()
    if len(sequence) > 243:
        return "Sequence is too long for program to process right now..."
    for i in range(0, len(sequence)):
        if sequence[i] not in lst_codons:
            return "Invalid Codon Sequence..."
    return "Processing..."

# # # Integration and example usage
# gene = csv_reader('Immunoglobulin G')
# pattern = 'GGA'
# bwt_string, suffix_array = bwt_and_suffix_array(gene)
# occ_ranks = rank_bwt(bwt_string)
# c_table = create_c_table(bwt_string)

# matches1 = suffix_search(gene, pattern)
# # Find pattern using FM-index
# matches = search_with_fm_index(pattern, bwt_string, c_table, occ_ranks, suffix_array)

# print(f"Pattern '{pattern}' found at positions: {matches}")
# print()
# print(f"Pattern '{pattern}' found at positions: {matches1}")

# print(gene_translation(gene))


# Overall Time Complexity O(n**2)
# Overall Space Complexity O(n)

'''
Citations or learning material:

Anderson, T., & Wheeler, T. J. (2021). An optimized FM-index library for nucleotide and amino acid search. Algorithms for Molecular Biology, 16(1). https://doi.org/10.1186/s13015-021-00204-6
Cheng, H., Wu, M., & Xu, Y. (2018). FMtree: a fast locating algorithm of FM-indexes for genomic data. Bioinformatics, 34(3), 416–424. https://doi.org/10.1093/bioinformatics/btx596
Cox, A. J., Bauer, M. J., Jakobi, T., & Rosone, G. (2012). Large-scale compression of genomic sequence databases with the Burrows-Wheeler transform. Bioinformatics, 28(11), 1415–1419. https://doi.org/10.1093/bioinformatics/bts173
Li, H. (2014). Fast construction of FM-index for long sequence reads. Bioinformatics, 30(22), 3274–3275. https://doi.org/10.1093/bioinformatics/btu541
'''