from SearchADT import *
import sys

pattern, sequence, bwt, suffix_array, c_table, occ_ranks = None, None, None, None, None, None
mutation_count = 0

def main():

    print("Welcome to the DNA Sequencer")
    while True:
        choice = input("Choose an option:\n1: Choose Protein From Library\n2: Add Your Own\n> ")
        if choice == '1':
            choose_from_library()
        elif choice == '2':
            add_your_own()
        elif choice == '3':
            print("Exiting program.")
            return
        else:
            print("Invalid choice. Please try again.")

def check_against_existing(): #O(n**2)
    # The below code performs an automatic internal search to see if the entered sequence was either found in the library in any protein.
    global pattern, sequence, bwt, suffix_array, c_table, occ_ranks
    
    library = ['Hemoglobin','Insulin','Collagen','Actin','Myosin','Fibrinogen','Melanin','Glucagon','Oxytocin','Immunoglobulin G']
    found = []
    if len(sequence) == 189:
        for pro in library:
            seq = csv_reader(pro)
            if seq == sequence:
                print()
                print(f'The entered sequence is the same as {pro} in our library!')
    
    elif len(sequence) > 189:   # Here it is looking for whether the larger sequence has any of the library sequences as a subseqeunce
        for pro in library:
            seq = csv_reader(pro)
            if seq in sequence:
                found.append(pro)

        if len(found) > 1:
            print()
            print(f'These DNA sequences from our library were found in your entered gene sequence: {found}')
            proceed = input('Do you want to see at what positions any one of them was found at or proceed with normal search?' +'\n' + 'Enter "Y" for find or "N" for normal search: ').upper().strip()
            while proceed != 'Y' and proceed != 'N':
                print()
                print('Invalid entry, try again!')
                proceed = input('Do you want to see at what positions any one of them was found at or proceed with normal search?' +'\n' + 'Enter "Y" for find or "N" for normal search: ').upper().strip()
            if proceed == 'Y':
                if len(found) > 1:
                    for v in range(len(found)):
                        print(f'{v+1}: {found[v]}')
                    pro = int(input('Enter the number for the protein that you want to find: '))
                    while pro < 1 or pro > v + 1:
                        print()
                        print('Invalid entry, try again!')
                        pro = int(input('Enter the number for the protein that you want to find: '))
                    name = found[pro-1]
                    pro = csv_reader(name)
                    bwt, suffix_array = bwt_and_suffix_array(sequence)
                    c_table = create_c_table(bwt)
                    occ_rank = rank_bwt(bwt)
                    coordinates = search_with_fm_index(pro, bwt, c_table, occ_rank, suffix_array)
                    print()
                    print(f'{name} was found at {coordinates} in the sequence')
                    print()
                    print('Now we proceed to normal search')

        elif len(found) == 1:
            pro = csv_reader(found[0])
            bwt, suffix_array = bwt_and_suffix_array(sequence)
            c_table = create_c_table(bwt)
            occ_rank = rank_bwt(bwt)
            coordinates = search_with_fm_index(pro, bwt, c_table, occ_rank, suffix_array)
            print()
            print(f'{found[0]} was found at {coordinates} in the sequence')
            print()
            print('Now we proceed to normal search')
    

    else:   # Here it is checking if the entered sequence is a subsequence for any library seqeunce.
        for pro in library:
            seq = csv_reader(pro)
            if sequence in seq:
                found.append(pro)

        if len(found) > 1:
            print(f'These DNA sequences from our library were found to contain your entered gene sequence: {found}')
            proceed = input('Do you want to see at what positions in any one of them it was found at or proceed with normal search?' +'\n' + 'Enter "Y" for find or "N" for normal search": ').upper().strip()
            while proceed != 'Y' and proceed != 'N':
                print('Invalid entry, try again!')
                proceed = input('Do you want to see at what positions in any one of them it was found at or proceed with normal search?' +'\n' + 'Enter "Y" for find or "N" for normal search": ').upper().strip()
            
            if proceed == 'Y':
                for v in range(len(found)):
                    print(f'{v+1}: {found[v]}')
                pro = int(input('Enter the number for the protein where you want to find its positions: '))
                while pro < 1 or pro > v + 1:
                    print('Invalid entry, try again!')
                    pro = int(input('Enter the number for the protein where you want to find its positions: '))
                name = found[pro-1]
                pro = csv_reader(name)
                bwt, suffix_array = bwt_and_suffix_array(pro)
                c_table = create_c_table(bwt)
                occ_rank = rank_bwt(bwt)
                coordinates = search_with_fm_index(sequence, bwt, c_table, occ_rank, suffix_array)
                print()
                print(f'Your sequence was found at {coordinates} in {name}')
                print()
                print('Now we proceed to normal search')
            else:
                pass

        elif len(found) == 1:
            pro = csv_reader(found[0])
            bwt, suffix_array = bwt_and_suffix_array(pro)
            c_table = create_c_table(bwt)
            occ_rank = rank_bwt(bwt)
            coordinates = search_with_fm_index(sequence, bwt, c_table, occ_rank, suffix_array)

            print()
            print(f'Your sequence was found at {coordinates} in {found[0]}')
            print()
            print('Now we proceed to normal search')
    process_sequence()
   

def choose_from_library():

    global pattern, sequence, bwt, suffix_array, c_table, occ_ranks

    proteins = {
        '1': "Hemoglobin",
        '2': "Insulin",
        '3': "Collagen",
        '4': "Actin",
        '5': "Myosin",
        '6': "Fibrinogen",
        '7': "Melanin",
        '8': "Glucagon",
        '9': "Oxytocin",
        '10': "Immunoglobulin G"
    }

    print("Choose any one protein from the existing library:")
    for key, value in proteins.items():
        print(f"{key}: {value}")
    selection = input("> ")
    while selection not in proteins:
        print("Invalid selection. Try again!")
        selection = input("> ")
    sequence = csv_reader(proteins[selection])
    process_sequence()


def add_your_own():

    global pattern, sequence, bwt, suffix_array, c_table, occ_ranks

    type = input("Enter 'T' for Template or 'C' for Coding Strand: ").upper().strip()
    while type != 'T' and type != 'C':
        print('Invalid entry! Try again')
        type = input("Enter 'T' for Template or 'C' for Coding Strand: ").upper().strip()

    sequence = input("Enter your gene sequence (Template Strand or Coding Strand): ").upper().strip()
    while validation_check(sequence) != "Processing..." or len(sequence) < 3:
        print('Invalid entry, try again')
        sequence = input("Enter your gene sequence (Template Strand or Coding Strand): ").upper().strip()
        
    if type == 'C':
        sequence = coding_to_template(sequence)
    check_against_existing()

def process_sequence():
    global pattern, sequence, bwt, suffix_array, c_table, occ_ranks

    result = gene_translation(sequence)
    # Check the length of the tuple returned.
    if len(result) == 3:
        polypeptide, _ , error = result
        if error:
            print(error)
            print(f"Polypeptide Chain: \n{polypeptide}")
    elif len(result) == 2:
        polypeptide, _ = result
        print(f"Polypeptide Chain: {polypeptide}")
    
    search_pattern()


def search_pattern(S = None):
    global pattern, sequence, bwt, suffix_array, c_table, occ_ranks, mutation_count

    if S is None:
        pattern = input("Enter a pattern to search for: ").upper().strip()
        while validation_check(pattern) != 'Processing...' or len(pattern) > len(sequence) or len(pattern) < 3:
            print('Pattern too long/short or of incorrect form, try again!')
            pattern = input("Enter a pattern to search for: ").upper().strip()

    bwt, suffix_array = bwt_and_suffix_array(sequence)
    c_table = create_c_table(bwt)
    occ_ranks = rank_bwt(bwt)
    coordinates = search_with_fm_index(pattern, bwt, c_table, occ_ranks, suffix_array)
    while coordinates == -1:
        c_or_e = input("Pattern not found. Try another pattern or exit. \n Enter 'C' for another pattern \n Enter 'E' to exit \n > ").upper().strip()
        while c_or_e != 'C' and c_or_e != 'E':
            print('Invalid entry, try again!')
            c_or_e = input("Pattern not found. Try another pattern or exit. \n Enter 'C' for another pattern \n Enter 'E' to exit \n > ").upper().strip()
        if c_or_e == 'C':
            pattern = input("Enter a pattern to search for: ").upper().strip()
            while validation_check(pattern) != 'Processing...' or len(pattern) > len(sequence):
                print('Pattern too long or incorrect form, try again!')
                pattern = input("Enter a pattern to search for: ").upper().strip()
            coordinates = search_with_fm_index(pattern, bwt, c_table, occ_ranks, suffix_array)
        else:
            print('END')
            sys.exit()

    print(f"Pattern found at positions: {coordinates}")
    print("You can perform mutations at the following coordinates:", coordinates)
    print()
    coords = lst_to_coord(coordinates)
    mutate_sequence(coords)

def lst_to_coord(coordinates): #O(n)

    if len(coordinates) > 1:
        is_tuple = True
        try:
            coord = input('Enter one of the coordinates to mutate: ').strip()
            a, b = map(int, coord.strip().split(",")) #O(n)
            coord = tuple((a, b))
            is_tuple = True
        except:
            is_tuple = False
        while coord not in coordinates or is_tuple == False:
            print('Invalid entry, try again! Enter in the form "a, b"')
            print()
            try:
                coord = input('Enter one of the coordinates to mutate: ')
                a, b = map(int, coord.strip().split(","))
                coord = tuple((a, b))
                is_tuple = True
            except:
                is_tuple = False
        return coord
    else:
        coord = coordinates[0]
        return coord

def mutate_sequence(coordinates):
    global pattern, sequence, bwt, suffix_array, c_table, occ_ranks, mutation_count

    
    mut_type = input("Choose mutation type (Insertion, Substitution, Deletion): ").capitalize().strip()

    while mut_type not in ['Insertion', 'Substitution', 'Deletion']:
        print('Incorrect entry, try again!')
        mut_type = input("Choose mutation type (Insertion, Substitution, Deletion): ").capitalize().strip()
        
    if mut_type in ["Insertion", "Substitution"]:
        sub_sequence = input("Enter the subsequence for mutation: ").upper().strip()

        while validation_check(sub_sequence) != "Processing...":
            print('Sub-seqeunce invalid!')
            sub_sequence = input("Enter the subsequence for mutation: ").upper().strip()
        sequence, bwt, suffix_array, c_table, occ_ranks = mutations(sequence, mut_type, coordinates, sub_sequence)
    elif mut_type == "Deletion":
        sequence, bwt, suffix_array, c_table, occ_ranks = mutations(sequence, mut_type, coordinates)

    translated_gene = gene_translation(sequence)

    if len(translated_gene) == 3:
        polypeptide, _, nonsense = translated_gene
        print(nonsense)
        print()
        print(f"Updated Polypeptide Chain: \n {polypeptide}")
    else:
        print('Updated Polypeptide Chain:')
        print(translated_gene[0])
    
    if mutation_count < 2:
        choice = input('Do you want to perform another mutation? Enter "Y" or "N": ').upper().strip()
        while choice != "Y" and choice != 'N':
            print('Invalid entry, try again!')
            choice = input('Do you want to perform another mutation? Enter "Y" or "N": ').upper().strip()
            
        if choice == 'Y':
            mutation_count += 1
            search_again()
        else:
            print('END')
            sys.exit()
    else:
        print(f'Your final sequence looks like this: \n {sequence}')
        sys.exit()


def search_again():
    global pattern, sequence, bwt, suffix_array, c_table, occ_ranks, mutation_count

    print()
    print(mutation_count, 'mutation(s) complete!')
    pt = input('Choose either one of the options:\n 1) Search with same pattern. \n 2) Search with new pattern \n >').strip()
    while pt not in ['1', '2']:
        print('Invalid entry, try again!')
        pt = input('Choose either one of the options:\n 1) Search with same pattern. \n 2) Search with new pattern \n >').strip()
    if pt == '1':
        search_pattern('1')
    else:
        search_pattern()


if __name__ == "__main__":
    main()

# Overall Time Complexity O(n**2)
# Overall Space Complexity O(n)
