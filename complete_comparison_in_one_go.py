#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from colorama import Fore, Style
import csv


def parse_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        identifier = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if identifier:
                    sequences[identifier] = ''.join(sequence).upper()
                identifier = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if identifier:
            sequences[identifier] = ''.join(sequence).upper()
    return sequences


# different replacement schemes, can be changed by commenting out/in
# BUT!!! needs to be also done in reverse unamsking with same scheme!!!!
def replace_amino_acids(sequences):
    # =============================================================================
    #     replacements = {
    #         'R': 'B', 'H': 'B', 'K': 'B',
    #         'D': 'J', 'E': 'J',
    #         'S': 'O', 'T': 'O', 'N': 'O', 'Q': 'O',
    #         'C': 'X', 'G': 'X', 'P': 'X',
    #         'A': 'Z', 'V': 'Z', 'I': 'Z', 'L': 'Z', 'M': 'Z', 'F': 'Z', 'Y': 'Z', 'W': 'Z'
    #     }
    # =============================================================================
    replacements = {
        'R': 'a', 'K': 'a',
        'H': 'b',
        'D': 'c', 'E': 'c',
        'S': 'd', 'T': 'd',
        'N': 'e', 'Q': 'e',
        'C': 'f',
        'G': 'g',
        'P': 'h',
        'A': 'i', 'V': 'i', 'I': 'i', 'L': 'i',
        'M': 'j',
        'F': 'k', 'Y': 'k', 'W': 'k'
    }
# =============================================================================
#     replacements = {
#         'R': 'a',
#         'K': 'l',
#         'H': 'b',
#         'D': 'c', 'E': 'c',
#         'S': 'd', 'T': 'd',
#         'N': 'e', 'Q': 'e',
#         'C': 'f',
#         'G': 'g',
#         'P': 'h',
#         'A': 'i', 'V': 'i', 'I': 'i', 'L': 'i',
#         'M': 'j',
#         'F': 'k', 'Y': 'k', 'W': 'k'
#     }
# =============================================================================
# =============================================================================
#     replacements = {
#         'R': 'a', 'K': 'a',
#         'H': 'b',
#         'D': 'c', 'E': 'c',
#         'S': 'd', 'T': 'd',
#         'N': 'e', 'Q': 'e',
#         'C': 'f',
#         'G': 'g',
#         'P': 'h',
#         'A': 'i', 'V': 'i', 'I': 'i', 'i': 'i',
#         'M': 'j',
#         'F': 'k',
#         'Y': 'l',
#         'W': 'm'
#     }
# =============================================================================
    refined_sequences = {}
    for identifier, sequence in sequences.items():
        refined_sequence = ''.join(replacements.get(
            residue, residue) for residue in sequence)
        refined_sequences[identifier] = refined_sequence
    return refined_sequences


def get_user_input_for_comparison(sequences):
    print(Fore.BLUE + "Available sequences:")
    print(Style.RESET_ALL)
    check = ""
    for key in sequences.keys():
        print(key)
    entries = input(
        Fore.BLUE + "Please enter the sequence identifiers you want to compare, separated by commas: " + Style.RESET_ALL)
    entry_list = [entry.strip() for entry in entries.split(',')]
    while(check == ""):
        counter = 0
        for entry in entry_list:
            if entry.strip() in sequences:
                counter += 1
            else:
                entries = input(
                    Fore.BLUE + "Unvalid Identifier. Please provide identifiers from the list." + Style.RESET_ALL)
                entry_list = [entry.strip() for entry in entries.split(',')]
                break
        if(counter == len(entry_list)):
            check = "moin"
        counter = 0
    return entry_list


def create_reference_sequence(sequences, entry_list):
    min_length = min(len(sequences[entry]) for entry in entry_list)
    reference_sequence = {'symbol': [], 'index': []}
    for i in range(min_length):
        current_symbols = [sequences[entry][i] for entry in entry_list]
        if all(symbol == current_symbols[0] for symbol in current_symbols):
            reference_sequence['symbol'].append(current_symbols[0])
            reference_sequence['index'].append(i)
    reference_sequence_filtered = {
        'symbol': [symbol for symbol in reference_sequence['symbol'] if symbol != '-'],
        'index': [reference_sequence['index'][i] for i, symbol in enumerate(reference_sequence['symbol']) if symbol != '-']
    }
    return reference_sequence_filtered


def define_motif_length():
    length = ""
    while length == "" or not isinstance(length, int):
        length = input(
            Fore.BLUE + "Please enter the minimum amount of amino acids to form a motif: " + Style.RESET_ALL)
        try:
            length = int(length)
            print(Fore.GREEN + "The motif length has been set to:", length)
            print(Style.RESET_ALL)
        except ValueError:
            print(Fore.RED + "Unvalid answer, please provide only integers!")
            print(Style.RESET_ALL)
            continue
    return length


def refine_reference_sequence(reference_sequence, length):
    indices = reference_sequence['index']
    symbols = reference_sequence['symbol']
    filtered_indices = []
    filtered_symbols = []
    i = 0
    counter = 0
    while i <= len(indices) - length:
        while indices[i + 1] == indices[i] + 1 and i <= len(indices) - length:
            i += 1
            counter += 1

        if (counter >= length):
            for j in range(counter):
                filtered_indices.append(indices[i - counter + j])
                filtered_symbols.append(symbols[i - counter + j])
        else:
            i += 1
        counter = 0
    motive_reference_sequence = {
        'symbol': filtered_symbols,
        'index': filtered_indices
    }
    return motive_reference_sequence


def create_table_from_string(s):
    return {
        'symbol': list(s),
        'index': list(range(len(s)))
    }


def filter_table_by_reference(table, reference_indices):
    filtered_table = {
        'symbol': [],
        'index': []
    }
    for symbol, index in zip(table['symbol'], table['index']):
        if index in reference_indices:
            filtered_table['symbol'].append(symbol)
            filtered_table['index'].append(index)
    return filtered_table


def compare_tables(reference_sequence, new_table, length):
    ref_indices = reference_sequence['index']
    ref_symbols = reference_sequence['symbol']
    new_indices = new_table['index']
    new_symbols = new_table['symbol']
    i = 0
    counter = 0
    while i <= len(new_indices) - length:
        while ref_indices[i + 1] == ref_indices[i] + 1 and i <= len(new_indices) - length:
            counter += 1
            i += 1
        for j in range(counter):
            if(ref_symbols[i - counter + j] == new_symbols[i - counter + j]):
                if(j == counter):
                    del ref_indices[i:i+counter]
                    del ref_symbols[i:i+counter]
            else:
                break
        if (counter == 0):
            i += 1
        counter = 0
    return {'symbol': ref_symbols, 'index': ref_indices}


def unmask(filename):
    # =============================================================================
    #     replacements = {
    #         'B': 'R H K',
    #         'J': 'D E',
    #         'O': 'S T N Q',
    #         'X': 'C G P',
    #         'Z': 'A V I L M F Y W'
    #     }
    # =============================================================================
    replacements = {
        'a': 'R K',
        'b': 'H',
        'c': 'D E',
        'd': 'S T',
        'e': 'N Q',
        'f': 'C',
        'g': 'G',
        'h': 'P',
        'i': 'A V I L',
        'j': 'M',
        'k': 'F Y W'
    }
# =============================================================================
#     replacements = {
#         'a': 'R',
#         'l': 'K',
#         'b': 'H',
#         'c': 'D E',
#         'd': 'S T',
#         'e': 'N Q',
#         'f': 'C',
#         'g': 'G',
#         'h': 'P',
#         'i': 'A V I L',
#         'j': 'M',
#         'k': 'F Y W'
#     }
# =============================================================================
# =============================================================================
#     replacements = {
#         'a': 'R K',
#         'b': 'H',
#         'c': 'D E',
#         'd': 'S T',
#         'e': 'N Q',
#         'f': 'C',
#         'g': 'G',
#         'h': 'P',
#         'i': 'A V I L',
#         'j': 'M'
#         'k': 'F',
#         'l': 'Y',
#         'm': 'W'
#     }
# =============================================================================
    with open(filename, mode='r', newline='', encoding='utf-8') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        data = list(reader)

    for row in data:
        if 'symbol' in row:
            original_symbol = row['symbol']
            for key, value in replacements.items():
                if key in original_symbol:
                    original_symbol = original_symbol.replace(key, value)
            row['symbol'] = original_symbol

    with open(filename, mode='w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

    return


def main():
    file_path = input(
        Fore.BLUE + "Please enter the path to your FASTA file: " + Style.RESET_ALL)
    sequences = parse_fasta(file_path)
    sequences = replace_amino_acids(sequences)
    entry_list = get_user_input_for_comparison(sequences)
    reference_sequence = create_reference_sequence(sequences, entry_list)
    length = define_motif_length()
    if reference_sequence:
        reference_sequence = refine_reference_sequence(
            reference_sequence, length)
        remaining_sequences = {key: val for key,
                               val in sequences.items() if key not in entry_list}
        for key in remaining_sequences.keys():
            new_table = create_table_from_string(
                remaining_sequences[key])
            new_table = filter_table_by_reference(
                new_table, reference_sequence['index'])
            reference_sequence = compare_tables(
                reference_sequence, new_table, length)

    filename = input(
        Fore.BLUE + "Please provide a name for the output file: " + Style.RESET_ALL)
    with open(filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)

        csvwriter.writerow(['index', 'symbol'])

        index = reference_sequence.get("index")
        value = reference_sequence.get("symbol")

        for element in range(len(index) - length):
            if(index[element + 1] == index[element] + 1):
                csvwriter.writerow([index[element], value[element]])
            else:
                csvwriter.writerow([index[element], value[element]])
                csvwriter.writerow(["", ""])
    unmask(filename)

    print(Fore.GREEN + "Output file", filename, "created successfully.")
    print(Style.RESET_ALL)


if __name__ == "__main__":
    main()
