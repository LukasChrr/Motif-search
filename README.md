# Motif-search
Python script to automise search for conserved motifs in amino acid sequences


How to run script:
Open terminal 
Navigate to folder that contains Motif_identification.py
Run with command: python Motif_identification.py

Allignment input file:
Provide an allignment created from all amino acid sequences of interest
Needs to be a FASTA file

Change of ambiguity codes:
Copy these into the function „replace_amino_acids“:
1
replacements = {
    'R': 'B', 'H': 'B', 'K': 'B',
    'D': 'J', 'E': 'J',
    'S': 'O', 'T': 'O', 'N': 'O', 'Q': 'O',
    'C': 'X', 'G': 'X', 'P': 'X',
    'A': 'Z', 'V': 'Z', 'I': 'Z', 'L': 'Z', 'M': 'Z', 'F': 'Z', 'Y': 'Z', 'W': 'Z'
}

2
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

3
replacements = {
    'R': 'a',
    'K': 'l',
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

4
replacements = {
    'R': 'a', 'K': 'a',
    'H': 'b',
    'D': 'c', 'E': 'c',
    'S': 'd', 'T': 'd',
    'N': 'e', 'Q': 'e',
    'C': 'f',
    'G': 'g',
    'P': 'h',
    'A': 'i', 'V': 'i', 'I': 'i', 'i': 'i',
    'M': 'j',
    'F': 'k',
    'Y': 'l',
    'W': 'm'
}


Copy the corresponding one (same #) into the function „unmask“:
1
replacements = {
    'B': 'R H K',
    'J': 'D E',
    'O': 'S T N Q',
    'X': 'C G P',
    'Z': 'A V I L M F Y W'
}
2
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
3
replacements = {
    'a': 'R',
    'l': 'K',
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
4
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
    'j': 'M'
    'k': 'F',
    'l': 'Y',
    'm': 'W'
}

User input reference sequence:
Provide identifiers (headers in alignment) of sequences for creation of reference sequence;
Need to be separated by commas

User input motif length:
Provide an integer 1 <= input <= sequence length

User input file name: 
Provide name for output file in the following way:
name_of_choice.csv (always provide .csv for file type!)
File is saved in the directory where the Python script is saved

Written in IPython 8.15.0; in Spyder IDE version 5.4.3
