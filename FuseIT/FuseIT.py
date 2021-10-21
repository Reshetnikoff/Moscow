import os

from modeller import *
from modeller.automodel import *


amino_acids_dict = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'ASX': 'B',
    'CYS': 'C',
    'GLU': 'E',
    'GLY': 'G',
    'GLN': 'Q',
    'QLX': 'Z',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    "TRP": 'W',
    'TYR': 'Y',
    'VAL': 'V'
}

appropriate_amino_acids = list(amino_acids_dict.values())

# This beta version of function which give sequence from pdb file
# But it's not work well with MODELLER. This function and it give 
# difference sequences

def parse_pdb(pdb_file, chain):
    try:
        f = open(pdb_file, 'r')
    except FileNotFoundError:
        print(f'File {pdb_file} does not exist')
    
    f = open(pdb_file, 'r')    
    pdb_data = [x[:-1].split() for x in f.readlines()]
    f.close()
    
    pdb_data_seq = [x for x in pdb_data if (x[0] == 'SEQRES')]
    pdb_data_seq = [x for x in pdb_data_seq if (x[2] == chain)]
    pdb_data_seq = [x[4:] for x in pdb_data_seq if (x[2] == chain)]
    
    seq = [[amino_acids_dict[i] for i in x if i in amino_acids_dict.keys()] for x in pdb_data_seq]
    seq = "".join(["".join([j for j in i]) for i in seq])
    return seq

# Function check that sequence constist of correct amino acids
def check_amino_acids_seq(seq, comment):
    for aa in seq:
        if aa not in appropriate_amino_acids:
            raise Exception(f"Sequence in {comment} contain appropriate symbols (\'{aa}\')")

# Return seq (also check it). If seq file is
# pdb_file - pdb file
# chain - name of chain
# seq file - file with sequence in  one line 
def return_seq(pdb_file, chain, seq_file):
    if seq_file != 'None':
        try:
            f = open(seq_file, 'r')
        except FileNotFoundError:
            print(f'File {seq_file} does not exist')
            raise Exception(FileNotFoundError)
        seq = f.read()
        seq = seq.rstrip()
        seq = seq.upper()
        check_amino_acids_seq(seq, seq_file)
        f.close()
    else:
        seq = parse_pdb(pdb_file, chain)
        if len(seq) == 0:
            raise Exception("Don't find seq in pdb. Maybe you make \
a mistake in chain name (A, B, C and etc - is usual name for it)")
    return seq

# Just check that file exists
def check_pdb(file):
    try:
        f = open(file, 'r')
    except FileNotFoundError:
        print(f'File {file} does not exist')
        raise Exception(FileNotFoundError)
    finally:
        f.close()

# Parse input file. Check that all is correct and get needed data
def parse_input(input_file):
    try:
        f = open(input_file, 'r')
    except FileNotFoundError:
        print(f'File {input_file} does not exist')
        raise Exception(FileNotFoundError)
    
    input_data = [x.rstrip() for x in f.readlines()]
    name = input_data[0]
    input_data = input_data[1:]
    f.close()
    if len(input_data) < 3:
        raise Exception("Lines of input file is less than 3")
    
    pdb_files = []
    sequences = []
    linkers = []
    chains = []
    
    for i in range(len(input_data)):
        input_one_line = input_data[i]
        if i % 2 == 0:
            input_one_line = input_one_line.split(',')
            if len(input_one_line) != 3: raise Exception("Description of pdb must be 3 name separated by comma")
            pdb_file, chain, seq_file = input_one_line
            seq = return_seq(pdb_file, chain, seq_file)
            
            check_pdb(pdb_file)
            pdb_files.append(pdb_file)
            sequences.append(seq)
            chains.append(chain)
        else:
            seq = input_one_line
            check_amino_acids_seq(seq, f"{int(i/2+1)} linker from input file")
            linkers.append(seq)
    if len(linkers) != len(pdb_files) - 1: raise Exception("Appropriate input format \
(numbers of linkers is not equal numbers of pdb files - 1)")
            
    return name, pdb_files, chains, sequences, linkers

# Create ali file that is needed for working modeller
# Input: 
# name - name of output pdb file and ali file
# pdb_files, chains, sequences and linkers is returned from parse input function
# Result of work - ali file (aligment file)

def create_ali_file(name, pdb_files, chains, sequences, linkers):
    f = open(f'{name}.ali', 'w')
    
    for i in range(len(pdb_files)):
        pdb_name = pdb_files[i].split('.')[0]
        chain = chains[i]
        sequence = sequences[i]

        f.write(f">P1;{pdb_name}.{chain}\n")
        f.write(f"StructureX:{pdb_name}::{chain}::::::\n")

        start_ = sum([len(x) for x in sequences[:i]]) + sum([len(x) for x in linkers[:i]])
        if len(pdb_files) == i+1:
            end_ = 0
        else:
            end_ = sum([len(x) for x in sequences[i+1:]]) + sum([len(x) for x in linkers[i:]])

        f.write(f"{'-'*start_}{sequence}{'-'*end_}*\n\n")
        
    f.write(f">P1;{name}\n")
    f.write(f"sequence:{name}::::::::\n")
    
    n_components = len(pdb_files) + len(linkers)
    for i in range(len(pdb_files)):
        if i + 1 == len(pdb_files):
            f.write(f"{sequences[i]}*\n")
        else:
            f.write(f"{sequences[i]}{linkers[i]}")
    f.close()

# General function that fuse protein
# name - name of output file (needed to equal name from create_ali_file
# pdb_files - pdb files ;)
# chains - names of chain
def fuse(name, pdb_files, chains):
    env = Environ()

    # env.io.hetatm = False

    knowns = [f"{i.split('.')[0]}.{j}"for i, j in zip(pdb_files, chains)]

    a = AutoModel(env, alnfile=f'{name}.ali',
                  knowns=knowns,
                  sequence=name, root_name=name)
    a.starting_model = 1
    a.ending_model = 1
    a.make()


if __name__ == "__main__":
    print("Parsing input...")
    name, pdb_files, chains, sequences, linkers = parse_input('input.txt')
    print("Create ali file...")
    create_ali_file(name, pdb_files, chains, sequences, linkers)
    print("Run modeller and fuse the proteins!")
    fuse(name, pdb_files, chains)

    print(f"Done! Your output file named {name}..B99990001.pdb!")
