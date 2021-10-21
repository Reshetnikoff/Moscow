# Moscow
Several programs for iGEM competition was created by Team Moscow

## Fuse IT

Aim of this program is fuse protein using linker (or without it) to this one. Base of this program module "Modeller" wrapped into user friendly shell. After fusing protein you can do whether you want! Optimize, minimize, modelling!

Input of this program file with any name of extention. This file should consists:

First line: name of file
Second line: name of pdb file, chain, name of file with sequence (three name separeted by comma!)
Third line: sequence of linker (can be emprty)
Fourt line: name of second pdb file, chain, name of file with sequence 
Etc...

Example:

Chimera
model.pdb,A,dcas9.seq
GGGGS
1ema.pdb,A,gfp.seq

### Usage

All program contain in "FuseIT directory"

python FuseIT [name of input file]

### Requirement:

- Python version 3 and higher
- Moddeller 

Result of the work of this program you can see in iGEM wiki modeling page:
https://2021.igem.org/Team:Moscow/Model

## Structure modeling with Chimera

We develope script which can create strcture of Cas9 Protein in double stranding DNA and single stranding DNA. We can change height of hairpin in out contruction and also length from second cas (which connect to ssDNA) to hairpin. 

Huge thanks to Roma Novikov who give script which served as a base for our. 

Result of the work of this program you can see in iGEM wiki modeling page:
https://2021.igem.org/Team:Moscow/Model
