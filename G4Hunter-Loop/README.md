# G4Assassin
JupyterNotebook to perform handsfree G4Hunter analysis on libraries of genomes. Original G4Hunter software is https://github.com/AnimaTardeb/G4Hunter

# Required libraries
- BioPython
- Matplotlib
- Numpy
- Pandas
- os
- subprocess
- glob

# Instructions
1. Rename GENOMES variable as the file path to your library of genomes. This will search all subfolders of varying depths, and save all directories of every FASTA files within the directory
```
# For example
GENOME = 'Genomes/'
# This is the filepath to a folder within the same folder as G4Assassin.ipynb
```
2. After the code has run, the results will be saved in the generated Results folder. Within this folder, each subfolder, containing the results for a specific genome, will be named with 'Results_OriginalFileName'
3. Within said folder, will be the raw data, merged data (the most useful) detailed with a '-Merged.txt' at the end of the file name, and the overall statistics such as Number of Putative Quadruplex Sequences, Number of GCs, Percentage of GCs, among others.
