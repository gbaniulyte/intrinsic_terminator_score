Intrinsic terminator strenght (T<sub>S</sub>) calculator
===========================

This tool was designed to score potential intrinsic terminators as a batch input.
It was used to score RNA sequences identified by Term-seq in:

```
Philip P Adams, Gabriele Baniulyte, Caroline Esnault, Kavya Chegireddy, Navjot Singh, Molly Monge, Ryan K Dale, Gisela Storz, Joseph T Wade
(2021) Regulatory roles of Escherichia coli 5' UTR and ORF-internal RNAs detected by 3' end mapping.
eLife 10:e62438. https://doi.org/10.7554/eLife.62438
```
The scoring method was written as described in:

```
Chen, Ying-Ja, Peng Liu, Alec A. K. Nielsen, Jennifer A. N. Brophy, Kevin Clancy, Todd Peterson, and Christopher A. Voigt.
(2013) Characterization of 582 Natural and Synthetic Terminators and Quantification of Their Design Constraints. 
Nature Methods 10 (7): 659–64. https://doi.org/10.1038/nmeth.2515.
```

Requirements
-----
**Prerequisite:**
- Python 3.7 or higher
- Kinefold executable download from http://kinefold.curie.fr/download.html
- ViennaRNA (https://www.tbi.univie.ac.at/RNA/)
- Tested on WSL/Linux
```
	conda install bioconda::viennarna
```
- Input is a tab-delimited file with multiple sequences where:

```
	Header is optional
	column1 = genome coordinates or id
	column2 = strand (+ or -) (will be reverse-complemented if "-")
	column3 = DNA or RNA seqeunce (we used 60 nt uptream of the 3' end)
```
Usage
-----
1. Install ViennaRNA
2. Unpack "kinefold_long_static" folder somewhere conevenient e.g. where this script is
3. Run python executable
```
python int_term_score.py
```
4. Enter input text file path
```
./Input.txt
```
5. Enter output text file path
```
./Output.txt
```
6. Enter path to kinefold fodler with example.dat file
```
./kinefold_long_static/TEST
```
7. Enter if the input file has headers ("y" or "n")

Be patient, this will take long time. ~800 sequences takes about 40 min to process. 
Test it with ~10 sequences first to make sure it runs. Test input and output files are also provided.

Algorithm
---------

1. RNA seqeunces are passed through Kinefold (for co-transcriptional folding) 20 times to get the most frequent structure;
2. Then Kinefold output is converted into a dot-bracket annotation and sructures like the hairpin, loop, U-tract, etc. are identified;
3. This information is passed through RNAfold and RNAeval from ViennaRNA to get free energies for all parts of the RNA structure;
4. Finally, those values are used to calculate terminator strength (T<sub>S</sub>) using the method from Chen et al, 2013 paper;
5. An output text file with coordinates/id, structures, free energies and Ts score is generated.

Output
------

A detailed text file with the following columns:

- **3' end** - original sequence name/id/coordinate
- **Strand** - strand
- **3' end_seq** - RNA sequence used
- **Kinefold_structure** - Kinefold structure in a dot-bratcket annotation
- **A-tract** - predicted A-tract sequence as defined by Chen et al, 2013
- **Hairpin** - predicted hairpin sequence
- **Hp_structure** - predicted hairpin structure in a dot-bratcket annotation
- **Loop_seq** - predicted hairpin loop sequence
- **U-Tract** - predicted U-tract sequence
- **dGU** -  free energy of binding between the U-tract and the template DNA
- **dGL** - free energy for the closure of the hairpin loop 
- **dGH** - free energy of hairpin folding
- **dGHA** -  free energy of the folded RNA from 8nt upstream to 8nt downstream of the hairpin
- **dGA** - free energy of the extended hairpin
- **dGB** - free energy of the hairpin base 
- **TS** - instrinsic terminator strength/score

Refer to Chen et al, 2013 for details on RNA structure dissection, annotations and free energie calculations.

Caveats
-------

- If "kinefold_fail" appears, it's frequently because kinefold predicts the the 'molecular time' to be too short for the RNA to finish folding.
- If it fails at RNAfold or RNAeval then it's most likely due to RNA structure with multiple hairpins (doesn't happen with proper terminators) 
   or if the loop is <3nt, this will be obvious from the dot-bracket structure and a warning in the output file.
