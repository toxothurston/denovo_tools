# denovo_tools
Scripts to make use of de novo sequences.

denovoToSSL.py
>>aim: Using the novor.cloud denovo.csv output file, this script will create an ssl file suitable for making a denovo sequence blib.  It also will make a fasta file containing the concatenated denovo sequences.

>>input: Using novor.cloud one can download a file called denovo.csv.  If its from a single LCMS run then use it as is, but if more than one raw file is used then the first column has to be changed to reflect the raw file names (instead of F1, etc designations).

>>output: A tab delineated text file (.ssl) is created and also a FASTA file that can contain a single protein with all peptides appended or just those that were not identified in the novor.cloud database search.
