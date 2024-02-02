# denovo_tools
Scripts to make use of de novo sequences.

denovoToSSL.py
>>aim: Using the novor.cloud denovo.csv output file, this script will create an ssl file suitable for making a denovo sequence blib.  It also will make a fasta file containing the concatenated denovo sequences.

>>input: Using novor.cloud one can download a file called denovo.csv.  If its from a single LCMS run then use it as is, but if more than one raw file is used then the first column has to be changed to reflect the raw file names (instead of F1, etc designations).

>>output: A tab delineated text file (.ssl) is created and also a FASTA file that can contain a single protein with all peptides appended or just those that were not identified in the novor.cloud database search.


modify_comet_pepxml.py
>>aim: If a Comet search is done using a fasta file that has denovo sequences appended as a concatenated 'protein', there is a competition between hits to a denovo sequence versus hits to the original fasta sequences.  Sometimes a denovo sequence is scored slightly
>>better than a FASTA-derived sequence even though the latter is correct.  This script locates such cases and reassigns their rankings in the pepxml file.
>>
>>input: A text file that lists the pepxml file names to analyze.
>>
>>output: Another pepxml file that is identical except for rankings 1 and 2 are reversed if the second ranked FASTA-derived sequence has a score very close to the first ranked denovo sequence.  This new pepxml file can then be subjected to peptideprophet in the TPP.
