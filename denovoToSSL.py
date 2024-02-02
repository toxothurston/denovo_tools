__author__ = 'rj8'
"""Input: novor.cloud denovo csv download.  Output: .ssl format file for making a blib spectral library.
1. Input the novor.cloud results file
2. If database result sequences available, determine FDR based on percent match between database/denovo 
3. Choose peptides and scan numbers based on mass accuracy, charge, peptide length, and FDR
4. Output a .ssl format file and a .fasta file

The fasta file can be created such that only denovo sequences are used that don't have a fasta match 
(denovo_only=True).  This can be appended to the fasta file to create a new file for additional database matches to
the denovo sequences.  If an ssl file is desired that contains all denovo sequences (even ones found in the fasta file)
then set this to False."""

import csv
from operator import itemgetter
from Bio import Align
aligner = Align.PairwiseAligner()

aligner.match_score = 1.0 #If pair of aa's in alignment are identical give it a score of 1
aligner.mismatch_score = 0 #If pair of aa's in alignment are not identical then score it as 0
aligner.gap_score = 0 #No penalties for gaps

#Parameters to set manually
input_filename = 'denovo_anais3ul.csv' #comes from novor.cloud
output_filename = input_filename[:-3] + 'ssl' #this plus the raw files are used to make a blib library
fasta_output = input_filename[:-3] + 'fasta' #need this to feed to Encyclopedia
experiment_name = 'NOVOR_NOVOR:anais' #just used to tag files and stuff
score_type = 'novor sequence score' #this is pretty meaningless, since skyline doesn't know what it is
score = 80 #The score can be specified here, or if 0 the score will be set based on the fdr
fdr_cutoff = 0.05 #Database and denovo sequences are compared and a score is set based on this fdr
fxn = 0.7 #Smith-Waterman alignment should have at least this much alignment agreement (denovo v datqbase seqs)
psm_fxn = 0.8 #Smith-Waterman alignment should have at least this much alignment agreement (denovo PSM to uniq)
max_mass_diff = 0.01 #max mass diff to be considered part of a PSM cluster of a single unique sequence (in Da)
scan_range_fxn = 20 #range of scan numbers divided by this number yields the max scan num diff between PSMs
find_unique = True #if True, only unique sequences are reported; otherwise, PSMs are reported
min_pep_length = 8 #peptides have to be at least this long to make it to the output
max_pep_length = 25 #peptides have to be no longer than this  to make it to the output
min_z = 2 #minimum charge to be in the output
max_z = 3 #maximum charge to be in the output
max_ppm_err = 10 #maximum precursor error to be allowed in the output
denovo_only = True #True returns only sequences that had no database match in novor.cloud, use False if there is comet_hits_filename
comet_hits_filename = '' #If not '', spectra that had Comet hits are excluded from the output files

def GetDenovoSequences():
    """GetDenovoSequences opens the csv output file from novovr.cloud, makes a list of dicts, where each dict is a
    denovo PSM that contains information like the scan number, denovo sequence score, and more.
        Args:
            none: the input file name is obtained from the header arguments
        Returns:
            sequences (list of dicts (N) containing N scans."""

    sequences = []
    with open(input_filename, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            sequence = {}
            sequence['file'] = row['Fraction'].strip()
            sequence['scan'] = int(row['Scan #'])
            sequence['mz'] = float(row['m/z'])
            sequence['z'] = int(row['z'])
            sequence['score'] = float(row['Score'])
            sequence['mw'] = float(row['Peptide Mass'])
            sequence['err_ppm'] = float(row['Error (ppm)'])
            sequence['length'] = int(row['Length'])
            sequence['denovo_ssl'] = row['De Novo Peptide'].replace('M(O)', 'M[+15.99]').replace('(Cam)', '[+57.00]').\
                replace('(Deamidated)', '[+0.984]').replace('(Pyro-Glu)', '[-18.015]')
            sequence['denovo'] = row['De Novo Peptide'].replace('(O)', '').replace('(Cam)', '').\
                replace('Q(Deamidated)', 'E').replace('N(Deamidated)', 'D').replace('(Pyro-Glu)', '')
            sequence['db'] = row['DB Sequence'].replace('I', 'L').replace('(O)', '').replace('(Cam)', '').\
                replace('Q(Deamidated)', 'E').replace('N(Deamidated)', 'D').replace('Pyro-Glu)', '')
            sequences.append(sequence)
        sorted_sequences = sorted(sequences, key=itemgetter('scan'), reverse=False)
    return sorted_sequences

def FDRtoScore(sequences):
    """FDRtoScore finds spectra for which both a database search match and denovo sequence exists.  The database match
    is presumed to be identified w/ high confidence (eg, FDR of 0.01).  The denovo score that corresponds to a desired
    denovo FDR is determined assuming the database sequence is correct.
            Args:
                sequences: list of dicts (N) containing N scans
            Returns:
                score (scalar containing the denovo sequence score that corresponds to fdr (the set FDR value desired
                for the denovo sequences)."""

    seqs_to_compare_lst = []

    #Get the sequences and make alignments
    for sequence in sequences:
        if sequence['denovo'] and sequence['db']:
            seqs_to_compare = {}
            seqs_to_compare['denovo'] = sequence['denovo']
            seqs_to_compare['db'] = sequence['db'].replace('I', 'L')
            align_score = aligner.score(seqs_to_compare['denovo'], seqs_to_compare['db'])
            seqs_to_compare['fxn_correct'] = align_score / len(seqs_to_compare['db'])
            if seqs_to_compare['fxn_correct'] > 1:
                seqs_to_compare['fxn_correct'] = 1.0
            seqs_to_compare['score'] = sequence['score']
            seqs_to_compare_lst.append(seqs_to_compare)
            sequence['fxn_correct'] = seqs_to_compare['fxn_correct']
        else:
            sequence['fxn_correct'] = -1

    #Sort the sequences by score (high to low), and then find FDR level and the corresponding score
    reverseSort_sequence_list = sorted(seqs_to_compare_lst, key=itemgetter('score'), reverse=True)
    fp = 0 #sum of false positives
    tp = 0 #sum of true positives
    for sequence in reverseSort_sequence_list:
        if sequence['fxn_correct'] >= fxn:
            tp += 1
        else:
            fp += 1
        sequence['fdr'] = fp / (fp + tp)
    forwardSort_sequence_list = sorted(reverseSort_sequence_list, key=itemgetter('score'), reverse=False)
    for i in range(len(forwardSort_sequence_list)):
        sequence_fdr = forwardSort_sequence_list[i]['fdr']
        sequence_score = forwardSort_sequence_list[i]['score']
        if sequence_fdr <= fdr_cutoff:
            score = sequence_score
            break

    print('Number of spectra that have both a denovo and database sequence: {0}'.format(len(forwardSort_sequence_list)))
    print('Number of true positives: {0}'.format(tp))
    print('Number of false positives: {0}'.format(fp))
    return score

def GetHighScoreSeqs(sequences, min_score):
    sorted_sequences = sorted(sequences, key=itemgetter('score'), reverse=True)
    high_score_sequences = []
    for sequence in sorted_sequences:
        if sequence['score'] >= min_score:
            if sequence['length'] <= max_pep_length and sequence['length'] >= min_pep_length:
                if sequence['z'] <= max_z and sequence['z'] >= min_z:
                    if sequence['err_ppm'] <= max_ppm_err:
                        high_score_sequences.append(sequence)
    return high_score_sequences

def FindUniqueSequences(sequences, min_score):
    """FindUniqueSequences tries to make a list of unique denovo sequences from the list of denovo PSMs
            Args:
                sequences: (list of dicts (N) containing N PSMs.
            Returns:
                uniq_seqs (list of dicts (M) containing M unique sequences."""

    #find range of scan numbers
    scans = sorted([d['scan'] for d in sequences if 'scan' in d])
    low = scans[0]
    high = scans[-1]
    max_scan_range = (high - low) / scan_range_fxn #scans plus/minus this range corresponds to fraction of the scan range

    #sort seqs by score high to low, then flag lower scoring seqs that seem to have similar seqs and mw
    sorted_sequences = sorted(sequences, key=itemgetter('score'), reverse=True)
    for i in range(len(sorted_sequences)):
        seq_i = sorted_sequences[i]
        if seq_i['score'] < min_score: break
        if seq_i['scan'] < 0: continue  # scan = -1 is flag to eliminate that seq
        for j in range(i + 1, len(sorted_sequences)):
            seq_j = sorted_sequences[j]
            if seq_j['score'] < min_score: break
            if seq_j['scan'] < 0: continue  # scan = -1 is flag to eliminate that seq
            seq_j = sorted_sequences[j]
            if seq_i['denovo'] == seq_j['denovo']:
                sorted_sequences[j]['scan'] = -1  # flag this PSM as being redundant and lower scoring
            else:
                scan_diff = abs(seq_i['scan'] - seq_j['scan'])
                if scan_diff < max_scan_range:
                    if seq_i['mw'] <= seq_j['mw']+max_mass_diff:
                        if seq_i['mw'] >= seq_j['mw']-max_mass_diff:
                            align_score = aligner.score(seq_i['denovo'], seq_j['denovo'])
                            align_fxn = align_score / len(seq_i['denovo'])
                            if align_fxn > psm_fxn:
                                sorted_sequences[j]['scan'] = -1 #flag this PSM as being redundant and lower scoring

    #Now that the redundant PSMs have been flagged, make the unique list of sequences
    uniq_seqs = []
    for sequence in sorted_sequences:
        if sequence['scan'] > 0:
            uniq_seqs.append(sequence)

    return uniq_seqs

def GetOnlyDenovos(high_score_sequences):
    denovos_only = []
    for sequence in high_score_sequences:
        if sequence['db'] == '':
            denovos_only.append(sequence)

    return denovos_only

def ReplaceDBHits(sequences):
    '''ReplaceDBHits is used when the Comet database search results are to be used instead of the Novor database
    search results.
    Args:
        sequences (list of dicts (N) containing N scans.
    Returns:
        sequences (list of dicts (N) containing N scans.'''

    #Get the comet hits
    comet_hits = []
    with open(comet_hits_filename, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            comet_hit = {}
            comet_hit['file'] = row['file'].strip()
            comet_hit['scan'] = int(row['scan'])
            comet_hit['sequence'] = row['sequence'].strip().replace('I', 'L')
            comet_hits.append(comet_hit)
    sorted_comet_hits = sorted(comet_hits, key=itemgetter('scan'), reverse=False)

    #remove the Novor hits
    for sequence in sequences:
        sequence['db_novor'] = sequence['db']
        sequence['db'] = ''

    #add the comet hits
    j_start = 0
    for i in range(len(sequences)):
        sequence = sequences[i]
        sequence_scan = sequences[i]['scan']
        for j in range(j_start, len(sorted_comet_hits)):
            comet_scan = sorted_comet_hits[j]['scan']
            if comet_scan > sequence_scan:
                j_start = j
                break
            elif comet_scan < sequence_scan:
                continue
            elif sorted_comet_hits[j]['file'] == sequences[i]['file']:
                sequences[i]['db'] = sorted_comet_hits[j]['sequence']

#debug
    '''keys = sequences[0].keys()
    with open('comet_added.txt', 'w', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, keys, dialect='excel-tab')
        dict_writer.writeheader()
        dict_writer.writerows(sequences)'''

    return(sequences)

def ExcludeCometHits(high_score_sequences):
    comet_hits = set()
    comet_sequences = set()
    denovo_only_sequences = []
    fasta_sequences = [] #for debugging
    with open(comet_hits_filename, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            comet_sequences.add(row['sequence'].strip().replace('I', 'L'))
            comet_hit = row['file'] + ' ' + str(row['scan'])
            comet_hits.add(comet_hit)
    for sequence in high_score_sequences:
        test_string = sequence['file'] + ' ' + str(sequence['scan'])
        if test_string not in comet_hits and sequence['denovo'] not in comet_sequences:
            denovo_only_sequences.append(sequence)
        else: #for debugging
            fasta_sequences.append(sequence)

    sorted_fasta_sequences = sorted(fasta_sequences, key=itemgetter('scan'), reverse=False) #for debugging


    return denovo_only_sequences

########################################################################################################################
########################################################################################################################

#Get the denovo sequence information from the input file
sequences = GetDenovoSequences()
print('Total number of PSMs, regardless of scores: {0}'.format(len(sequences)))

#Figure out which score corresponds to the desired FDR
if score == 0:
    score = FDRtoScore(sequences)
    print('For an FDR of {0}, the score cutoff is {1}'.format(fdr_cutoff, score))
elif score < 100 and score > 0:
    print('No FDR cutoff set, but the score was manually set to {0}'.format(score))
else:
    print('The score must be 0 or a number less than 100')
    exit(1)

#Remove the low scoring sequences below the FDR cutoff (or the specified score)
high_score_sequences = GetHighScoreSeqs(sequences, score)
print('Number of high scoring PSMs with proper z, mass error, and length: {0}'.format(len(high_score_sequences)))

#Remove sequences that have been identified in a fasta file from within novor.cloud, or by Comet
if denovo_only:
    high_score_sequences = GetOnlyDenovos(high_score_sequences)
    print('Only spectra that did not have a novor database match will be included in the output.')
    print('Number of denovo-only PSMs: {0}'.format(len(high_score_sequences)))
elif comet_hits_filename:
    high_score_sequences = ExcludeCometHits(high_score_sequences)
    print('Only spectra that did not have a Comet database match will be included in the output.')
    print('Number of denovo-only PSMs: {0}'.format(len(high_score_sequences)))
else:
    print('Spectra with database hits will be included in the output.')

#Find unique sequences
if find_unique:
    final_list = FindUniqueSequences(high_score_sequences, score)
    print('Number of unique sequences: {0}\n'.format(len(final_list)))
else:
    final_list = high_score_sequences

#Make the stupid fasta file for Encyclopedia' background or to append to an existing fasta file
header = '>' + experiment_name + ' ' + 'FDR=' + str(fdr_cutoff) + ' ' + 'score=' + str(score) + '\n'
seqs = ''.join([d['denovo'] for d in final_list if 'denovo' in d])
with open(fasta_output, "w") as text_file:
    text_file.write(header)
    text_file.write(seqs)

#Based on the score, create a new list of sequences to be used to make the ssl output file
ssl_sequences = []
for sequence in final_list:
    if sequence['score'] >= score:
        ssl_sequence = {}
        ssl_sequence['file'] = sequence['file']
        ssl_sequence['scan'] = sequence['scan']
        ssl_sequence['charge'] = sequence['z']
        ssl_sequence['sequence'] = sequence['denovo_ssl']
        ssl_sequence['score-type'] = score_type
        ssl_sequence['score'] = sequence['score']
        ssl_sequences.append(ssl_sequence)

keys = ssl_sequences[0].keys()
with open(output_filename, 'w', newline='') as output_file:
    dict_writer = csv.DictWriter(output_file, keys, dialect='excel-tab')
    dict_writer.writeheader()
    dict_writer.writerows(ssl_sequences)