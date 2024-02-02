from __future__ import division
__author__ = 'rj8'

"""This script looks at each spectrum query to determine if a Novor peptide is ranked #1.  If a database sequence is
also assigned a rank of 1, then the Novor ranking is deprecated to a rank of 2.  It also deals with the issue of when a
database sequence has an xcorr less than the de novo sequence, but its close enough to be considered correct."""

from lxml import objectify
from lxml import etree
import numpy as np

min_prec_z = 2 # only analyze spectra whose precursor charge state is as low as this
max_prec_z = 4 # only analyze spectra whose precursor charge state is as high as this
min_peptide_length = 8 # min peptide length to analyze
novor_fract = 0.99 # Based on decoy model, the percent of Novor peptides to capture w/ xcorr > fasta xcorr. Use 0.99

# Extract the spectra from lxml object into a list of dicts
def ExtractSpectraFromRoot(root):
    spectrums = []
    for e in root.msms_run_summary.getchildren():
        if 'spectrum_query' in e.tag:
            spectrum = {}
            spectrum['scan'] = int(e.attrib['start_scan'])
            spectrum['mw'] = float(e.attrib['precursor_neutral_mass'])
            spectrum['charge'] = int(e.attrib['assumed_charge'])
            spectrum['rt'] = float(e.attrib['retention_time_sec'])
            spectrum['index'] = int(e.attrib['index'])
            spectrum['matches'] = []
            for f in e.search_result.getchildren():
                match = {}
                match['rank'] = int(f.attrib['hit_rank'])
                match['peptide'] = f.attrib['peptide']
                match['protein'] = f.attrib['protein']
                for g in f.getchildren():
                    if 'search_score' in g.tag:
                        if g.attrib['name'] == 'xcorr':
                            match['xcorr'] = float(g.attrib['value'])
                spectrum['matches'].append(match)
            spectrums.append(spectrum)
    return spectrums

# Extract the xcorr score differences between the first and second highest DECOY hits
def GetDecoyDiffs(spectrums):
    decoy_diffs = []
    for spectrum in spectrums:
        if spectrum['charge'] < min_prec_z or spectrum['charge'] > max_prec_z: continue
        decoy_1 = -1
        decoy_2 = -1
        for match in spectrum['matches']:
            if 'DECOY' in match['protein'] and len(match['peptide']) >= min_peptide_length and decoy_1 == -1:
                decoy_1 = match['xcorr']
                continue
            if decoy_1 > 0:
                if 'DECOY' in match['protein'] and len(match['peptide']) >= min_peptide_length and decoy_2 == -1:
                    decoy_2 = match['xcorr']
                    break
        if decoy_1 > 0 and decoy_2 > 0:
            ratio = (decoy_1 - decoy_2) / spectrum['mw']
            decoy_diffs.append(ratio)

    return decoy_diffs #These are all normallized by their molecular weights

#From the decoy_diffs find a mw-normalized xcorr diff that captures novor_fract of the Novor peptides
def ModelFromDecoys(decoy_diffs):
    fract = 1 - novor_fract
    mw_norm_xcorr_diff = -1
    decoy_diffs.sort(reverse=True) #sort high to low
    percent_total = [i / len(decoy_diffs) for i in range(1, len(decoy_diffs)+1)]
    for i in range(len(percent_total)):
        if percent_total[i] > fract:
            mw_norm_xcorr_diff = decoy_diffs[i]
            break
    if mw_norm_xcorr_diff < 0:
        print ('DECOY modeling failed.  Most distressing!')
        exit(1)
    return mw_norm_xcorr_diff

#======================================================================================================================
#======================================================================================================================

# Get the list of pepxml file names to analyze
file_name_list = "input.txt" #each line is the name of a pep.xml file to be modified
# open the input file name list
with open(file_name_list, 'r') as InFile:
    lines = InFile.readlines()
    file_names = [x.strip() for x in lines]

# Loop through each pepxml file
for pepxml_file_name in file_names:
    print ('\nStarting on {0}'.format(pepxml_file_name))
    output_pepxml_file_name = pepxml_file_name[:-8] + '_unNovored.pep.xml'

    # Parse the xml file into a list of dicts
    with open(pepxml_file_name) as InFile:
        root = objectify.parse(InFile).getroot()
    spectrums = ExtractSpectraFromRoot(root)

    # Find mw-normalized xcorr differences between first and second ranked DECOY's
    decoy_diffs = GetDecoyDiffs(spectrums)
    print ('Number of decoy pairs: {0}'.format(str(len(decoy_diffs))))

    # Determine the xcorr diff that captures novor_fract of the Novor sequences that match the fasta matches
    mw_norm_xcorr_diff = ModelFromDecoys(decoy_diffs)

    # Now look through the pepxml file and change the rankings between Novor and fasta hits
    delta_xcorr_count = 0
    less_than_threshold_count = 0
    for e in root.msms_run_summary.getchildren():
        if 'spectrum_query' in e.tag:
            mw = float(e.attrib['precursor_neutral_mass'])
            xcorr_threshold = mw * mw_norm_xcorr_diff
            novor_xcorr = -1
            fasta_xcorr = -1
            try:
                #NOVOR_NOVOR is present, but not DECOY NOVOR
                if 'NOVOR_NOVOR' in e.search_result.search_hit[0].attrib['protein'] and \
                        'DECOY_nv' not in e.search_result.search_hit[0].attrib['protein']:
                    novor_xcorr = float(e.search_result.search_hit[0].search_score[0].attrib['value'])
                    """if e.search_result.search_hit[0].attrib['peptide_prev_aa'] != "K" and \
                        e.search_result.search_hit[0].attrib['peptide_prev_aa'] != "R":
                        novor_xcorr = 0 #mark this as being a non-tryptic cleavage of NOVOR_NOVOR which is BS
                    if e.search_result.search_hit[0].attrib['peptide'][-1] != "K" and \
                        e.search_result.search_hit[0].attrib['peptide'][-1] != "R":
                        novor_xcorr = 0"""
                    for i in range(1, len(e.search_result.search_hit[0])):
                        if 'NOVOR_NOVOR' not in e.search_result.search_hit[i].attrib['protein']:
                            if 'DECOY' not in e.search_result.search_hit[i].attrib['protein']:
                                fasta_xcorr = float(e.search_result.search_hit[i].search_score[0].attrib['value'])
                                fasta_index = i
                                break
            except:
                continue #probably an empty spectrum query

            # if denovo is first ranked and there is a lower ranked fasta match that is "close enough" then swap ranks
            if novor_xcorr >= 0 and fasta_xcorr >= 0:
                delta_xcorr = novor_xcorr - fasta_xcorr
                delta_xcorr_count += 1
                if delta_xcorr < xcorr_threshold:
                    less_than_threshold_count += 1
                    rank_swap = e.search_result.search_hit[fasta_index].attrib['hit_rank']
                    if rank_swap == '1':
                        rank_swap = '2' #If the denovo and fasta sequences are both ranked '1', need to force to '2'
                    e.search_result.search_hit[0].attrib['hit_rank'] = rank_swap
                    e.search_result.search_hit[fasta_index].attrib['hit_rank'] = "1"
    print ('Number of cases where NOVOR_NOVOR exceeded FASTA: {0}'.format(str(delta_xcorr_count)))
    print ('Number of cases where the FASTA was close enough to be re-ranked: {0}'.format(str(less_than_threshold_count)))
    et = etree.ElementTree(root)
    et.write(output_pepxml_file_name, pretty_print=True, xml_declaration=True, encoding="utf-8")







