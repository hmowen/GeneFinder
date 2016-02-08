# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: HARPER OWEN

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
	""" Returns the complementary nucleotide

	    nucleotide: a nucleotide (A, C, G, or T) represented as a string
	    returns: the complementary nucleotide
	>>> get_complement('A')
	'T'
	>>> get_complement('C')
	'G'
	>>> get_complement('R')

	>>> get_complement(2)
	"""
	# TODO: implement this
	letter = str(nucleotide)																# set letter = parameter (make sure it's a string)
	if letter == 'A':																		# check if letter is A
		return 'T'																			# return T
	elif letter == 'T':																		# check if letter is T
		return 'A'																			# return A
	elif letter == 'G':																		# check if letter is G
		return 'C'																			# return C
	elif letter == 'C':																		# check if letter is C
		return 'G'																			# return G
	else:
		return None

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("AAAAAAAA")
    'TTTTTTTT'
    """
    rdna = dna[::-1]																		# reverses input

    rev_dna = ""																			# initializes empty string
    index = 0																				# intitialize index
    while index < len(rdna):																# while loop, ends at len(dna)-1
    	reverse_letter = get_complement(rdna[index])										# gets the complement for the string
    	rev_dna = rev_dna + reverse_letter													# adds the new letter to the string
    	index += 1																			# indexes up 1
    return rev_dna 																			# returns string


def rest_of_ORF(dna):
	""" Takes a DNA sequence that is assumed to begin with a start
	codon and returns the sequence up to but not including the
	first in frame stop codon.  If there is no in frame stop codon,
	returns the whole string.

	dna: a DNA sequence
	returns: the open reading frame represented as a string
	>>> rest_of_ORF("ATGTGAA")
	'ATG'
	>>> rest_of_ORF("ATGAGATAGG")
	'ATGAGA'
	>>> rest_of_ORF("ATGTAGTAGTAG")					# Edge case
	'ATG'
	"""
	# TODO: implement this
	index = 0
	while index < len(dna):
		codon = dna[index:index + 3]														# groups 3 nucleotides to make codon, index by 3
		if codon in ["TAA", "TAG", "TGA"]:													# checks if stop codon is in string
			return dna[:index]																# if True, returns string up to stop codon
		index += 3																			# indexes by 3
	return dna

#print rest_of_ORF("ATGAGATAGG")

def find_all_ORFs_oneframe(dna):
	""" Finds all non-nested open reading frames in the given DNA
		sequence and returns them as a list.  This function should
		only find ORFs that are in the default frame of the sequence
		(i.e. they start on indices that are multiples of 3).
		By non-nested we mean that if an ORF occurs entirely within
		another ORF, it should not be included in the returned list of ORFs.

		dna: a DNA sequence
		returns: a list of non-nested ORFs
	>>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
	['ATGCATGAATGTAGA', 'ATGTGCCC']
	"""
	index = 0																				# unitiate index at 0
	segment = []																			# intitializes empty list
	while index < len(dna):																	# given that index is smaller than length of string
		codon = dna[index:index + 3]														# creates codons from nucleotides starting with first 3 letters
		if codon == "ATG":																	# if the codon is ATG
			segment.append(rest_of_ORF(dna[index:]))										# run rest_of_ORF starting at index
			index = index + len(rest_of_ORF(dna[index:]))									# new index starting at old + length of ORF output string
		else:
			index += 3																		# if ATG isn't found, move to next codon and check again

	return segment																			# returns list of ORFs

#print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    index = 0																				# initialize index
    all_ORFs = []																			# intitalize empty list
    while index < 3:																		# switches throught the 3 frames
    	all_ORFs = all_ORFs + find_all_ORFs_oneframe(dna[index:])	# adds all ORFs oneframe to the list
    	index += 1																			# indexes +1
    return all_ORFs 																		# returns list of all ORFs in all 3 frames



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    all_ORFs_both_strands = []																# initialize empty list
    all_ORFs_both_strands = all_ORFs_both_strands + find_all_ORFs(dna)						# adds all ORFs one strand to list
    reverse_complement_dna = get_reverse_complement(dna)									# finds 2nd strand
    all_ORFs_both_strands = all_ORFs_both_strands + find_all_ORFs(reverse_complement_dna)	# finds all ORFs second strand
    return all_ORFs_both_strands 															# returns all ORFs both strands


#print find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest_ORF = []																		# initialize empty list
    longest_ORF = longest_ORF + find_all_ORFs_both_strands(dna)								# adds all ORFs both strands to list

    for n in range(len(longest_ORF)):														# searches through the list
	    if longest_ORF[n] == max(longest_ORF, key=len):										# finds longest string
	    	return longest_ORF[n]															# returns longest string

#print longest_ORF("ATGCGAATGTAGCATCAAA")


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_random_ORF = []																	# initialize empty list
    for i in range(num_trials):																# repeats for num_trials
    	randomized_dna = shuffle_string(dna)												# finds strings of random dna
    	longest_random_ORF.append(longest_ORF(randomized_dna))								# adds random dna to list
	return max(longest_random_ORF, key=len)													# returns longest dna string

#print longest_ORF_noncoding("ATGCGAATGTAGCATCAAACGAATAGGTAGCATATGGCGTATGCTAGA", 10)       

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    s = ""																					# intitialize empty list																		
    for i in range(0, len(dna)-2, 3):														# for range of length of dna, indexes w/ step 3 (to isolate codons)
    		amino_acid = aa_table[dna[i:i+3]]												# translates each codon to an amino acid
    		s = s + amino_acid 																# adds amino acid to list
    return s 																				# returns list of amino acids

#print coding_strand_to_AA("ATGCCCGCTTTACAT")

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    viable_strings = []																		# intitialize empty list (for strings)
    viable_amino_acids = []																	# intitialize empty list (for amino acids)
    threshold = longest_ORF_noncoding(dna, 1500)											# sets threshold to longest random dna string
    real_dna = list(find_all_ORFs_both_strands(dna))										# sets real_dna equal to all the ORFs, both strands
    for i in range(len(real_dna)):															# searches through all the elements in list real_dna
    	if len(real_dna[i]) > len(threshold):												# compares real string to random string
    		viable_strings.append(real_dna[i])												# if real string is longer, adds it to list
    for i in range(len(viable_strings)):													# searches through all elements in viable_strings
    	a = coding_strand_to_AA(viable_strings[i])											# translates each string to amino acid sequence
    	viable_amino_acids.append(a)														# adds amino acids to list
    return viable_amino_acids																# returns list




dna = load_seq("./data/X73525.fa")

print gene_finder(dna)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
