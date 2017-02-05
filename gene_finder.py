# -*- coding: utf-8 -*-
"""
MiniProject 1: GeneFinder

@author: Kerry McConnaughay

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
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return "C"
    else:
        print('This is not possible', nucleotide)


# get_complement('A')
# get_complement('C')
# get_complement('H')


def get_reverse_complement(dna):

    """ Computes the reverse complementary sequence of DNA for the specified DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    complementary_dna = ''
    for nucleotide in dna:
        different_nucleotide = get_complement(nucleotide)
        # print(different_nucleotide)
        complementary_dna = complementary_dna + different_nucleotide
    # print(complementary_dna)

    new_dna = ''
    index = len(complementary_dna) - 1
    while index >= 0:
        reverse = complementary_dna[index]
        index = index - 1
        new_dna = new_dna + reverse
        # print(new_dna)
    return new_dna


# get_reverse_complement('ATGCCCGCTTT')
# get_reverse_complement('CCGCGTTCA')


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
    """
    r = len(dna) % 3
    new_length_dna = len(dna) - r
    index = 3  # where the index will start looking at 3 for frame stop codons
    while index < new_length_dna:
        if dna[index] == 'T':
            if dna[index + 1] == 'G':
                if dna[index + 2] == 'A':
                    return dna[:index]
            if dna[index + 1] == 'A':
                if dna[index + 2] == 'A' or 'G':
                    return dna[:index]
        index = index + 3
    return dna


# rest_of_ORF('ATGTGAA')
# rest_of_ORF("ATGAGATAGG")
    # i = 3  # where the index will start looking at 3 for frame stop codons
    # while i < len(dna):
        # if dna[i:i+3] == 'TGA' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TAG':
            # return dna[:i]
    # return dna


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
    orfs = []
    index = 0
    while index < len(dna):
        if dna[index:index + 3] == 'ATG':
            new_ORF = rest_of_ORF(dna[index:])
            orfs.append(new_ORF)
            index = index + len(new_ORF)
        else:
            index = index + 3
    return(orfs)


# find_all_ORFs_oneframe('ATGCATGAATGTAGATAGATGTGCCC')


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
    first = find_all_ORFs_oneframe(dna)
    second = find_all_ORFs_oneframe(dna[1:])
    third = find_all_ORFs_oneframe(dna[2:])
    all_orfs = first + second + third
    return all_orfs


# find_all_ORFs('ATGCATGAATGTAG')


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    first_strand = find_all_ORFs(dna)
    second_strand = find_all_ORFs(get_reverse_complement(dna))
    all_orfs_both_strands = first_strand + second_strand
    return all_orfs_both_strands


# find_all_ORFs_both_strands('ATGCGAATGTAGCATCAAA')


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """

    all_ORFS = find_all_ORFs_both_strands(dna)
    longest_ORFrame = max(all_ORFS, key=len)
    for orf in all_ORFS:
        if len(orf) == len(longest_ORFrame):
            return orf


# longest_ORF('ATGCGAATGTAGCATCAAA')


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specified DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    lengths_of_longest = []
    first_ORF = longest_ORF(dna)
    first_ORF_length = len(first_ORF)
    lengths_of_longest.append(first_ORF_length)
    for i in range(num_trials - 1):
        dna = shuffle_string(dna)
        longest_shuffled_ORF = longest_ORF(dna)
        length_shuffled_ORF = len(longest_shuffled_ORF)
        lengths_of_longest.append(length_shuffled_ORF)
    THE_LONGEST = max(lengths_of_longest)
    return THE_LONGEST


# no test units for this function


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
    amino_acids = ''
    r = len(dna) % 3
    for i in range(0, (len(dna) - r), 3):
        amino_acid = aa_table[dna[i: i + 3]]
        amino_acids += amino_acid
    return amino_acids


# coding_strand_to_AA("ATGCGA")


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    list_amino_acids = []
    threshold = longest_ORF_noncoding(dna, 1500)
    open_reading_frames = find_all_ORFs_both_strands(dna)
    # print(open_reading_frames)
    for i in open_reading_frames:
        if len(i) >= threshold:
            protein = coding_strand_to_AA(i)
            list_amino_acids.append(protein)
    return list_amino_acids


dna = load_seq("./data/X73525.fa")
print(gene_finder(dna))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
