import collections

import matplotlib.pyplot as plt
from matplotlib.pyplot import show


def nucleo_counter_2(seq):
    """ counts the number of each nucleotide in a string and returns a dictionary"""

    return dict(collections.Counter(seq))


def readFASTQ(file):
    """reads a FASTQ file and returns a list of dna strings and a list  of strings
        of quolity encoding  symbols"""

    sequence = []
    quality = []

    with open(file, 'r') as file:
        while True:
            file.readline()
            seq = file.readline().rstrip()
            file.readline()
            quol = file.readline().rstrip()
            if len(seq) == 0:
                break
            else:
                sequence.append(seq)
                quality.append(quol)
    return sequence, quality


def to_sting(reads):
    """join FASTQ reads to  a single DNA string"""

    return ''.join(i for i in reads)


def to_reverse_complement(string):
    """create a reverse compelement string"""

    compelement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    reverse_complement = ''

    for i in string:
        reverse_complement = compelement[i] + reverse_complement

    return reverse_complement

def get_frequency_GC(reads):
    """ get the ratio of GC pairs to the total number of  pairs from a list of DNA strings"""
    
    GC = []
    AT = []
    for read in reads:
        for i  in read:
            if i == 'C' or i == 'G':
                GC.append(1)
            else:
                AT.append(1)
    total = sum(GC+AT)
    return float(sum(GC)/total)


def get_frequency_GC_str(string):
    """ get the ratio of GC pairs to the total number of  pairs from a  DNA strings """
    
    GC = []
    total = len(string)
    for i in string:
            if i == 'C' or i == 'G':
                GC.append(1)
    return float(sum(GC)/total)


def find_exacte_matching(pattern, string):
    """takes a DNA string and the sought pattern returns a list of the pattern's locations in the string"""

    location_indexes = []

    for i in range(len(string) - len(pattern)+1):
        match = True
        for n in range(len(pattern)):
            if string[i+n] != pattern[n]:
                match = False
                break
        if match:
            location_indexes.append(i)
    return location_indexes


def phred33ToQ(quolity):
    """convert character to quolity score  accroding to ASCII table"""
    
    return ord(quolity) - 33


def quality(quality_str):
    """ Create a list of quality scores"""

    quality = []
    for read in quality_str:
        for i in read:
            i = phred33ToQ(i)
            quality.append(i)
    return quality


def plot_quality(quality_list):
    """ Plot the histogram of quality scores"""

    plt.plot(range(len(quality_list)), quality_list)
    plt.xlabel('index')
    plt.ylabel('quqlity score')
    plt.show()
