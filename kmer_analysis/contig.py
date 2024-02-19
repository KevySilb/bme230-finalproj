# Usage python contig.py -i <inFile> -o <outFile> -k <kmer size integer>

import pandas as pd
import itertools as it
from fastaReader import FastAreader
from copy import deepcopy
import pickle
from commandLine import commandLine

class contig:
    """
    A class to calculate the probability of any kmer appearing in the sequence. 
    
    Attributes:
        head(str): FASTA header 
        seq(str): the contiguous sequence
        kmer_dict(dict): a copy of the count dictionary for counting kmers
        k(int): kmer size
    
    """
        
    def __init__(self, head, seq, templateDict, kmerMap) -> None:
        self.head = head
        self.seq = seq
        self.kmer_dict = deepcopy(templateDict)
        self.k = len(next(iter(kmerMap)))
        # count the occurance of each 4-mer in the sequence
        for i in range(len(self.seq) - self.k + 1):
            kmer = self.seq[i:i + self.k]
            if kmer not in kmerMap:
                continue
            key = tuple(sorted([kmer, kmerMap[kmer]]))
            if key in self.kmer_dict:
                self.kmer_dict[key] += 1.0
        # sum all the counts and normalize the counts into probabilities
        value_sum = float(0.0)
        for value in self.kmer_dict.values():
            if value > 0:
                value_sum += value
        if value_sum > 0:
            for key in self.kmer_dict.keys():
                self.kmer_dict[key] /= value_sum

def rev_compliment(kmer, nucMap):
    """
    Computes the reverse compliment of a kmer.

    Args:
        kmer (str): k-letter word from the alphabet {AGCT}
        nucMap (dictionary): unique mappings of nuclotides A->T, C->G

    Returns:
        str: reverse compliment
    """
    string = []
    for char in kmer:
        string.append(nucMap[char])
    return ''.join(string[::-1])
        
def make_temp_dict(k):
    """
    Constructs a template dictionary to be passed and copied to every contig class.

    Returns:
        dict: count dictionary of all possible kmers and their reverse compliments as tuples. {(kmer, revKmer) : 0}
        dict: a dictionary of kmers and their reverse compliments to be used for constructing tuples in the contig class
    """
    alphabet = ['A', 'T', 'G', 'C']
    nucMap = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    kmerPermutations = [''.join(p) for p in it.product(alphabet, repeat=k)]
    templateDict = dict()
    kmerMap = dict()
    # make kmerMap
    for kmer in kmerPermutations:
        revKmer = rev_compliment(kmer, nucMap)
        kmerMap[kmer] = revKmer
        key = tuple(sorted([kmer, revKmer]))
        if key not in templateDict:
            templateDict[key] = 0
    return templateDict, kmerMap
    
def main():
    cmd = commandLine()
    contigList = []
    templateKmerDict, kmerMap = make_temp_dict(cmd.args.kmerSize)
    myReader = FastAreader(cmd.args.inFile)
    for head, seq in myReader.readFasta(): 
        obj = contig(head, seq, templateKmerDict, kmerMap)
        contigList.append(obj)
    # build DataFrame
    matrix = []
    row_index = []
    col_index = [key for key in templateKmerDict.keys()]
    for obj in contigList:
        row_index.append(obj.head)
        # ensures column index lines up with correct values
        row_values = [obj.kmer_dict[key] for key in col_index]
        matrix.append(row_values)
    df = pd.DataFrame(matrix, index = row_index, columns = col_index)
    with open(cmd.args.outFile, 'wb') as file:
        pickle.dump(df, file, protocol=pickle.HIGHEST_PROTOCOL)
    
    print('Saved .pkl file as {}'.format(cmd.args.outFile))
    
if __name__ == "__main__":
    main()