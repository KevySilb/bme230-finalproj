import argparse

class commandLine():
    
    def __init__(self, inputs = None) -> None:
        self.parser = argparse.ArgumentParser(description = 'Command Line Usage')
        
        #required arguments 
        self.parser.add_argument('-k', '--kmerSize', type = int, default = 4, help = 'Specifies the kmer size')
        self.parser.add_argument('-i', '--inFile', type = str, default = None, help = 'fasta file')
        self.parser.add_argument('-o', '--outFile', type = str, default = None, help = 'name of .pkl file to be saved')

        if inputs is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inputs)