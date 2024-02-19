import sys
import io

class FastAreader:
    """
    Class for reading a FASTA file
    
    This class provides methods to open and read a FASTA file, 
    either from a specified path or from standard input.
    """
    
    def __init__ (self, fname=None):
        """
        Constructor for FastAreader class.
        
        Args:
            fname: Name of file to read. If None, reads from stdin.
        """
        self.fname = fname
            
    def doOpen (self):
        """
        Open the file specified during object creation, or stdin if no file was specified.
        
        Returns:
            File handle to the opened file or stdin.
        """
        if self.fname == None:
            return io.TextIOWrapper(sys.stdin.buffer, encoding='utf-8')
        else:
            return open(self.fname, 'r',  encoding= 'utf-8')
        
    def readFasta (self):
        """
        Read FASTA format sequences from the opened file or stdin.
        
        Yields:
            tuple: Contains the header and sequence of each FASTA record.
        """
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence