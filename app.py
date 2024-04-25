#!usr/bin/python3
import sys
import getopt

# get and read file
file_name = sys.argv[1]
try:
    fasta_string = open(file_name).read()
except IOError:
    print(f"{file_name} doesn't exist")
    
# create dictionary of sequences
seqs = {}
lines = fasta_string.split('\n')
for i in range(len(lines)):
    line = lines[i]
    # if header...
    if line.startswith('>'):
        name = line[1:]
        seqs[name] = ''  
    # if sequence (not header) 
    else:
        seqs[name] = seqs[name] + line
 
def num_seqs():       
    return len(seqs)



"""
COMMAND LINE ARGS AND USAGE
"""

# define command line args
def usage():
    print(
        """
        app.py : accepts a FASTA file, builds a 
        dictionary with all sequences, and offers
        processing commands
       
        <filename>          file to process. must be in FASTA format. must be first arg.
        -h                  print this help message
        -l                  returns number of sequences in file
        -i <identifier>     specify a sequence identifier
        """
    )
    
# create list of optional (o) and required (a) arguments
o, a = getopt.getopt(sys.argv[2:], 'lhi:')

opts = {}
seqlen=0

# add key and val of optional args to to opts dictionary
for k,v in o:
    opts[k] = v
    
# if fasta file name not included, alert
if len(sys.argv) < 2:
    usage(); sys.exit("input FASTA file is missing")

# if -h, print usage message 
if '-h' in opts.keys():
    usage(); sys.exit()

    
if '-l' in opts.keys():
    # ensure a negative value was not specified for length
    print(f"Number of sequences in file: {num_seqs()}.")
    
if '-i' in opts.keys():
    # check if there is a seq id that starts with this. false by default.
    contains_id = False
    for id in seqs:
        if id.startswith(opts['-i']):
            contains_id = True
            seq_id = opts['-i']
            break
    if contains_id == False:
        print('sequence id not found. try wrapping id in quotes.')



    

# close file
# fasta_string.close()