#!usr/bin/python3
import sys
import getopt
import operator

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
 
# get number of sequences in file
def num_seqs():       
    return len(seqs)

# get a specific id's sequence length
def seq_len(id):
    return len(seqs[id])

# find shortest or longest sequence in entire file
def shortest_longest(sl):
    # initialize id and length to return
    id = ''
    length = 0
    # initialize list for tie values
    ties = []
    
    # choose greater or lesser than sign based on shortest or longest
    if sl == 'shortest':
        comparison = operator.lt
    else: # sl == 'longest'
        comparison = operator.gt
    
    for k, seq in seqs.items():
        # if looking for shortest length, re-initialize length to equal first sequence's length
        # since no seq will ever be less than 0
        if sl == 'shortest': 
            length = len(seq)
        if comparison(len(seq), length):
            # update length and id to shortest or longest length and id
            length = len(seq)
            id = k
            ties = []
            
        # if already equal...
        if len(seq) == length:
            if k not in ties:
                ties.append(k)
            ties.append(id)

    # if there are ties, return that list and length    
    if ties != []:
        return ties, length
    else: # otherwise, return the superlative id and its length
        return id, length
    
# handle what is printed for -L or -S args depending on which flag (longeset or shortest)    
def handle_option(option):
    type = 'longest' if option == '-L' else 'shortest'
    id, length = shortest_longest(type)
    if isinstance(id, str):
        print(f"The {type} sequence is {length} characters. Its id begins with {id[:30]}.")
    if isinstance(id, list):
        ids_string = '\n    '.join(i[:30] for i in id) # here because f-strings can't use backslash
        print(f"""
    There is a {len(id)}-way tie for {type} sequence, at length {length}.
    The beginnings of each id of {type} seq length are:
    
    {ids_string}
              """
            )

# create reading frame 1, 2 or 3 for a sequence
def create_reading_frame(n, seq):
   
    # ensure n is 1, 2 or 3 for valid reading frame
    if n not in [1, 2, 3]:
        raise ValueError("n for reading frame must == 1, 2 or 3")
    
    # initialize list for reading frame
    rf = []
    
    # change n from position to index
    n -= 1
    
    # fill reading frame with 3-nucleotide codons, starting at n position of sequence
    rf = [seq[i+n:n+i+3] for i in range(0, len(seq), 3)]
    
    # if last codon is an empty string, delete it
    if rf[len(rf)-1] == '':
        del rf[len(rf)-1]
    
    return rf
    
    

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
        -h                  print this message
        -n                  returns number of sequences in file
        -i <identifier>     specify a sequence identifier
        -l                  get length of identifier. requires -i <identifier>
        -L                  get id and length of longest sequence in entire file
        -S                  get id and length of shortest sequence in entire file
        
        """
    )
    
# create list of optional (o) and required (a) arguments
o, a = getopt.getopt(sys.argv[2:], 'hni:lLS')

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

if '-n' in opts.keys():
    # ensure a negative value was not specified for length
    print(f"Number of sequences in file: {num_seqs()}.")
    
if '-L' in opts.keys():
    handle_option('-L')
    
if '-S' in opts.keys():
    handle_option('-S')
    
if '-i' in opts.keys():
    # check if there is a seq id that starts with this. false by default.
    contains_id = False
    full_id = ''
    for id in seqs:
        if id.startswith(opts['-i']):
            id = id
            contains_id = True
            seq_id = opts['-i']
            break
    if contains_id == False:
        print('sequence id not found. try wrapping id in quotes.')
    else: 
        # return length of identifier sequence if -l flag included
        if '-l' in opts.keys():
            print(f"Length of this identifier is: {seq_len(id)}.")
            
