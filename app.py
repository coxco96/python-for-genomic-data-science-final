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

# get comparison operator depending on whether searching for shortest or longest value
def comparison_operator(sl):
    if sl == 'shortest':
        comparison = operator.lt
    else: # sl == 'longest'
        comparison = operator.gt
    return comparison

# find shortest or longest sequence in entire file
def shortest_longest(sl, dict):
    # initialize id and length to return
    id = ''
    length = 0
    # if looking for shortest length, re-initialize length to equal first sequence's length
    # since no seq will ever be less than 0
    if sl == 'shortest': 
        length = len(next(iter(dict.values())))
    # initialize list for tie values
    ties = []
    
    # choose greater or lesser than sign based on shortest or longest
    comparison = comparison_operator(sl)
    
    for k, seq in dict.items():

        if comparison(len(seq), length):
            # update length and id to shortest or longest length and id
            length = len(seq)
            id = k
            ties = []
            
        # if already equal...
        elif len(seq) == length:
            if k not in ties:
                ties.append(k)
            ties.append(id)
    # if there are ties, return that list and length tied for    
    if ties != []:
        return ties, length
    else: # otherwise, return the superlative id and its length
        return id, length
    
# handle what is printed for -L or -S args depending on which flag (longeset or shortest)    
def handle_LS(option):
    type = 'longest' if option == '-L' else 'shortest'
    id, length = shortest_longest(type, seqs)
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
    
    n = int(n)
   
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

# identify ORFs within a reading frame
def find_open_reading_frames(rf):
    
    # define start and stop codons
    start_codon = 'ATG'
    stop_codons = ('TAA','TAG','TGA')
    
    # initialize lists of locatoins for start codons
    start_codon_locations = []
    stop_codon_locations = []
    
    # find them and fill the locations lists
    for i, codon in enumerate(rf):
        if codon == start_codon: start_codon_locations.append(i)
        if codon in stop_codons: stop_codon_locations.append(i)
    
    # create list of ORF locations within this reading frame
    orf_locations = []
    for i in start_codon_locations:
        for j in stop_codon_locations:
            # start codon must come before stop codon and cannot overlap
            if i < j and j - i >= 3:
                orf_loc = [i, j]
                orf_locations.append(orf_loc)
                
   
                
    return orf_locations
    

# print(len(find_open_reading_frames(create_reading_frame(3, seqs['gi|142022655|gb|EQ086233.1|91 marine metagenome JCVI_SCAF_1096627390048 genomic scaffold, whole genome shotgun sequence']))))
# print(find_open_reading_frames(create_reading_frame(3, seqs['gi|142022655|gb|EQ086233.1|16 marine metagenome JCVI_SCAF_1096627390048 genomic scaffold, whole genome shotgun sequence'])))


# find shortest or longest ORF in entire file at specified reading frame n
def shortest_longest_orf(sl, n, id='', id_only=False):
    # get < or > operation depending on 'shortest' or 'longest' as sl
    comparison = comparison_operator(sl)
    # initialize dictionary of shortest or longest sequences for each id
    sl_dict = {}
    
    # find shortest or longest ORF of each individual sequence
    for id, seq in seqs.items():
        # initialize list for case of ties
        seq_ties = []
        rf = create_reading_frame(n, seq)
        orfs = find_open_reading_frames(rf)
        if orfs:
            if sl == 'longest': length = 0 # initialize with 0 if looking for longest
            else: length = orfs[0][1] - orfs[0][0] + 1 # initialize with first length if looking for shortest
            for i in orfs:
                orf_len = i[1] - i[0] + 1 + 2 # +1 to make inclusive and +2 for the 2 nucleotides that follow beginning of stop_codon_location
                if comparison(orf_len, length):
                    length = orf_len
                    seq_ties = [] # empty ties list if something else was longer/shorter
                if orf_len == length:
                    seq_ties.append(id)
            if len(seq_ties) <= 1:        
                # print(f"the {sl} length in {id[:50]} is {orf_len} characters")
                sl_dict[id] = {'len': orf_len, 'count': 1}
            else: 
                sl_dict[id] = {'len': orf_len, 'count': len(seq_ties)}
                # print(f"There was a {len(seq_ties)}-way tie for sequence id {id[:30]} at length {orf_len}")
    
    # if id was specified, return shortest or longest ORF for that id only, then break
    if id != '':
        if len(id) < 30:
            raise ValueError('Specified ID must be at least 30 characters long.')
        matching_keys = [key for key in sl_dict if key.startswith(id)]
        if not matching_keys:
            raise KeyError(f"No keys found starting with {id}.")
        elif len(matching_keys) > 1:
            raise KeyError(f"Multiple keys found starting with {id}. Try extending the identifier.")
        else: # identifer found in sl_dict
            if sl_dict[id]['count'] == 1:
                print(f"The {sl} ORF for specified key is {sl_dict[id]['len']} characters long.")
            elif sl_dict[id]['count'] > 1:
                print(f"The identifier you specified had a {sl_dict[id]['count']}-way tie for {sl} ORF, at {sl_dict[id]['len']} characters long.")
            
            # if -e or -t flags were used but neither -r nor -g were used, don't print overall longest/shortest for file
            if id_only: return 
                
    # find longest or shortest among ALL sequences, 
    # #using dictionary or shortest or longest of EACH sequence
    
    winner_ties = [] # initialize for case of tie
    winner_id = ''
    if sl == 'longest': winner_length = 0 
    else: # if shortest, initialize to first length that will be looped through
        winner_length = next(iter(sl_dict.values()))['len']

    # loop through longest or shorteset sequence ORFs
    for id, v in sl_dict.items():
        this_length = v['len']
        if comparison(this_length, winner_length):
            winner_length = this_length
            winner_ties = []
            winner_id = id
        # in instance of tie
        if this_length == winner_length:
            winner_ties.append(id)
            
    if len(winner_ties) <= 1:
        print(f"The {sl} ORF in the file is {winner_length} characters. Its sequence ID is {winner_id}.")
    else: # in the case of ties
        print(f"There was a {len(winner_ties)}-way tie for {sl} ORF at a length of {winner_length} characters.") 
        answer = input('Would you like a list of the ids (first 31 characters)? Enter YES or NO.')
        if answer == 'YES':
            for id, v in sl_dict.items():
                if id in winner_ties:
                    print(id[:31])
        elif answer != 'NO': print('invalid response'); 
        
def repeats(n):
    n = int(n)
    # initialize list of substrings extracted from all sequences of length n
    substr_list = []
    
    # loop through each sequence to add its substrings to substr_list
    for id, seq in seqs.items():
        substr_list += [seq[i:i+n] for i in range(0, len(seq)) if len(seq[i:i+n]) == n]

    # number of times the most frequently repeated string is repeated
    max_count = substr_list.count(max(substr_list, key=substr_list.count))
    # what that string is:
    most_repeated = [item for item in set(substr_list) if substr_list.count(item) == max_count]

    print(f"Among all sequences, the most frequently repeated string of {n} characters is {most_repeated}, occuring {max_count} times.")
    
    # ask if user wants to check if a specific string is repeated
    check_string = input("Want to check if a string repeats anywhere in the file? Enter it here. (If not, type 'n'.)")
    if check_string == 'n':
        return
    elif len(check_string) != n:
        raise ValueError(f"Your string is an invalid length of {len(check_string)}. You entered {n} as the substring length.")
    else:
        # get list of all substrings that repeat of length n
        # repeats_list = [substr for substr in set(substr_list) if substr_list.count(substr) > 1]
        count = substr_list.count(check_string)
        print(f"Your string {check_string} occurs {count} times in the file.")
     
     
        
        
        
        
    
    

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
       
        <filename>                file to process. must be in FASTA format. must be first arg.
        -h                        print this message
        -n                        returns number of sequences in file
        -i <identifier>           specify a sequence identifier
        -l                        get length of identifier. requires -i <identifier>
        -L                        get id and length of longest sequence in entire file
        -S                        get id and length of shortest sequence in entire file
        -o <readingFrameNumber>   prints number of open reading frames in specified sequence identifier for specified reading frame. requires arg 1, 2 or 3 for reading frame number to analyze. also requires -i <identifier>.
        -r <readingFrameNumber>   prints shortest ORF id and length among all sequences in file.
        -g <readingFrameNumber>   prints longest ORF id and length among all sequencse in file.
        -t <readingFrameNumber>   prints shortest ORF in specified identifer. requires -i <identifier>
        -e <readingFrameNumber>   prints longest ORF in specified identifier. requires -i <identifier>
        -p <repeatedSubstrLength> prints number of repeats of length <repeatedSubstrLength> in all sequences.
        """
    )
    
# create list of optional and required arguments
o, a = getopt.getopt(sys.argv[2:], 'hni:lLSo:r:g:t:e:p:')

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
    handle_LS('-L')
    
if '-S' in opts.keys():
    handle_LS('-S')
    
if ('-o' in opts.keys() or '-l' in opts.keys() or '-t' in opts.keys() or '-e' in opts.keys()) and '-i' not in opts.keys():
    raise ValueError('must include -i <SequenceIdentifier>')
    
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
        
        # return count of open reading frames on specified reading frame if -o flag is included
        if '-o' in opts.keys():
            rf = create_reading_frame(opts['-o'], seqs[id])
            orfs = find_open_reading_frames(rf)
            print(f"The specified sequence contains {len(orfs)} ORFs on reading frame no. {opts['-o']}")
        
        # return shortest ORF in identifier
        if '-t' in opts.keys():  
            if '-r' not in opts.keys() and '-g' not in opts.keys():
                id_only = True
            else: id_only = False
            
            shortest_longest_orf('shortest', opts['-t'], opts['-i'], id_only)
        
        # return longest ORF in identifier
        if '-e' in opts.keys():  
            
            if '-r' not in opts.keys() and '-g' not in opts.keys():
                id_only = True
            else: id_only = False
            
            shortest_longest_orf('longest', opts['-e'], opts['-i'], id_only)
           
           
            
if '-r' in opts.keys():
    shortest_longest_orf('shortest', opts['-r'])
    
if '-g' in opts.keys():
    shortest_longest_orf('longest', opts['-g'])
    
if '-p' in opts.keys():
    repeats(opts['-p'])
            
            
