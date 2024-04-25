file_name = input('Enter file name: ')

try:
    fasta_string = open(file_name).read()
except IOError:
    print(f"{file_name} doesn't exist")
    
print(fasta_string[:20])