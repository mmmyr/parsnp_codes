from Bio import SeqIO
from Bio import AlignIO

input_file = "parsnp.xmfa"

with open (input_file,'r') as file:
    contents = file.readlines()

modified_input = "parsnp_modified.xmfa"
with open (modified_input,'w') as file:
    for line in contents:
        if line.startswith('>'):
            line = '> ' + line[1:]
            file.write(line)
        else:
            file.write(line) 

alignment = list(AlignIO.parse(modified_input, "mauve"))

output_file = "parsnp.maf"
AlignIO.write(alignment, output_file, "maf") 