import os
from Bio import SeqIO
import argparse
# solve reference problem  - copy reference - alternative: make this as input - figure out the --flag / -flag thing 
# optimize the code
# convert it into mat format 

def compare(str1,str2):
    if len(str1) == len(str2):
        return False
    else: #
        return all(ch1 == ch2 for (ch1, ch2) in filter(lambda pair: '-' not in pair, zip(str1, str2)))

#parse CLI 
parser = argparse.ArgumentParser()
parser.add_argument('xmfa', type=str, help="the path to the xmfa file you are trying to verify")
parser.add_argument('ref', type=str, help="the path to the reference fna file")
parser.add_argument('fna', type=str, nargs='+', help="the path(s) to the fna file(s)")
args = parser.parse_args()

print(args.xmfa)
print(args.ref)
print(args.fna[0]) # this is a list 
directory_path = os.path.dirname(args.fna[0])
print(directory_path)

#parse headers 
files = {}
with open(args.xmfa,'r') as file:
    contents = file.readlines()

idx = 2 
line = contents[idx]
while line.startswith("##"):
    if line.startswith("##SequenceIndex"): 
        seq_num = int(line[16:])
        idx += 1
        line = contents [idx]
        file_path = line[15:]
        files[seq_num]=file_path.rstrip('\n')
    idx += 1 
    line = contents[idx]
print(files)

#parse xmfa output 
#improve this part afterwards 
alignment = SeqIO.parse(args.xmfa, "fasta")
for record in alignment:
   # print(record.id) #5:2605645-2608743
   # print(record.description) # + cluster1382 s46:p391742
   #print(record.seq)   sequence
    idx=0 
    while idx < len(record.id):
        if record.id[idx] == ':':
            num_index=idx
        if record.id[idx] == '-':
            coord1_index=idx
        idx+=1
    coord2_index=len(record.id)
    seq_num = int(record.id[0:num_index])
    coord1 = int(record.id[num_index+1:coord1_index])
    coord2 = int(record.id[coord1_index+1:coord2_index])
    length = coord2 - coord1 + 1

    idx=0
    while idx < len(record.description):
        if record.description[idx] == "s":
            idx1=idx
        if record.description[idx] == ":":
            idx2=idx
        idx+=1
    mul_num = int(record.description[idx1+1:idx2])
    start_coord = int(record.description[idx2+2:])
    print(seq_num, coord1, coord2,mul_num,start_coord,length)

    #read records from input file and make comparison 
    if seq_num == 1: 
        input_file = args.ref 
        continue
    else: 
        input_file = directory_path + '/' + files[seq_num]
    rec_num=0
    for mul in SeqIO.parse(input_file, "fasta"):
        mul_num-=1
        if mul_num==0:
            subsequence = mul.seq[start_coord:start_coord+length]
            if compare(subsequence, record.seq) == True:
                print("match")
            else: 
                print("unmatch")
