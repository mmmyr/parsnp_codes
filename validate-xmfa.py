import os
from Bio import SeqIO
import argparse
import re

#python validate-xmfa.py parsnp.xmfa ref/GCF_000219945.1_ASM21994v2_genomic.fna test/*.fa

def compare(str1,str2):
    if len(str1) == len(str2):
        return False
    else: 
        return all(ch1 == ch2 for (ch1, ch2) in filter(lambda pair: '-' not in pair, zip(str1, str2)))

#parse CLI 
parser = argparse.ArgumentParser()
parser.add_argument('xmfa', type=str, help="the path to the xmfa file you are trying to verify")
parser.add_argument('ref', type=str, help="the path to the reference fna file")
parser.add_argument('fna', type=str, nargs='+', help="the path(s) to the fna file(s)")
args = parser.parse_args()

print(args.xmfa)
print(args.ref)
print(args.fna[0]) 
directory_path = os.path.dirname(args.fna[0])
print(directory_path)

#read the file 
with open(args.xmfa,'r') as file:
    contents = file.readlines()

#parse headers 
files = {}
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
alignment = SeqIO.parse(args.xmfa, "fasta")
for record in alignment:
    #delete '=' from the end of sequence 
    if record.seq.endswith('='):
        record.seq = record.seq[:-1] 

    numbers = re.split(r'\D+', record.description)
    seq_num = int(numbers[0])
    coord1 = int(numbers[1])
    coord2 = int (numbers[2])
    mul_num = int(numbers[4])
    start_coord = int(numbers[5])
    length = coord2 - coord1 + 1
    #print(seq_num, coord1, coord2,mul_num,start_coord,length)

    #read records from input file and make comparison 
    if seq_num == 1: 
        input_file = args.ref 
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

