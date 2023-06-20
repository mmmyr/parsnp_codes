import re
from Bio.Seq import Seq
from Bio import SeqIO
my_seq = Seq("AGTACACTGGT")
rev_comp = my_seq.reverse_complement()
print(rev_comp)

files={}

log_file = xmfa_file =  'parsnp/P_2023_06_15_011137218904/parsnpAligner.log'
with open(log_file,'r') as file:
	contents = file.read()
	lines = contents.splitlines()
	for line in lines: 
		if len(line)!=0 and line[0]=="S":
			split_parts = line.split(':')
			path = split_parts[1].strip()
			#print(path)
			parts = path.split("/")
			real_path="/".join(parts[4:])
			match = re.search(r'\b(\d+)\b', line)
			number = match.group(1)
			files[number]=path
print(files)

xmfa_file =  'parsnp/P_2023_06_15_011137218904/parsnp.xmfa'

alignment_blocks = SeqIO.parse(xmfa_file, "fasta")

index = SeqIO.index(xmfa_file, "fasta")

# Access sequences by their IDs
for seq_id in index:
	sequence = index[seq_id].seq
	if sequence[-1] == "=":
		sequence = sequence[:-1]
	#print("sequence:",sequence)
	#print("seqid", seq_id)	
	idx=1
	while idx < len(seq_id):
		if seq_id[idx] == ':':
                	num_index=idx
		if seq_id[idx] == '-':
			coord1_index=idx
		idx+=1
	coord2_index=len(seq_id)
	seq_num = seq_id[0:num_index]
	coord1 = int(seq_id[num_index+1:coord1_index])
	coord2 = int(seq_id[coord1_index+1:coord2_index])
	print(seq_num)
	print(coord1)
	print(coord2)
	#if seq_num == '1':
	#	continue
	input_file = "parsnp/" + files[seq_num]
	for record in SeqIO.parse(input_file, "fasta"):
		seq = record.seq
       		# Extract the subsequence using slicing
		subsequence = seq[coord1:coord2+1]
       		# Print the subsequence
		#print("Subsequence:",subsequence)
		#print("Sequence:",sequence)
		if subsequence.upper() == sequence.upper():
			print("match")
		else: 
			print("unmatch")

