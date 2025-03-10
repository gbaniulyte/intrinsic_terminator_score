import os
import sys
#path to tab-delimited table, Col1 = coordinate, col2 = strand, last column = sequence (DNA in the 5' to 3' direction)
term_table = input('Enter full path to a tab-delimited table where col1 = coordinate, col2 = strand (+ or -), last column = sequence: ')
if os.path.isfile(term_table) == False:
	print('Exited. Check your path (e.g. /home/user/Documents/3end_sequences.txt')
	sys.exit()
#path to output file with results
calculated_ts = input('Enter full path to the results file to be created (including file name, will be tab-delimited format): ')
#path to kinefold parameters file
path_to_kinefold_folder = input("Enter full path to kinefold folder that has files named 'example.req' etc. (home/user/.../kinefold_long_static/TEST): ")
if os.path.isdir(path_to_kinefold_folder) == False:
	print('Exited. Check your path (e.g. /home/user/Documents/kinefold_long_static/TEST')
	sys.exit()
headers = input('Input file has headers (y/n)?')

from collections import defaultdict
from decimal import Decimal
from collections import Counter
import random
import re
from math import e
import subprocess
from datetime import datetime
print('Started at:', datetime.now())
#write the header in the results
with open(calculated_ts, 'w') as ter_strenght:
	ter_strenght.write("3' end\tStrand\t3' end_seq\tKinefold_structure\tA-tract\tHairpin\tHp_structure\tLoop_seq\tU-Tract\tdGU\tdGL\tdGH\tdGHA\tdGA\tdGB\tTS\n")
#change paths and parameters in the default kinefold parameters file
with open(path_to_kinefold_folder + '/example.req', 'r') as p:
	read_req = p.readlines()
	read_req[1] = path_to_kinefold_folder + '/example.p\n'
	read_req[2] = path_to_kinefold_folder + '/example.e\n'
	read_req[3] = path_to_kinefold_folder + '/example.rnm\n'
	read_req[4] = path_to_kinefold_folder + '/example.rnms\n'
	read_req[5] = path_to_kinefold_folder + '/example.rnml\n'
	read_req[6] = path_to_kinefold_folder + '/example.rnm2\n'
	read_req[7] = path_to_kinefold_folder + '/example.dat\n'
	read_req[12] = '0\t\t# pseudoknots   1=yes 0=no \n'
	read_req[13] = '1\t\t# entanglements\t1=yes 0=no\n'
	read_req[14] = '2 20\t\t# simulation type: 1=renaturation; 2 20=cotrans. @ 20msec/nt \n'
	read_req[19] = '<SEQNAME>\n'
	read_req[20] = '<BASE>\n'
	read_req[21] = '<SEQUENCE>\n'
	read_req[22] = '<ZIPFILE>\n'
with open(path_to_kinefold_folder + '/example.req', 'w') as w:
	w.write(''.join(read_req))
#reverse complement sequence
def DNA_to_RNA(sequence, strand):
	sequence = sequence.lower()
	if strand == '+':
		sequence = sequence.replace('a', 'A') #don't need this if input is already reversed complement
		sequence = sequence.replace('c', 'C')
		sequence = sequence.replace('g', 'G')
		sequence = sequence.replace('t', 'U')
	elif strand == '-':
		sequence = sequence[::-1]
		sequence = sequence.replace('a', 'U')
		sequence = sequence.replace('c', 'G')
		sequence = sequence.replace('g', 'C')
		sequence = sequence.replace('t', 'A')
	elif strand == 'RNA':
		sequence = sequence.upper()
	return sequence
#write kinefold data file
def kinefold_dat_write(sequence):
	with open(path_to_kinefold_folder + '/example.dat', 'w') as f:
		f.write('< example\n' + sequence)
#rewrite req file with random seed and specific molecular time
def change_req(time, text):
	with open(path_to_kinefold_folder + '/example.req', 'w') as k:
		text[0] = str(random.randrange(1000,100000)) + '\t\t# random seed\n'
		text[11] = str(time) + '\t\t# folding time requested in msec\n'
		k.write(''.join(text))
#pass sequence to kinefold and return kinefold output
def run_kinefold(): #probably don't need input because those won't change
	kinefold_exe = path_to_kinefold_folder + '/kinefold_long_static'
	kinefold_req_file = path_to_kinefold_folder + '/example.req'
	subprocess.run([kinefold_exe, kinefold_req_file], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
#detect molecular time required for kinefold
def get_mol_time():
	with open(path_to_kinefold_folder + '/example.i', 'r') as f: #check if enough time for folding (can be different results)
		text1 = f.read().splitlines()
		try:
			mol_time = int(text1[8][find_last(text1[8], ':')+2:find_last(text1[8], ')')])
			return mol_time
		except IndexError:
			return 'kinefold_fail'
#find last character in a string
def find_last(sequence, character):
	position = len(sequence) - sequence[::-1].find(character) -1
	return position
#calculate delta GU based on Sugimoto 1996 paper
def dGU_calc(sequence):
	term_par = {'AA': -1.0, 'AC': -2.1, 'AG': -1.8, 'AU': -0.9, 'CA': -0.9, 'CC': -2.1, 'CG': -1.7, 'CU': -0.9, 'GA': -1.3, 'GC': -2.7, 'GG': -2.9, 'GU': -1.1, 'UA': -0.6, 'UC': -1.5, 'UG': -1.6, 'UU': -0.2}
	sum_term = 0
	for i in range(0, 7):
		dinucleotide = sequence[i:i+2]
		sum_term += term_par[dinucleotide] #can merge into one line
	G_value = 3.1 + sum_term  # can merge into one line
	return G_value
#Get raw dot-bracket annotation from kinefold structure
def get_structure_from_kinefold():
	with open(path_to_kinefold_folder + '/example.p', 'r') as f:
		text = f.read().splitlines() #turn all lines into a list
		for ind, string in enumerate(text): #need position and string of each line
			if str(lenght) + ' over ' + str(lenght) + ' bases' in string:  #need one specific line from the whole text
				seq1 = text[ind]
				kine_seq = seq1[0:seq1.index('|')].strip() #get spaced out sequence with kine brackets (so it matched helix annotation in terms of list position)
				helix1 = text[ind + 1]
				helix = helix1[0:helix1.index('H')-1].strip() #get spaced out kine helix annotation for sequence
	try: #if kinefold file never picked up variables, return fail
		kine_seq
	except UnboundLocalError:
		return 'Fail'
	pos = 0
	dot_bracket = ''
	while pos < len(kine_seq): #read brackets sequencially and identify which ones are in stem or loop or bulge
		if kine_seq[pos] == 'A' or kine_seq[pos] == 'U' or kine_seq[pos] == 'G' or kine_seq[pos] == 'C':
			dot_bracket += '.'
			pos += 1
		elif kine_seq[pos] == ' ':
			pos += 1
		elif kine_seq[pos] == '[':  #start of the stem
			pos += 1
			stem_side = helix[pos + 1:helix.find('-', pos + 1)].count("'")  # if 0 then it's ( if 1 then it's )
			while kine_seq[pos] != '[' and kine_seq[pos] != ']':  #loop until reach end of teh current stem
				if kine_seq[pos] == 'A' or kine_seq[pos] == 'U' or kine_seq[pos] == 'G' or kine_seq[pos] == 'C':
					if stem_side == 0:
						dot_bracket += '('
					elif stem_side == 1:
						dot_bracket += ')'
					pos += 1
				elif kine_seq[pos] == ' ':
					pos += 1
				elif kine_seq[pos] == '^':
					pos += 1
					stem_side = helix[pos + 1:helix.find('-', pos + 1)].count("'")  #needed if change in direction in the middle of continous stem
		elif kine_seq[pos] == ']':
			pos += 1
	return dot_bracket
#find U-tract
def parse_kinefold_structure(term_region, kine_dot_bracket):
	loop_list = re.finditer('\(\.+\)', kine_dot_bracket)  #calleable iterator, use next(l) to iterate through each item
	#get loop start and end position for all loops
	loop_dict = defaultdict(dict)
	for value in loop_list:
		start = value.start() + 1
		end = value.end() - 2
		loop_dict.update({'loop' + str(start): [start, end]})

	#check if loop stem is at least 4nt long, if no, delete that entry
	for loop in list(loop_dict.keys()):
		boundary = loop_dict[loop]
		right_stem_end = kine_dot_bracket.find('(', boundary[1] + 1)
		if right_stem_end >= 0:
			bracket_count = kine_dot_bracket[boundary[1] + 1:right_stem_end].count(')') #use string until the next loop
		else:
			bracket_count = kine_dot_bracket[boundary[1] + 1:].count(')') #use string until the end
		if bracket_count < 4:
			del loop_dict[loop]

	#find all Utracts and keep the one with the highest GU value
	GU_dict = defaultdict(dict)
	for loop, boundary in loop_dict.items(): #check for dGU in every loop
		GU_dict_temp = defaultdict(dict)
		loop_pos = boundary[1]
		closest_next_loop = kine_dot_bracket.find('(', loop_pos + 1)
		if closest_next_loop < 0: #if it's the last loop
			stem_end = kine_dot_bracket.rfind(')', loop_pos + 1)
		else:
			stem_end = kine_dot_bracket.rfind(')', loop_pos + 1, closest_next_loop) #find last 3' stem base between loop and next stem
		pos = loop_pos + 6
		if stem_end < 0:
			print('stem issue at line 127')
		while pos <= stem_end + 1 and pos < len(kine_dot_bracket): #only look  from 6nt in the stem until the end of the predicted stem
			if term_region[pos] == 'U' and kine_dot_bracket[pos-1] != '.': #kine_dot_bracket[pos-3:pos] == ')))':
				u_tract = term_region[pos:pos + 8] #take 8nt each time
				if len(u_tract) >= 8: #if it's shorter, can't calculate GU correctly
					GU_dict_temp.update({pos: [dGU_calc(u_tract), loop]}) #save u start, GU and which loop it belongs to
				pos += 1
			else:
				pos += 1 #if no U got o next position
		if bool(GU_dict_temp) == False: #check if any Utrack starting with U is found (not emtry dictionary)
			bulge_list = re.finditer('\)\.+\)', kine_dot_bracket)#when no Utrack tract with U...
			bulge_dict = defaultdict(dict)
			for val in bulge_list: #find all regions after a set of brackets...
				start = val.start() + 1
				end = val.end() - 2
				if start >= loop_pos and end <= stem_end: #...inbetween loop and final stem end #EDIT: was loop_pos + 6
					bulge_dict.update({'bulge' + str(start): [start, end]})
			for bulge, boundary in bulge_dict.items():
				u_tract = term_region[boundary[0]:boundary[0] + 8]  #if not, utrack is 8bp immediately after the stem
				if len(u_tract) >= 8: #if it's shorter, can't calculate GU correctly
					GU_dict_temp.update({boundary[0]: [dGU_calc(u_tract), loop]})
		if bool(GU_dict_temp) == False:
			u_tract = term_region[stem_end + 1:stem_end + 9]
			if len(u_tract) >= 8: #if it's shorter, can't calculate GU correctly
				GU_dict_temp.update({stem_end + 1: [dGU_calc(u_tract), loop]})
		GU_dict.update(GU_dict_temp)

	if bool(GU_dict) == False:
		return 'Fail'
	max_GU_pos = max(GU_dict, key=lambda y: GU_dict[y][0])
	U_track_12 = term_region[max_GU_pos:max_GU_pos + 12]
	dGU_val = GU_dict[max_GU_pos][0]
	temp_start = loop_dict[GU_dict[max_GU_pos][1]][0] #find the start and end positions of the loop
	temp_end = loop_dict[GU_dict[max_GU_pos][1]][1] #end of the loop
	loop = term_region[temp_start:temp_end+1]
	right_bracket_count = kine_dot_bracket.count(')', temp_end + 1, max_GU_pos) #count right stem base pair no
	left_pos = temp_start - right_bracket_count #ideally would be the same number
	left_bracket_count = kine_dot_bracket[left_pos:temp_start].count('(')
	if kine_dot_bracket.count(')', left_pos, temp_start) > 0:
		while kine_dot_bracket.count(')', left_pos, temp_start) > 0: #recalculate where stem ends based on right stem
			left_pos +=1
		left_bracket_count = kine_dot_bracket[left_pos:temp_start].count('(')
		right_bracket_count = kine_dot_bracket.count(')', temp_end+1, temp_end+1+left_bracket_count)
		move_right = 0
		while True:
			if left_bracket_count == right_bracket_count:
				hairpin_str = kine_dot_bracket[left_pos:temp_end+1+left_bracket_count + move_right]
				hairpin_seq = term_region[left_pos:temp_end+1+left_bracket_count + move_right]
				u_tract = term_region[temp_end+1+left_bracket_count + move_right:temp_end+1+left_bracket_count+ move_right+8]
				U_track_12 = term_region[temp_end+1+left_bracket_count + move_right:temp_end+1+left_bracket_count+ move_right+12]
				dGU_val = dGU_calc(u_tract)
				break
			elif left_bracket_count != right_bracket_count:
				move_right += 1
				right_bracket_count = kine_dot_bracket.count(')', temp_end+1, temp_end+1+left_bracket_count + move_right)
	elif kine_dot_bracket.count(')', left_pos, temp_start) == 0:
		while True: #use this when don't need to change left side brackets
			if left_bracket_count == right_bracket_count:
				hairpin_str = kine_dot_bracket[left_pos:max_GU_pos]
				hairpin_seq = term_region[left_pos:max_GU_pos]
				break
			elif left_bracket_count != right_bracket_count:
				left_pos -= 1
				left_bracket_count = kine_dot_bracket.count('(', left_pos, temp_start)
	if left_pos >= 8: #make sure it's not too close to 5' end
		A_track = term_region[left_pos-8:left_pos]
	elif left_pos < 8:
		A_track = term_region[0:left_pos]
	return A_track, hairpin_seq, hairpin_str, loop, U_track_12, dGU_val
#put structures through RNAeval
def run_RNAeval(hp_sequence, structure):
	RNAeval_input = path_to_kinefold_folder + '/RNAeval_input.txt'
	RNAeval_out = path_to_kinefold_folder + '/RNAeval_output.txt'
	with open(RNAeval_input, 'w') as f:
		f.write(hp_sequence + '\n' + structure)
	result = subprocess.run(['RNAeval', RNAeval_input, '-d0', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)#.stdout.decode('utf-8')
	#error = result.stderr.decode('utf-8')
	with open(RNAeval_out, 'w') as o:
		o.write(result.stdout.decode('utf-8'))
#get values from RNAeval
def get_dG_from_RNAeval():
	RNAeval_out = path_to_kinefold_folder + '/RNAeval_output.txt'
	with open(RNAeval_out, 'r') as f:
		RNAeval_lines = f.read().splitlines()
	base_count = 1
	dGB_val = 0.0
	for line in RNAeval_lines:
		if line.startswith('Interior') and base_count <= 3:
			G_temp = float(line[line.index(':') + 1:].lstrip())/100
			if base_count <= 3:
				dGB_val += G_temp
				base_count += 1
		elif line.startswith('Hairpin'):
			dGL_val = float(line[line.index(':') + 1:].lstrip()) / 100
		elif line.startswith('.') == True or line.startswith('(') == True:
			dGH_val = float(line[find_last(line, '(') + 1:find_last(line, ')')])
	try:
		dGL_val
	except UnboundLocalError:
		return 'Fail'
	return dGB_val, dGL_val, dGH_val
#get dGA from RNAfold
def get_dGHA_from_RNAfold(sequence, structure):
	RNAfold_input = path_to_kinefold_folder + '/RNAfold_input.txt'
	with open(RNAfold_input, 'w') as f:
		f.write(sequence + '\n' + structure)
	result = subprocess.run(['RNAfold', RNAfold_input, '-d0', '-C'], stdout=subprocess.PIPE).stdout.decode('utf-8')
	#dGHA_val = float(result.split()[2][1:-1])
	dGHA_val = float(result[find_last(result, '(')+1:find_last(result, ')')])
	return dGHA_val
#loop through terminator sequences and process each individually
with open(term_table, 'r') as a, open(calculated_ts, 'a') as ter_strenght:
	if headers == 'y': #skip header if present
		next(a)
	for line in a:
		splitlines = line.split('\t')
		name = splitlines[1]+splitlines[0] #format has to be: col1,2 - some ID(e.g. coord and strand)
		term_region = splitlines[-1].rstrip() #last column has to be sequence
		#transform to RNA
		term_region = DNA_to_RNA(term_region, '+') #check if some sequences are reversed (-), or already 'RNA'
		lenght = len(term_region)
		#pass termn region through Kinefold
		kinefold_dat_write(term_region)
		run_kinefold() #default seed, could be modified in .req file
		#get molecular time
		molecular_time = get_mol_time()
		if molecular_time == 'kinefold_fail':
			ter_strenght.write(name[1:] + '\t' + name[0] + '\t' +term_region + '\t' + '\tNA\tNA\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n')
			continue
		#process kinefold output n times
		structure_list = []
		for i in range(0, 20):
			#change parameters for kinefold
			change_req(molecular_time, read_req)
			run_kinefold()
			kine_dot_bracket = get_structure_from_kinefold()
			structure_list.append(kine_dot_bracket)
		freq_structure = Counter(structure_list)
		kine_dot_bracket = max(freq_structure, key=freq_structure.get)
		#get  values based on kinefold
		parsed = parse_kinefold_structure(term_region, kine_dot_bracket)
		if parsed == 'Fail':
			ter_strenght.write(name[1:] + '\t' + name[0] + '\t' +term_region + '\t' + kine_dot_bracket + '\tNo_8nt_U_tract\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n')
			continue
		Atrack = parsed[0]
		hairpin_seq = parsed[1]
		hairpin_str = parsed[2]
		loop_seq = parsed[3] #migth include first base of the stem
		Utrack12 = parsed[4]
		dGU = parsed[5]
		#pass hairpin through RNAeval
		run_RNAeval(hairpin_seq, hairpin_str)
		parsed1 = get_dG_from_RNAeval()
		if parsed1 == 'Fail':
			ter_strenght.write(name[1:] + '\t' + name[0] + '\t' +term_region + '\t' + kine_dot_bracket +'\t' + Atrack + '\t' + hairpin_seq + '\t' + hairpin_str + '\t' + loop_seq + '\t' + Utrack12 + '\t' + str(dGU) + '\tRNAeval_fail\n')
			continue
		dGB = parsed1[0]
		dGL = parsed1[1]
		dGH = parsed1[2]
		#Run RNAfold and get dGHA value
		seq_for_RNAfold = Atrack + hairpin_seq + Utrack12[0:8]
		stru_for_RNAfold = len(Atrack) * '.' + hairpin_str + len(Utrack12[0:8]) * '.'
		dGHA = get_dGHA_from_RNAfold(seq_for_RNAfold, stru_for_RNAfold)
		dGA = dGHA - dGH
		#calculate Ts based on Chen et al 2013 paper
		try:
			Ts = 1 + (1 / (0.005 * e ** (0.6 * dGL) + (6.0 * e ** (0.45 * (dGB + dGA - dGU))) * (1 + 0.005 * e ** (0.6 * dGL))))
		except OverflowError:
			Ts = 'Impossible'
		ter_strenght.write(name[1:] + '\t' + name[0] + '\t' + term_region + '\t' + kine_dot_bracket + '\t' + Atrack + '\t' + hairpin_seq + '\t' + hairpin_str + '\t' + loop_seq + '\t' + Utrack12 + '\t' + str(dGU) + '\t' + str(dGL) + '\t' + str(dGH) + '\t' + str(dGHA) + '\t' + str(dGA) + '\t' + str(dGB) + '\t' + str(Ts) + '\n')
os.remove(path_to_kinefold_folder + '/RNAeval_output.txt')
os.remove(path_to_kinefold_folder + '/RNAeval_input.txt')
os.remove(path_to_kinefold_folder + '/RNAfold_input.txt')
print('Finished!', datetime.now())