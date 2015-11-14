###############################################################################
###############################################################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
###############################################################################
###############################################################################


def error_check(line, length):
	"""
	   Description:	Error check to make sure Annotation Table is being filled out
	   Input:		Annotation Table Line, Line Length 
	"""
	errors = {
		 	   3:'Row length error in init_table.py - add_fasta() TranscriptID: ',
			   5:'Row length error in init_table.py - add_swissprot_top_blastx() TranscriptID: ',
			   7:'Row length error in init_table.py - add_protid() TranscriptID: ',
			   8:'Row length error in init_table.py - add_swissprot_top_blastp() TranscriptID: ',
			   9:'Row length error in init_table.py - add_pfam() TranscriptID: ',
			  10:'Row length error in init_table.py - add_signalp() TranscriptID: ',
			  11:'Row length error in init_table.py - add_tmhmm() TranscriptID: ',
			  12:'Row length error in init_table.py - add_rnammer() TranscriptID: ',
			  13:'Row length error in init_table.py - add_trembl_top_blastx() TranscriptID: ',
			  14:'Row length error in init_table.py - add_trembl_top_blastp() TranscriptID: ' }

	if len(line) != length:
		error_message = 'Annotation Table line incorrect size\n' + errors[length] + line[0] + '\nLineLength: ' + str(len(line)) + ' ' + '\t'.join(line)
		raise ValueError(error_message)
	return


def add_space(annot_table):
	"""
	   Description:	Space generator for empty files
	   Input:		Annotation Table
	   Return:		Annotation Table
	"""
	for index in range(len(annot_table)):
		annot_table[index].append(".")
	return annot_table


def compare(transcriptid_file, transcriptid_table):
	"""
	   Description:	Compares if transcriptid from a file is > transcriptid in the row of the annot_table
	   Input:		TranscriptID, TranscriptID
	   Return:		Boolean
	"""
	transcriptid_file = [int(x[1:]) for x in (transcriptid_file.split('_'))]
	transcriptid_table = [int(x[1:]) for x in (transcriptid_table.split('_'))]
	for x, y in zip(transcriptid_file, transcriptid_table):
		if x > y:
			return True
	return False


def conversion_transcriptid_geneid(conversion_file):
	"""
	   Description:	Creates dictionary for TranscriptID -> GeneID conversion
	   Input:		GeneTransMap File
	   Return:		Dictionary 
	"""
	convDt = {}
	file = open(conversion_file)
	for line in file:
		line = line.strip('\n').split()
		first = line[0]
		second = line[1]
		convDt[second] = first
	file.close()
	return convDt


def add_fasta(fasta_file, transmap_file, nr_file):
	"""
	   Description:	Initiates Table with Columns 0,1,2 - TranscriptID, GeneID, TranscriptLength
	   Input:		Fasta File, GeneTransMap File
	   Return:		Annotation Table
	"""
	annot_table = []
	file = open(fasta_file)
	transcriptid = ''
	length = 0
	geneid = conversion_transcriptid_geneid(transmap_file)
	lengths = {}
	for line in file:
		if '>' == line[0]:
			if length > 0:
				annot_table += [[transcriptid, geneid[transcriptid], str(length)]]
				lengths[transcriptid] = length
			length = 0
			transcriptid = line.split(' ')[0][1:]
		else:
			length += len(line.strip())
	annot_table += [[transcriptid, geneid[transcriptid], str(length)]]
	file.close()
	if nr_file == 'NONE':
		return annot_table, 0
	return annot_table, lengths


def add_swissprot_top_blastx(annot_table, swissprotx_file):
	"""
	   Description:	Adds Columns 3,4 - SwissProtBlastX, SwissProtBlastXLength
	   Input:		Annotation Table, SwissProtBlastX File
	   Return:		Annotation Table
	"""
	file = open(swissprotx_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		line = line.strip('\n').split('\t')
		if prev == line[0]:
			continue
		prev = line[0]
		while (True):
			error_check(annot_table[index], 3)
			if line[0] == annot_table[index][0]:
				annot_table[index].append(line[12])
				annot_table[index].append(line[14])
				index += 1
				break
			if compare(line[0], annot_table[index][0]):
				annot_table[index].append('.')
				annot_table[index].append('.')
				index += 1
			else:
				break
	while (index < length):
		error_check(annot_table[index], 3)
		annot_table[index].append('.')
		annot_table[index].append('.')
		index += 1
	file.close()
	return annot_table


def add_protid(annot_table, trans_file):
	"""
	   Description:	Adds Columns 5,6 - ProteinID, ProteinCoordinates
	   Input:		Annotation Table, Transdecoder File
	   Return:		Annotation Table
	"""
	file = open(trans_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		if line[0] != '>':
			continue
		line = line.strip('\n').split()
		transcriptid = line[8].split(':')[0]
		if prev == transcriptid:
			continue
		prev = transcriptid
		while (index < length):
			error_check(annot_table[index], 5)
			if transcriptid == annot_table[index][0]:
				proteinid = line[0][1:]
				protein_coordinate = line[8].split(':')[1].replace('(', '[').replace(')', ']')
				annot_table[index].append(proteinid)
				annot_table[index].append(protein_coordinate)
				index += 1
				break
			if compare(transcriptid, annot_table[index][0]):
				annot_table[index].append('.')
				annot_table[index].append('.')
				index += 1
			else:
				break
	while (index < length):
		error_check(annot_table[index], 5)
		annot_table[index].append('.')
		annot_table[index].append('.')
		index += 1
	file.close()
	return annot_table


def add_swissprot_top_blastp(annot_table, swissprotp_file):
	"""
	   Description:	Adds Column 7 - SwissProtBlastP
	   Input:		Annotation Table, SwissProtBlastP File
	   Return:		Annotation Table
	"""
	file = open(swissprotp_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		line = line.strip('\n').split('\t')
		transcriptid = line[0].split('|')[0]
		if prev == transcriptid:
			continue
		prev = transcriptid
		while (True):
			error_check(annot_table[index], 7)
			if transcriptid == annot_table[index][0]:
				annot_table[index].append(line[1]) ###could change to 12###
				index += 1
				break
			if compare(transcriptid, annot_table[index][0]):
				annot_table[index].append('.')
				index += 1
			else:
				break
	while (index < length):
		error_check(annot_table[index], 7)
		annot_table[index].append('.')
		index += 1
	file.close()
	return annot_table


def add_pfam(annot_table, pfam_file):
	"""
	   Description:	Adds Column 8 - PFAM
	   Input:		Annotation Table, PFAM File
	   Return:		Annotation Table
	"""
	file = open(pfam_file)
	index = 0
	length = len(annot_table)
	first = True
	for line in file:
		if index == length:
			break
		if '#' == line[0]:
			continue
		line = line.strip('\n').split()
		transcriptid = line[3].split('|')[0]
		while (True):
			if first:
				error_check(annot_table[index], 8)
			if transcriptid == annot_table[index][0]:
				if first:
					#temp = line[1] + '^' + line[0] + '^' + ' '.join(line[22:]) + '^' + line[15] + '-' + line[16] + '^E:' + line[10]
					temp = line[1] + '^' + ' '.join(line[22:])
					annot_table[index].append(temp)
					first = False
					break
				else:
					#temp = line[1] + '^' + line[0] + '^' + ' '.join(line[22:]) + '^' + line[15] + '-' + line[16] + '^E:' + line[10]
					if line[1] not in annot_table[index][-1]:
						temp = line[1] + '^' + ' '.join(line[22:])
						annot_table[index][-1] = annot_table[index][-1] + '`' + temp
					break
			if compare(transcriptid, annot_table[index][0]):
				#prevents extra '.' from being added
				if first:
					annot_table[index].append('.')
				first = True
				index += 1
			else:
				break
	if len(annot_table[index]) > 8:
		index += 1
	while (index < length):
		error_check(annot_table[index], 8)
		annot_table[index].append('.')
		index += 1
	file.close()
	return annot_table


def add_signalp(annot_table, signalp_file):
	"""
	   Description:	Adds Column 9 - SignalP
	   Input:		Annotation Table, SignalP File
	   Return:		Annotation Table
	"""
	file = open(signalp_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		if '#' == line[0]:
			continue
		line = line.strip('\n').split('\t')
		transcriptid = line[0].split('|')[0]
		if prev == transcriptid:
			continue
		prev = transcriptid
		while (True):
			error_check(annot_table[index], 9)
			if transcriptid == annot_table[index][0]:
				temp = 'sigP:' + line[3] + '^' + line[4] + '^' + line[5] + '^' + line[8]
				annot_table[index].append(temp)
				index += 1
				break
			if compare(transcriptid, annot_table[index][0]):
				annot_table[index].append('.')
				index += 1
			else:
				break
	while (index < length):
		error_check(annot_table[index], 9)
		annot_table[index].append('.')
		index += 1
	file.close()
	return annot_table


def add_tmhmm(annot_table, tmhmm_file):
	"""
	   Description:	Adds Column 10 - TMHMM
	   Input:		Annotation Table, TMHMM File
	   Return:		Annotation Table
	"""
	file = open(tmhmm_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		line = line.strip('\n').split('\t')
		transcriptid = line[0].split('|')[0]
		if prev == transcriptid:
			continue
		prev = transcriptid
		while (True):
			error_check(annot_table[index], 10)
			if transcriptid == annot_table[index][0]:
				temp = line[2] + '^' + line[4] + '^' + line[5]
				annot_table[index].append(temp)
				index += 1
				break
			if compare(transcriptid, annot_table[index][0]):
				annot_table[index].append('.')
				index += 1
			else:
				break
	while (index < length):
		error_check(annot_table[index], 10)
		annot_table[index].append('.')
		index += 1
	file.close()
	return annot_table


def add_rnammer(annot_table, rnammer_file):
	"""
	   Description:	Adds Column 11 - RNAMMER
	   Input:		Annotation Table, RNAMMER File
	   Return:		Annotation Table
	"""
	file = open(rnammer_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		line = line.strip('\n').split('\t')
		if prev == line[0]:
			continue
		prev = line[0]
		while (True):
			error_check(annot_table[index], 11)
			if line[0] == annot_table[index][0]:
				temp = line[7] + '^' + line[4] + '-' + line[5]
				annot_table[index].append(temp)
				index += 1
				break
			if compare(line[0], annot_table[index][0]):
				annot_table[index].append('.')
				index += 1
			else:
				break
	file.close()
	while (index < length):
		error_check(annot_table[index], 11)
		annot_table[index].append('.')
		index += 1
	return annot_table


def add_trembl_top_blastx(annot_table, unirefx_file):
	"""
	   Description:	Adds Column 12 - Uniref90BlastX
	   Input:		Annotation Table, Uniref90BlastX File
	   Return:		Annotation Table
	"""
	file = open(unirefx_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		line = line.strip('\n').split('\t')
		if prev == line[0]:
			continue
		prev = line[0]
		while (True):
			error_check(annot_table[index], 12)
			if line[0] == annot_table[index][0]:
				annot_table[index].append(line[1]) ###could change to 12###
				index += 1
				break
			if compare(line[0], annot_table[index][0]):
				annot_table[index].append('.')
				index += 1
			else:
				break
	while (index < length):
		error_check(annot_table[index], 12)
		annot_table[index].append('.')
		index += 1
	file.close()
	return annot_table


def add_trembl_top_blastp(annot_table, unirefp_file):
	"""
	   Description:	Adds Column 12 - Uniref90BlastP
	   Input:		Annotation Table, Uniref90BlastP File
	   Return:		Annotation Table
	"""
	file = open(unirefp_file)
	index = 0
	length = len(annot_table)
	prev = ''
	for line in file:
		if index == length:
			break
		line = line.strip('\n').split('\t')
		transcriptid = line[0].split('|')[0]
		if prev == transcriptid:
			continue
		prev = transcriptid
		while (True):
			error_check(annot_table[index], 13)
			if transcriptid == annot_table[index][0]:
				annot_table[index].append(line[12])
				index += 1
				break
			if compare(transcriptid, annot_table[index][0]):
				annot_table[index].append('.')
				index += 1
			else:
				break
	while (index < length):
		error_check(annot_table[index], 13)
		annot_table[index].append('.')
		index += 1
	file.close()
	for index in range(len(annot_table)):
		error_check(annot_table[index], 14)
		annot_table[index] += ['.'] * 15
	return annot_table