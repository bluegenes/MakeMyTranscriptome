#initialize Table

"""Space generator for empty files"""
def addSpace(annot_table):
    for index in range(len(annot_table)):
        annot_table[index].append(".")
    return annot_table

"""Compares if trID from a file is > trID in the row of the annot_table"""
def compare(a, b):
    a = [int(x[1:]) for x in (a.split("_"))]
    b = [int(x[1:]) for x in (b.split("_"))]
    for x, y in zip(a, b):
        if x > y:
            return True
    return False

"""Creates dictionary for transcript ID -> gene ID conversion"""
def createConversionTrIDtoGID(conversion_file): #conversion file from gene_trans_map
    convDt = {}
    file = open(conversion_file)
    for line in file:
        line = line.strip("\n").split()
        first = line[0]
        second = line[1]
        convDt[second] = first
    file.close()
    return convDt

"""Columns 0,1,2 - fasta & gene_trans_map  TranscriptID  GeneID   Transcript Length"""
def addFasta(fastaFile, transMapFile):
    annot_table = []
    file = open(fastaFile)
    trID = ""
    length = 0
    geneID = createConversionTrIDtoGID(transMapFile)
    for line in file:
        if ">" == line[0]:
            if length > 0:
                annot_table += [[trID, geneID[trID], str(length)]]
            length = 0
            trID = line.split(" ")[0][1:]
        else:
            length += len(line.strip())
    annot_table += [[trID, geneID[trID], str(length)]]
    file.close()
    return annot_table

"""Column 3,4 - swissprot_blastX    swissprot Length"""
def addSPTopBlastX(annot_table, swissprotXFile):
    file = open(swissprotXFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        line = line.strip("\n").split("\t")
        if prev == line[0]:
            continue
        prev = line[0]
        while (index < length):
            if line[0] == annot_table[index][0]:
                annot_table[index].append(line[12])
                annot_table[index].append(line[14])
                index += 1
                break
            #if sum([ord(c) for c in line[0]]) > sum([ord(c) for c in annot_table[index][0]]):
            if compare(line[0], annot_table[index][0]):
                annot_table[index].append(".")
                annot_table[index].append(".")
                index += 1
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table

"""Columns 5,6 - ProtID     ProtCoords"""
def addProtID(annot_table, transFile):
    file = open(transFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        if line[0] != ">":
            continue
        line = line.strip("\n").split()
        trID = line[8].split(":")[0]
        if prev == trID:
            continue
        prev = trID
        while (index < length):
            if trID == annot_table[index][0]:
                protID = line[0][1:]
                protCoord = line[8].split(":")[1].replace("(", "[").replace(")", "]")
                annot_table[index].append(protID)
                annot_table[index].append(protCoord)
                index += 1
                break
            #if sum([ord(c) for c in trID]) > sum([ord(c) for c in annot_table[index][0]]):
            if compare(trID, annot_table[index][0]):
                annot_table[index].append(".")
                annot_table[index].append(".")
                index += 1
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table

"""Column 7 - swissprot_blastP"""
def addSPTopBlastP(annot_table, swissprotPFile):
    file = open(swissprotPFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        line = line.strip("\n").split("\t")
        trID = line[0].split("|")[0]
        if prev == trID:
            continue
        prev = trID
        while (index < length):
            if trID == annot_table[index][0]:
                annot_table[index].append(line[1]) ###could change to 12###
                index += 1
                break
            if compare(trID, annot_table[index][0]):
                annot_table[index].append(".")
                index += 1
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table

"""Column 8 - PFAM --- needs to be checked"""
def addPfam(annot_table, pfamFile):
    file = open(pfamFile)
    index = 0
    length = len(annot_table)
    first = True
    for line in file:
        if index == length:
            break
        if "#" == line[0]:
            continue
        line = line.strip("\n").split()
        trID = line[3].split("|")[0]
        while (index < length):
            if trID == annot_table[index][0]:
                if first:
                    #temp = line[1] + "^" + line[0] + "^" + " ".join(line[22:]) + "^" + line[15] + "-" + line[16] + "^E:" + line[10]
                    temp = line[1] + "^" + " ".join(line[22:])
                    annot_table[index].append(temp)
                    first = False
                    break
                else:
                    #temp = line[1] + "^" + line[0] + "^" + " ".join(line[22:]) + "^" + line[15] + "-" + line[16] + "^E:" + line[10]
                    temp = line[1] + "^" + " ".join(line[22:])
                    annot_table[index][-1] = annot_table[index][-1] + "`" + temp
                    break
            if compare(trID, annot_table[index][0]):
                annot_table[index].append(".")
                index += 1
                first = True
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table

"""Column 9 - SignalP"""
def addSignalP(annot_table, signalpFile):
    file = open(signalpFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        if "#" == line[0]:
            continue
        line = line.strip("\n").split("\t")
        trID = line[0].split("|")[0]
        if prev == trID:
            continue
        prev = trID
        while (index < length):
            if trID == annot_table[index][0]:
                temp = "sigP:" + line[3] + "^" + line[4] + "^" + line[5] + "^" + line[8]
                annot_table[index].append(temp)
                index += 1
                break
            if compare(trID, annot_table[index][0]):
                annot_table[index].append(".")
                index += 1
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table

"""Column 10 - TMHMM"""
def addTmHMM(annot_table, tmhmmFile):
    file = open(tmhmmFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        line = line.strip("\n").split("\t")
        trID = line[0].split("|")[0]
        if prev == trID:
            continue
        prev = trID
        while (index < length):
            if trID == annot_table[index][0]:
                temp = line[2] + "^" + line[4] + "^" + line[5]
                annot_table[index].append(temp)
                index += 1
                break
            if compare(trID, annot_table[index][0]):
                annot_table[index].append(".")
                index += 1
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table

"""Column 11 - RNAMMER"""
def addRNAMMER(annot_table, rnammerFile):
    file = open(rnammerFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        line = line.strip("\n").split("\t")
        if prev == line[0]:
            continue
        prev = line[0]
        while (index < length):
            if line[0] == annot_table[index][0]:
                temp = line[7] + "^" + line[4] + "-" + line[5]
                annot_table[index].append(temp)
                index += 1
                break
            if compare(line[0], annot_table[index][0]):
                annot_table[index].append(".")
                index += 1
            else:
                break
    file.close()
    while (index < length):
        annot_table[index].append(".")
        index += 1
    return annot_table

"""Column 12 - uniref90_blastX"""
def addTrEMBLTopBlastX(annot_table, unirefXFile):
    file = open(unirefXFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        line = line.strip("\n").split("\t")
        if prev == line[0]:
            continue
        prev = line[0]
        while (index < length):
            if line[0] == annot_table[index][0]:
                annot_table[index].append(line[1]) ###could change to 12###
                index += 1
                break
            if compare(line[0], annot_table[index][0]):
                annot_table[index].append(".")
                index += 1
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table

"""Column 13 - uniref90_blastP"""
def addTrEMBLTopBlastP(annot_table, unirefPFile):
    file = open(unirefPFile)
    index = 0
    length = len(annot_table)
    prev = ""
    for line in file:
        if index == length:
            break
        line = line.strip("\n").split("\t")
        trID = line[0].split("|")[0]
        if prev == trID:
            continue
        prev = trID
        while (index < length):
            if trID == annot_table[index][0]:
                annot_table[index].append(line[12])
                index += 1
                break
            if compare(trID, annot_table[index][0]):
                annot_table[index].append(".")
                index += 1
            else:
                break
    while (index < length):
        annot_table[index].append(".")
        index += 1
    file.close()
    return annot_table



"""
out = open("testOut.txt", "w")
for row in Table:
    out.write("\t".join(row))
    out.write("\n")
out.close()
"""
