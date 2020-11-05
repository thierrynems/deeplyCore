#Premiere etape : il faut charger la librairie Pandas
import pandas
import os
import sys
import requests
#import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess # module pour executer les commandes Linux
#par convenance pure, nous modifions le nombre de lignes
# afficher dans les prints. L'idee est d'eviter que le tutoriel
#se resume a de multples affichages de longs tableaux
#vous pouvez modifier cette option a votre guise
pandas.options.display.max_rows = 10

datasetE_aa="Dataset/essential/degaa-p.dat"
datasetE_nt="Dataset/essential/degseq-p.dat"
datasetNE_aa="Dataset/nonEssential/degaa-np.dat"
datasetNE_nt="Dataset/nonEssential/degseq-np.dat"

#init de folder name
baseDirEssential="Dataset/Train/essential"
#test if baseDir exists
if os.path.isdir(baseDirEssential): 
	print("le repertoire existe")
else: 
	print("Creation du repertoire "+ baseDirEssential)
	os.mkdir(baseDirEssential);
#init de folder name
baseDirNonEssential="Dataset/Train/nonEssential"
#test if baseDir exists
if os.path.isdir(baseDirNonEssential): 
	print("le repertoire existe")
else: 
	print("Creation du repertoire "+ baseDirNonEssential)
	os.mkdir(baseDirNonEssential);

idListE_aa=list()
idListE_nt=list()
idListNE_aa=list()
idListNE_nt=list()
fastaSequences = SeqIO.parse(open(datasetE_aa), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	idListE_aa.append(name)

fastaSequences = SeqIO.parse(open(datasetE_nt), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	idListE_nt.append(name)

fastaSequences = SeqIO.parse(open(datasetNE_aa), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	idListNE_aa.append(name)

fastaSequences = SeqIO.parse(open(datasetNE_nt), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	idListNE_nt.append(name)

#charger le dataset de SRB 
dirEssentialSeq = sys.argv[1]
dirNonEssentialSeq = sys.argv[2]

datasetE_aa=os.path.join(dirEssentialSeq, "protein.fasta")
datasetE_nt=os.path.join(dirEssentialSeq, "gene.fasta")
datasetNE_aa=os.path.join(dirNonEssentialSeq, "protein.fasta")
datasetNE_nt=os.path.join(dirNonEssentialSeq, "gene.fasta")

seqaaDict_E=dict()
seqntDict_E=dict()
seqaaDict_NE=dict()
seqntDict_NE=dict()

fastaSequences = SeqIO.parse(open(datasetE_aa), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	seqaaDict_E[name]=sequence

fastaSequences = SeqIO.parse(open(datasetE_nt), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	seqntDict_E[name]=sequence

fastaSequences = SeqIO.parse(open(datasetNE_aa), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	seqaaDict_NE[name]=sequence

fastaSequences = SeqIO.parse(open(datasetNE_nt), 'fasta')
for fasta in fastaSequences:
	name, sequence = fasta.id, str(fasta.seq)
	seqntDict_NE[name]=sequence

e_EssentialFileAAPath= os.path.join(baseDirEssential, "degaa-p.dat")
e_EssentialFileNTPath= os.path.join(baseDirEssential, "degseq-p.dat")
e_EssentialFileAA=open(e_EssentialFileAAPath, 'a')
e_EssentialFileNT=open(e_EssentialFileNTPath, 'a')
#Essential
index=0
for cle in seqntDict_E.keys():
	if cle in seqaaDict_E:
		print(cle)
		sequenceProt=seqaaDict_E[cle].upper()
		sequenceGene=seqntDict_E[cle].upper()
		gene_Locus=idListE_nt[index]
		recordAA=SeqRecord(Seq(sequenceProt), id=str(gene_Locus),description="")
		recordNT=SeqRecord(Seq(sequenceGene), id=str(gene_Locus),description="")
		SeqIO.write(recordAA, e_EssentialFileAA, "fasta")
		SeqIO.write(recordNT, e_EssentialFileNT, "fasta")
		index=index+1
#nonEssential
e_EssentialFileAAPath= os.path.join(baseDirNonEssential, "degaa-np.dat")
e_EssentialFileNTPath= os.path.join(baseDirNonEssential, "degseq-np.dat")
e_EssentialFileAA=open(e_EssentialFileAAPath, 'a')
e_EssentialFileNT=open(e_EssentialFileNTPath, 'a')
index=0
for cle in seqntDict_NE.keys():
	if cle in seqaaDict_NE:
		print(cle)
		sequenceProt=seqaaDict_NE[cle].upper()
		sequenceGene=seqntDict_NE[cle].upper()
		gene_Locus=idListNE_nt[index]
		recordAA=SeqRecord(Seq(sequenceProt), id=str(gene_Locus),description="")
		recordNT=SeqRecord(Seq(sequenceGene), id=str(gene_Locus),description="")
		SeqIO.write(recordAA, e_EssentialFileAA, "fasta")
		SeqIO.write(recordNT, e_EssentialFileNT, "fasta")
		index=index+1

print("Fin du procesus de dataming")
