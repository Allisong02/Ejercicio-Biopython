from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Archivo de la carpeta en escritorio, se debe de cambiar a la dirección del nuevo archivo .gbk o .fasta a leer
filename =  "/mnt/c/Users/world/OneDrive/Escritorio/Allison/Ejercicio-Biopython/data/m_cold.fasta" 

def summarize_contents(filename):
	lista = []
	lista = os.path.split(filename)
	cadena = " "
	cadena = ("\nfile: "+ lista[1] + "\npath: " + lista[0])
	File_Extension = []
	File_Extension = os.path.splitext(filename)
	if(File_Extension[1]==".gbk"):
		type_file="genbank"
	else:
		type_file="fasta"
		pass
	records = list(SeqIO.parse(filename, type_file))
	cadena += ("\nnum_records: " + str(len(records)))
	cadena += ("\nrecord(s): ")
	for seq_record in SeqIO.parse(filename, type_file):
		cadena += ("\n- id: " + str(seq_record.id))
		cadena += ("\nname: " + seq_record.name)
		cadena += ("\ndescription: " + str(seq_record.description))
	return cadena
if __name__=="__main__":
	resultado = summarize_contents(filename)
	print(resultado)

#Comienzo de la nueva función agregada.
from Bio.Seq import Seq
seq1 = Seq("GTCAGCATA")
seq2 = Seq("GACTCATCA") 
def concatenate_and_get_reverse_complement(seq1, seq2):
    seqcon = seq1 + seq2 
    print(seqcon.reverse_complement())
concatenate_and_get_reverse_complement(seq1, seq2)
