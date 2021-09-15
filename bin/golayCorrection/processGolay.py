from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *
from golay import *
import commands
import sys

out_handle = open(str(sys.argv[1]) + ".corrected","w")

for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
  res = decode(str(seq_record.seq))
  if res[0] != None:
    SeqIO.write(SeqRecord(Seq(res[0], SingleLetterAlphabet()), id=seq_record.id, description=""), out_handle, "fasta")
