```python
from Bio import SeqIO
import pandas as pd

DNA_Mus= [repr(seq_record.seq) for seq_record in SeqIO.parse("Mus_musculus_Actc1_sequence.txt",("fasta"))]

Mus_Acnt = (sum([round(((i.count("A"))/len(DNA_Mus)),0)for i in DNA_Mus])) # count A Mus

##print("{0}%is the A proportion in Mus musculus" .format(Mus_Acnt))


DNA_Rat= [repr(seq_record.seq) for seq_record in SeqIO.parse("Rattus_norvegicus_Acta1_sequence.txt",("fasta"))]

Rat_Acnt = (sum([round(((i.count("A"))/len(DNA_Rat)),0)for i in DNA_Rat])) # count Rat


DNA_Homo_sap= [repr(seq_record.seq) for seq_record in SeqIO.parse("Homo_sapiens_ACTA1_sequence.txt",("fasta"))]

Homo_Acnt = (sum([round(((i.count("A"))/len(DNA_Homo_sap)),0)for i in DNA_Homo_sap])) # count Human


DNA_Gorilla= [repr(seq_record.seq) for seq_record in SeqIO.parse("Gorilla_gorilla_ACTA1_sequence.txt",("fasta"))]

Gorilla_Acnt = (sum([round(((i.count("A"))/len(DNA_Gorilla)),0)for i in DNA_Gorilla])) # count Gorilla


DNA_F= [repr(seq_record.seq) for seq_record in SeqIO.parse("Felis_catus_ACTA1_sequence.txt",("fasta"))]

Felis_Acnt = (sum([round(((i.count("A"))/len(DNA_F)),0)for i in DNA_F])) # count Felis


DNA_Cavia_p= [repr(seq_record.seq) for seq_record in SeqIO.parse("Cavia_porcellus_ACTA1_sequence.txt",("fasta"))]

Cavia_Acnt = (sum([round(((i.count("A"))/len(DNA_Cavia_p)),0)for i in DNA_Cavia_p])) # count Cavia


DNA_Bos_t= [repr(seq_record.seq) for seq_record in SeqIO.parse("Bos_taurus_ACTA1_sequence.txt",("fasta"))]

Bos_Acnt = (sum([round(((i.count("A"))/len(DNA_Gorilla)),0)for i in DNA_Bos_t])) # count Bos


All_A_proportion = [Mus_Acnt, Rat_Acnt, Homo_Acnt, Gorilla Acnt, Felis_Acnt, Cavia_Acnt, Bos_Acnt]
##print(All_A_proportion)
SpNames = ["Mus_sp", "Rattus_sp", "Homo_sp", "Gorilla_sp", "Felis_sp", "Cavia_sp", "Bos_sp"]
SMR = [1.59, 0.84, 0.228, 0.20, 0.50, 2.13, 0.127] #weight specific metabolic rates

MyDict = {"SpName":SpNames,
         "SMR": SMR,
         "A_prop%":All_A_proportion}

df= pd.DataFrame(MyDict)
df.set_index("SpName", inplace = True)

print(df)

correl = round(df["SMR"].corr(df["A_prop%"]),2)
print("correlation between SMR and A_prop% is \n:{0}".format(correl))
```

    correlation between SMR and A_prop% is 
    :nan

    
