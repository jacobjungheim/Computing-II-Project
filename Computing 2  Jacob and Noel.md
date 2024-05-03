```python
from Bio import SeqIO
import pandas as pd

DNA_Mus= [repr(seq_record.seq) for seq_record in SeqIO.parse("Mus_musculus_Actc1_sequence.txt",("fasta"))]

Mus_Acnt = (sum([round(((i.count("A"))/len(DNA_Mus)),0)for i in DNA_Mus])) # count A Mus

##print("{0}%is the A proportion in Mus musculus" .format(Mus_Acnt))


All_A_proportion = [Mus_Acnt]
##print(All_A_proportion)
SpNames = ["Mus_sp"]
SMR = [1.59] #weight specific metabolic rates

MyDict = {"SpName":SpNames,
         "SMR": SMR,
         "A_prop%":All_A_proportion}

df= pd.DataFrame(MyDict)
df.set_index("SpName", inplace = True)

##print(df)

correl = round(df["SMR"].corr(df["A_prop%"]),2)
#print("correlation between SMR and A_prop% is \n:{0}".format(correl))
```

    correlation between SMR and A_prop% is 
    :nan
    
