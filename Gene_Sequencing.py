## Gene Sequencing

### Importing packages and reading in data.
### Genome sequencing data is stored in the dataframe pn.
### Notice how each column has separate sequencing data.
### Our goal will be to find similarities between the sequences and characterize it.

import pandas as pd
import numpy as np

pn = pd.read_csv('SARS_CORONAVIRUS_NC_045512_sequence.fasta', header = None)[1:][0]
pn = pd.DataFrame(pn)
pn.columns = ['Genome Sequence']

## Find the length of the sequence in each row and store the results in a new
## column called ['Len']. Then, find the median length of these genome sequences and store it 
## in a variable called 'median_len'. Delete all rows that are not of this length.

def SeqLen(column):
    return len(column)

pn["Len"]=pn["Genome Sequence"].apply(SeqLen)
median_len = np.median(pn["Len"])
pn=pn[pn["Len"]==median_len]

## Create a new column titled 'Match?' with a 1/0 corresponding to whether or not the sequence
## contains the substring 'TAATTTAGGCATGCCTT'.

def match(column):
    if 'TAATTTAGGCATGCCTT' in column :
        return 1
    else :
        return 0
    
pn["match?"]=pn["Genome Sequence"].apply(match)

## Create a matrix called 'differences' storing the number of differences in the sequences between the rows.

def dif_seq(column1, column2):
    count = 0
    for i in range(len(column1)):
        if column1[i] != column2[i]:
            count+=1
        else : 
            count = count
    return count
differences = []
for j in range(1,len(pn["Genome Sequence"])):
    differences.append(dif_seq(pn["Genome Sequence"][1],pn["Genome Sequence"][j]))

## Find the two genome sequences that are most similar. Store them in a tuple called 'similar' with the format (row1, row2).

table_sim=np.zeros((len(pn["Genome Sequence"]),len(pn["Genome Sequence"])))
differences2=[]
min_val = 1000
simm=(0,0)

for i in range(1,len(pn["Genome Sequence"])+1):
        for j in range(1,len(pn["Genome Sequence"])+1):
                differences2.append(dif_seq(pn["Genome Sequence"][i],pn["Genome Sequence"][j]))
        table_sim[i-1]=differences2
        differences2=[]

for i in range(0,len(pn["Genome Sequence"])):
        for j in range(0,len(pn["Genome Sequence"])-1):
                if table_sim[i][j]==0:
                        table_sim[i][j]=1000

for i in range(0,len(pn["Genome Sequence"])):
        for j in range(0,len(pn["Genome Sequence"])-1):
                if table_sim[i][j]<min_val:
                        min_val=table_sim[i][j]
                        simm=(i,j)







