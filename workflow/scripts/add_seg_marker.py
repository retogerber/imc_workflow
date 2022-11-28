import pandas as pd
import re

df_name_in = 'Cirrhosis-TMA-5_01062022_A5_5_summary.csv'
df_name_out = "test.csv"
regex_used = "E-cadherin"

df_name_in = snakemake.input[0]
df_name_out = snakemake.output[0]
regex_used = snakemake.params[0]


df=pd.read_csv(df_name_in)

newdeepcell=df["deepcell"].fillna('-1').astype(int).astype(str)
segmatches = [ re.search(regex_used,stri) for stri in df["name"] ]
for i in range(len(segmatches)):
    if segmatches[i] != None:
        newdeepcell[i] = "2"
    if newdeepcell[i] == "-1":
        newdeepcell[i] = ""
print(newdeepcell)

df["deepcell"] = newdeepcell

df.to_csv(df_name_out,index=False)

