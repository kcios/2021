import os
import pandas as pd
import sys
from t3_findvdjum import build_df
import concurrent.futures
from math import floor

dir = '/work/pi_gblanck/Konrad/DLBCL_RNA/DLBCLRNAbams_Results/' #str(sys.argv[1])
max_size = 40_000.0
receptor_list = 'TRA TRB TRD TRG IGH IGK IGL TRA_UM TRB_UM TRD_UM TRG_UM IGH_UM IGK_UM IGL_UM'.split()
receptor_counts = {}
print("Counting reads")

def findfiles_andmakedf(path,receptor):
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    print(f"Found {len(files)} {receptor} files in {path}")
    return build_df(path=path,files=files,receptor=receptor)

with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
    futures= []
    for receptor in receptor_list:
        print(f"Working on {receptor}")
        path = dir+receptor+'/'
        
        futures.append(executor.submit(findfiles_andmakedf, path=path,receptor=receptor))
      
    for future in concurrent.futures.as_completed(futures):
        fulldf = future.result()
        rid = fulldf.receptor_id
        df = fulldf.drop_duplicates(["Read"], keep='first')
        df = df.reset_index(drop=True)
        
        receptor_counts[rid] = len(df)
        print('finished->' + str(rid))
    
print(receptor_counts)
to_split = []
for rec in receptor_counts:
    ct = receptor_counts[rec]
    if  ct > max_size: # if the receptor has more reads than we want
        if (ct / max_size) < 1.45: # if the number of extra reads is less than 1/3 the max_size
            pass# dont bother splitting
        else:
            to_split.append(rec) # mark for splitting
to_process = []
for receptor in receptor_counts:
    if receptor not in to_split:
        to_process.append(receptor)
        

for rec in to_split:
    ct = receptor_counts[rec]
    splits = floor(ct/max_size)+1

    for i in range(splits):
        recpart = rec+'@'+str(splits)+'@'+str(i)
        to_process.append(recpart)
        
print(to_process)
process_string = '#'.join(to_process)

with open('t3_rectodo.txt',mode='w') as f:
    f.write(process_string)

numjobs = len(to_process) -1 
print(f"You need to use 0-{numjobs} for the array setting.")