import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import qth_survival_time
import os
import time
# This has the sub-sequence stuff, with the motifs
#########
class TooFewCases(Exception):
    pass
##########
#%% loading all the data and formatting it. coped from t6_vdj

start= time.time()
print('loading clinical data')

rna_cpd = pd.read_csv('RNA_clinical.tsv',sep='\t',
                      usecols=['case_id','days_to_last_follow_up','vital_status',
                               'gender','age_at_diagnosis',"ann_arbor_clinical_stage",
                               "international_prognostic_index",
                               "site_of_resection_or_biopsy",
                               'tissue_or_organ_of_origin',"tumor_stage"],na_values=('\'--'))
dna_cpd = pd.read_csv('DNA_clinical.tsv',sep='\t',
                      usecols=['case_id','days_to_last_follow_up',
                               'vital_status','gender',
                               'age_at_diagnosis',"ann_arbor_clinical_stage",
                               "international_prognostic_index",
                               "site_of_resection_or_biopsy",'tissue_or_organ_of_origin',
                               "tumor_stage"],na_values=('\'--'))
print(f'loaded clinical - {time.time()-start:.5f}s')

## alive dead changed!!! original used was alive  = 1
rna_cpd['alive'] = 1
rna_cpd.loc[rna_cpd['vital_status']=='Alive','alive'] = 0

dna_cpd['alive'] = 1
dna_cpd.loc[dna_cpd['vital_status']=='Alive','alive'] = 0
print(f'processed clinical- {time.time()-start:.4f}s')

print(f'loading vdj - {time.time()-start:.4f}s')
#vdjdata  = pd.read_excel('VDJ_recoveries_fix.xlsx',sheet_name=1,engine='openpyxl')
dna_vdjdata = pd.read_csv('DNA_csv_vdj_data.tsv',sep='\t')
rna_vdjdata  = pd.read_csv('RNA_VDJ_Rec_2.csv')

print(f'loaded vdj - {time.time()-start:.5f}s')
dna_vdjdata['VIDc'] = dna_vdjdata['VID'].str.split('|').str[1]
dna_vdjdata['JIDc'] = dna_vdjdata['JID'].str.split('|').str[1]
dna_vdjdata['File ID'] = dna_vdjdata['Filename'].str.split('_').str[0]

rna_vdjdata['VIDc'] = rna_vdjdata['VID'].str.split('|').str[1]
rna_vdjdata['JIDc'] = rna_vdjdata['JID'].str.split('|').str[1]
rna_vdjdata['File ID'] = rna_vdjdata['Filename'].str.split('_').str[0]

#removing alleles###
print(f'Creating allele free column - {time.time()-start:.5f}s')
dna_vdjdata['VIDa'] =  dna_vdjdata['VIDc'].str.split('*').str[0]
dna_vdjdata['JIDa'] =  dna_vdjdata['JIDc'].str.split('*').str[0]
rna_vdjdata['VIDa'] =  rna_vdjdata['VIDc'].str.split('*').str[0]
rna_vdjdata['JIDa'] =  rna_vdjdata['JIDc'].str.split('*').str[0]

print(f'finished cleaning V/J - {time.time()-start:.4f}s')
##########

dna_alldata = (dna_cpd.merge(dna_vdjdata, right_on='Case ID',left_on='case_id')
    .loc[:,['Case ID','Receptor','CDR3','VIDc','JIDc',
            'days_to_last_follow_up',
            'gender','alive','File ID','VIDa','JIDa']
         ]
    )
rna_alldata = (rna_cpd.merge(rna_vdjdata, right_on='Case ID',left_on='case_id')
    .loc[:,['Case ID','Receptor','CDR3','VIDc','JIDc',
            'days_to_last_follow_up','gender',
            'alive','File ID',
            'VIDa','JIDa']
         ]
    )
#loading HLA data.
print(f'loading HLA data - {time.time()-start:.4f}s')
hladata = pd.read_excel('hla.xlsx',sheet_name=0 ,engine='openpyxl' ).drop(columns=[1,2])
print(f'Read hla data - {time.time()-start:.5f}s')

hladata.loc[:,'File ID'] =  hladata.loc[:,'File ID'].str.strip()

dna_data= pd.read_csv('DNA_csv_vdj_data.tsv',sep='\t')
dna_data['File ID'] = dna_data['Filename'].str.split('_').str[0]
dna_data = dna_data.merge(hladata, on='File ID', how='inner')
dna_data = dna_data.loc[:,'Case ID':].drop(columns ='File ID').drop_duplicates()
dna_alldatahla = dna_alldata.merge(dna_data, on ='Case ID',how='inner')


rna_alldatahla = rna_alldata.merge(dna_data, on ='Case ID',how='inner')
print(f'Merged hla data - {time.time()-start:.3f}s')

del hladata, dna_data, dna_cpd, rna_cpd, dna_vdjdata, rna_vdjdata


# print(f'loading physiochem data - {time.time()-start:.3f}s')
#pkdata  = pd.read_excel('VDJ_recoveries_fix.xlsx',sheet_name=3,engine='openpyxl')
#pkdata.to_csv('csv_vdj_pkdata.tsv',sep='\t')

#dna_pkdata  = pd.read_csv('DNA_csv_vdj_pkdata.tsv', sep = '\t').iloc[:,3:].drop(columns=["Receptor","length","Case ID"]) # drop useless columns
#dna_pkdata = dna_pkdata.drop_duplicates(subset="CDR3")
#dna_alldatapk = dna_alldata.merge(dna_pkdata, on="CDR3",how="left")

#rna_pkdata = pd.read_csv('RNA_VDJ_Recoveries_pk.csv')
#rna_pkdata = rna_pkdata.drop_duplicates(subset="CDR3")
#rna_alldatapk = rna_alldata.merge(rna_pkdata,on="CDR3",how="left")

# print(f'physiochem data loaded - {time.time()-start:.4f}s')

#del dna_pkdata, rna_pkdata

print(f'done - {time.time()-start:.4f}s')  
#cleaning up
del start


#%%

#%% run to init the function, once.


def motifReplace(x):
    
    replacements =  {'P' : 'KR',
    'N' : 'DE',
    'A' : 'WFY',
    'H' : 'AVLIM',
    'X': 'GHCSTNQP'}
    
    cdr = x['CDR3']
    new = []
    for aa in cdr:
        for rep in replacements.keys():
            if aa in replacements[rep]:
                new.append(rep)
    return ''.join(new)

#%%

#%% looking for sub-sequences. This finds al thes sub-sequences in the chosen receptor.
v = 'TRBV13'
j = 'TRBJ2-5'
seq_mode= 'rna'

if seq_mode =='rna':
    procdat = rna_alldata
elif seq_mode =='dna':
    procdat = dna_alldata

#procdat =  procdat.query('(Receptor == "TRA")  or (Receptor == "TRB") ')
procdat =  procdat.query('(Receptor == "TRB") ')

msk = (procdat['VIDa'] ==v ) &  (procdat['JIDa'] == j) 
#msk = (procdat['VIDa'] == v) 
mm = procdat[msk].drop_duplicates(subset=['VIDc','JIDc','Case ID','CDR3'])\
                .dropna(subset=['days_to_last_follow_up'],axis=0)
                
mm['reduced'] = mm.apply(motifReplace,axis=1)
sslist = []
for i in range(5,15): # Range of lengths to test here!
    SS_LEN =i
    def get_all_ss(seq,l):
        slen = len(seq)
        l = int(l)
        maxstart = slen- l +1
        ss = []
        for i in range(0,maxstart):
            
            i = int(i)
            ss.append(seq[i:i+l])
        return ss
    
    tt = mm.reduced.apply(get_all_ss,args=(SS_LEN,))
    
    allsubseq = np.array(())
    for grp in tt:
        allsubseq= np.append(allsubseq,grp)
    
    subseqs = pd.Series(allsubseq).value_counts()[:2000] # Note value counts seems to return ina  random order each time.
    subseqs = subseqs.rename('from_CID_count')
    sslist.append((subseqs,i))
#%% Needs tobe run once to load things up.
data = procdat
data = data.dropna(axis=0,subset=['days_to_last_follow_up'])
data  = data.drop_duplicates(subset=["Case ID",'CDR3'])
data['reduced'] = data.apply(motifReplace,axis=1)

data.loc[:,'Case ID'] =  data.loc[:,'Case ID'].astype('category')
data_dedup = data.drop_duplicates(subset="Case ID")
#%%

def ssClean(s,df,mode='s'):
    
    if mode=='s':
        df2 = df.loc[(df['CDR3'].str.contains(s)) ,]
    elif mode =='m': # motif contains method
        df2 = df.loc[(df['reduced'].str.contains(s)) ,]
    elif mode=='match':
        df2 = df.loc[(df['reduced'] == s) ,]
    df2 = df2.drop_duplicates(subset='Case ID')
    
    targetT = df2.loc[:,'days_to_last_follow_up']
    targetE  = df2.loc[:,'alive']
    
    tlen =  len(targetT) 
    
    if(tlen <15):
        raise TooFewCases()
    
    #dfdd = df.drop_duplicates(subset='Case ID') #consider making global to increase perf
    dfdd = data_dedup
    mask = dfdd['Case ID'].isin(df2['Case ID'])
    
    df3 = dfdd.loc[~mask,:].loc[:,['days_to_last_follow_up','alive']]

    restT = df3['days_to_last_follow_up']
    restE = df3['alive']
    
    #printing out input for B
    if s == 'XHXXHXXAXX':
        print('printing example inputs')
        toprintout = pd.concat([dfdd.loc[~mask,:],df2])
        toprintout.to_excel('RNA_motif_example_XHXXHXXAXX.xlsx')
    
    if(len(restT)<15):
        raise TooFewCases()
    del df2
    return ((targetT,targetE,restT,restE),tlen)
   

def ssRank(s,df,mode):
    ret = ssClean(s,df,mode)
    tT,tE,rT,rE = ret[0]

    return (logrank_test(tT,rT,tE,rE),ret[1])


def ssCat(ss,df,mode='s'):
    targetList = []
    pvals = np.zeros(len(ss))
    for idx,s in enumerate( ss.index):
        
        try:
            logrank = ssRank(s ,df,mode)
        except TooFewCases:
            pvals[idx] = 1 # ie, it will get ignored
            targetList.append(1)
        else:
            pvals[idx] = logrank[0].p_value
            targetList.append(logrank[1])
        if idx %75 ==0:
            print(f'{idx/len(ss):.1%}')
    #print(targetList)
    pvalues = pd.DataFrame(ss)
    pvalues['pval'] = pvals
    pvalues['unique-CID-ct'] =  pd.Series(targetList).values
    #pvalues['goodness']  =    (pvalues['unique-CID-ct']**2) / pvalues['pval']**2.5
    return pvalues



resl = []
res_fl = []
PROCESS_MODE ='m'
for sitem in sslist:

    tocheck = sitem[0]
    res = ssCat( tocheck,data,mode=PROCESS_MODE) # change mode!
    
    res_f = res.query('pval<.05')
    sid = sitem[1]
    resl.append((res,sid))
    res_fl.append((res_f,sid))
    print(f'↟↟↟↟↟↟↟ {sid} done ↟↟↟↟↟↟↟↟')
#%% plotting subsewqs

plt.style.use('grayscale')

ssdir= f'{seq_mode}_motif_subsequence'
try:
    os.mkdir(ssdir)
except OSError as e:
    print(e)
survs ={}
def ssMakeKM(s,df,mode):# category ie VID/JID , the VID ID, df - alldata
    #logrank = ssRank(s, df,mode)[0]
    parts = ssClean(s, df,mode)[0]
    
    fig,ax = plt.subplots(figsize=(9.5,6),dpi=200)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
      
    kmf = KaplanMeierFitter()
    lab = s +' '+' '+' mode:' + mode
    
    #tl = range(0,6201,7) # RECENTLY ADDED
    kmf.fit(parts[0],event_observed=parts[1],label=lab)
    kmf.plot(ax=ax,show_censors=True,ci_show=False,
             censor_styles={'ms': 8, 'marker': 2,'zorder':2},linewidth=3,zorder=2)
    st15 = qth_survival_time(.5,kmf)
    st13 = qth_survival_time(.33,kmf)  
    st17 = qth_survival_time(.70,kmf) 
    
    kmf.fit(parts[2],event_observed=parts[3],label=f'All but {lab}')
    kmf.plot(ax=ax,show_censors=True,ci_show=False,
             linestyle=(0,(1,.3)),
             censor_styles={'ms': 8, 'marker': 2,'zorder':1},linewidth=3,zorder=1)

    st25  = qth_survival_time(.5,kmf)
    st23  = qth_survival_time(.33,kmf)
    st27 = qth_survival_time(.70,kmf)

    survs[s] = {'Trg-50th':st15,'Rst-50th':st25,'Trg-33':st13,'Rst-33':st23,
                'Trg-70':st17, 'Rst-70':st27,'Motif_Len':SS_LEN}    

    lab = lab.replace('*','+')
    lab = lab.replace('/','.')
    lab = lab.replace(':','$')

    try:
        os.mkdir(f'{ssdir}/{mode}')
    except OSError as e:
        print(e)


    ax.set_xlabel(None)
    ax.set_xticks([1000,3000,5000])
    ax.set_yticks([0,.5,1])
    ax.set_yticklabels([0,50,100])
    
    ax.tick_params(axis='x',labelsize=45,width=5,length=22,labelrotation=-35)
    ax.tick_params(axis='y',labelsize=45,width=5,length=22)
   
    ax.set_ylim(0,1.02)
    ax.set_xlim(1,6050)    
    
    '''
    ax.set_ylabel('Survival')
    ax.set_xlabel('Days Since Last Follow-Up')
    ax.text( .05,.16 ,f'P: {logrank.p_value:.3f}'  , transform=ax.transAxes)
    #ax.text( .1,.28 ,f'Test stat: {logrank.test_statistic:.3f}'  , transform=ax.transAxes)
    ax.text( .05, .1 , f'Pt#: {len(parts[0])}' ,transform=ax.transAxes)
    ax.text( .01, .02 , f'50|33th% {st13:n}|{st15:n} --- All: {st23:n}|{st25:n}' ,transform=ax.transAxes)
    ax.legend(loc=1)
    '''
    ax.get_legend().remove()
    fig.savefig(f'{ssdir}/{mode}/{SS_LEN}_{s}.png',bbox_inches='tight',dpi=400)

for g in res_fl:
    
    for idx, s in g[0].iterrows():
        print(s)
        SS_LEN = g[1]
        ssMakeKM(s.name,data,PROCESS_MODE)

survdf = pd.DataFrame.from_dict(survs,orient='index')\
                    .reset_index(level=0).rename(columns={'index':'seq'})
survdf['50diff'] =  survdf['Trg-50th'] -  survdf['Rst-50th']
survdf['33diff'] =  survdf['Trg-33'] -  survdf['Rst-33']

v  = v.replace('/','.')
tocct = []
for x in resl:
    tocct.append(x[0])
allres = pd.concat(tocct).reset_index(level=0).rename(columns={'index':'seq'})
allres = allres.merge(survdf, on='seq', how='left')
allres.to_excel(f'{ssdir}/{ssdir}_{v}_{j}.xlsx', index=False)


#%%
##############################################################################################
#%% Excel sheet combinator to check all combos. Not used.
folders_to_check = [ 'rna_ss_m' , 'dna_ss_m'  ]

from pathlib import Path

sstype = []
for fold in folders_to_check:
    path = Path(fold)
    
    excel_f =  list(path.glob('*.xlsx'))
    print(excel_f)
    groups = []
    for p in excel_f:
        if p.is_file():
            file = pd.read_excel(p)
            groups.append(file)
    sstype.append((groups  ,fold.split('_')[0]   ))

combined = []
for ty in sstype:
    name = ty[1]
    dfs = ty[0]
    combdf = pd.concat(dfs)
    combdf['type'] = name
    combdf = combdf.rename(columns={0:'seq'})
    #combdf = combdf.drop(columns=[combdf.columns[0]])
    combined.append(combdf)
    
fdf = pd.concat(combined).drop_duplicates().query('pval < .05')

matched = fdf.groupby('seq', as_index=False).count().query('type==2')

passing = matched.merge(fdf, on='seq')
