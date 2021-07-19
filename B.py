import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import qth_survival_time, median_survival_times
import os
import time
# This has teh V-J combo, VJ-HLA code
#########
class TooFewCases(Exception):
    pass
##########
#%% loading all the data and formatting it.

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


#%% Individual V or J segement impact on surival

def cleanSplitOut(category,value,df):
    df = df.dropna(axis=0,subset=['days_to_last_follow_up'])
    
    df2 = df.loc[df[category]==value,]
    df2 = df2.drop_duplicates(subset='Case ID')
    
    targetT = df2.loc[:,'days_to_last_follow_up']
    targetE  = df2.loc[:,'alive']

    if(len(targetT)<15):
        raise TooFewCases()
    
    restT = df.loc[~df['Case ID'].isin(df2['Case ID']),:]\
                .drop_duplicates(subset='Case ID').loc[:,'days_to_last_follow_up']
                
    restE = df.loc[~df['Case ID'].isin(df2['Case ID']),:]\
                .drop_duplicates(subset='Case ID').loc[:,'alive']
    
    
    return (targetT,targetE,restT,restE)

def calculateLogRankForTarget(category,value, df ):#calculates log rank for a VJID vs everything but that thing
    p = cleanSplitOut(category, value, df)
    #if p == 'fail':
    #    return 'fail'
    
    targetT = p[0]
    targetE  = p[1]
    
    restT = p[2]
    restE = p[3]
    
    return logrank_test(targetT, restT, targetE, restE )


#a = calculateLogRankForTarget('VIDc', value='TRAV12-1*01', df=alldata)    
#print(a.p_value)

#%%  generating pvalues for each V/J ID alone and its impact on surival
#adjust these to either c or a to get the respective set
###########################
Vver = 'VIDc'
Jver = 'JIDc'
data = rna_alldata

##########################
allV = data[Vver].unique()
allJ = data[Jver].unique()


def doCat(cat, IDs,df):
    size = IDs.shape[0]
    pvals = np.zeros((size))
    
    for idx,ID in enumerate(IDs):
        try:
            logrank = calculateLogRankForTarget(cat, ID, df)
        except TooFewCases:
            pvals[idx] = 1 # ie, it will get ignored
        else:
            pvals[idx] = logrank.p_value
        if idx % 15 == 0:
            print(f'{idx/size:.1%}')
    pvalues = pd.DataFrame({ 'ID':IDs, 'pval': pvals,'type':cat}  )
    
    return pvalues
    
v = doCat(Vver, allV, data)
j = doCat(Jver, allJ, data)

indivpvals = pd.concat([v,j])
good_pvals = indivpvals[indivpvals['pval']<.05]

#%% making kapmier graphs for the good pvalued individual VJIDs
resultsdir = 'rna_individual'
try:
    os.mkdir(resultsdir)
except OSError as e:
    print(e)

def makeKM(cat,val,df):# category ie VID/JID , the VID ID, df - alldata
    logrank = calculateLogRankForTarget(cat, val, df)
    parts = cleanSplitOut(cat, val, df)
    
    fig,ax = plt.subplots()
    
    kmf = KaplanMeierFitter()
    
    kmf.fit(parts[0],event_observed=parts[1],label=val)
    kmf.plot(ax=ax,show_censors=True)
    st15 = qth_survival_time(.5,kmf)
    st13 = qth_survival_time(.33,kmf)
    
    kmf.fit(parts[2],event_observed=parts[3],label=f'All but {val}')
    kmf.plot(ax=ax)
    st25  = qth_survival_time(.5,kmf)
    st23  = qth_survival_time(.33,kmf)
    
    val = val.replace('*','.')
    val = val.replace('/','.')
    
    subdir = val[0:3] # IGL, TRA, etc...
    try:
        os.mkdir(f'{resultsdir}/{subdir}')
    except OSError as e:
        print(e)
    
    ax.set_ylabel('Survival')
    ax.set_xlabel('Days Since Last Follow-Up')
    ax.text( .1,.2 ,f'P:{logrank.p_value:.3f}'  , transform=ax.transAxes)
    ax.text( .1,.28 ,f'Test stat: {logrank.test_statistic:.3f}'  , transform=ax.transAxes)
    ax.text( .1, .1 , f'# pt: {len(parts[0])}' ,transform=ax.transAxes)
    ax.text( .35, .08 , f'50|33th% {st13:n}|{st15:n} --- All: {st23:n}|{st25:n}' ,transform=ax.transAxes)
    fig.savefig(f'{resultsdir}/{subdir}/{val}.png',bbox_inches='tight',dpi=500)

    
for idx, ID in good_pvals.iterrows():
    print(idx,'--',ID)
    makeKM(ID['type'],ID['ID'],data)

counts = data.drop_duplicates(subset = ["Case ID",Vver, Jver])
counts = counts.melt( id_vars=['Case ID'],value_vars=[Jver,Vver],var_name='receptN', value_name='ID').drop_duplicates().groupby('ID',as_index=False).count()
c2 = indivpvals.merge(counts, on='ID', how='left').drop(columns=['Case ID'])

c2.to_excel(f'{resultsdir}/{resultsdir}.xlsx',index=False)

#%% ######################################################################
#%% survival differency by vid-jid combo

#adjust these to either c or a(no allele) to get the respective set
Vver = 'VIDc'
Jver = 'JIDc'
data = rna_alldata

combos = data.drop_duplicates(subset=['Case ID',Vver,Jver])

counts = combos.groupby([Vver,Jver],as_index=False).count()

good_combos = counts[counts['alive']>=14].reset_index(drop=True) # minimum limit, note can include those without survival data

#%%
#%% prun -T out6.txt
def comboClean(v,j,df):
    """ Removes dupes, ones without survival data. TO be used when searching for V/J combos. """
    #df = df.dropna(axis=0,subset=['days_to_last_follow_up'])
    
    df2 = df.loc[(df[Vver]==v) & (df[Jver]==j)  , ]
    df2 = df2.drop_duplicates(subset='Case ID')
    
    targetT = df2.loc[:,'days_to_last_follow_up']
    targetE  = df2.loc[:,'alive']

    if(len(targetT)<15):
        raise TooFewCases()
    
    #df0 = df.drop_duplicates(subset='Case ID')
    df0 = data_dedup
    mask = df0['Case ID'].isin(df2['Case ID'])
    
    df3 = df0.loc[~mask,:].loc[:,['days_to_last_follow_up','alive']]

    restT = df3['days_to_last_follow_up']
    restE = df3['alive']
    '''
    if(('V9' in v) and  ('2-2' in j) ):# printing inputs for fig 1
        print(f"printing {v} {j}")
        a = df0.loc[~mask,:]
        a['pos'] =  False
        df2['pos'] = True
        toprint = pd.concat([a,df2])
        v = v.replace("*","+")
        j = j.replace("*","+")
        toprint.to_excel(f'RNA_VJ_fig1_input_{v}_{j}.xlsx',index=False)
    '''
    del df2
    if(len(restT)<15):
        raise TooFewCases()

    return (targetT,targetE,restT,restE)
#just testing the func.
#b = comboClean(good_combos.loc[0,'VIDc'],good_combos.loc[0,'JIDc'],alldata)

def comboCalculateLogRank(v,j,df):
    tT,tE,rT,rE = comboClean(v,j,df)
    return logrank_test(tT,rT,tE,rE)


def doComboCat(combs,df):
    
    pvals = np.zeros(len(combs))
    for idx, row in combs.iterrows():
        #print(row)
        try:
            logrank = comboCalculateLogRank(row[Vver],row[Jver], df)
        except TooFewCases:
            pvals[idx] = 1 # ie, it will get ignored
        else:
            pvals[idx] = logrank.p_value
        if idx % 70 == 0:
            print(f'{idx/len(combs):.1%}')
    #pvalues = pd.DataFrame({ 'ID':IDs, 'pval': pvals,'type':cat.split('c')[0]}  )
    pvalues = combs.copy()
    pvalues['pval'] = pvals
    return pvalues

'''
def calculate(tocalc):
    
    pvals = np.zeros(len(tocalc))
    for idx, row in tocalc.iterrows():
        #print(row)
        try:
            logrank = comboCalculateLogRank(row[Vver],row[Jver], data)
        except TooFewCases:
            pvals[idx] = 1 # ie, it will get ignored
        else:
            pvals[idx] = logrank.p_value
        if idx % 55 == 0:
            print(f'{idx/len(tocalc):.1%}')
    #pvalues = pd.DataFrame({ 'ID':IDs, 'pval': pvals,'type':cat.split('c')[0]}  )
    pvalues = tocalc.copy()
    pvalues['pval'] = pvals
    return pvalues

import concurrent.futures

def doComboCat(combs,df):
    
    df_split = np.array_split(combs, 6)
    pool = concurrent.futures.ProcessPoolExecutor(max_workers = 2)
    
    return pd.concat(pool.map(calculate, df_split))
    
'''
##################
data  = data.dropna(axis=0,subset=['days_to_last_follow_up'])
data.loc[:,'Case ID'] =  data.loc[:,'Case ID'].astype('category') # gives nice speed boost for some reason
data.loc[:,Vver] =  data.loc[:,Vver].astype('category')
data[Jver] =  data[Jver].astype('category')


#For avoiding recalculating this over and over in the main function
data_dedup = data.drop_duplicates(subset='Case ID')

good_combos_f = good_combos.copy() # data in here is meaningless due to not filtering out no surv data pts
good_combos_f['Case ID'] = good_combos_f['Case ID'].astype('category')

comboPvals = doComboCat(good_combos_f,data)

g_comboPvals = comboPvals.query('pval <.05')

#%% makes KM graphs for the good V/J combos we found above.
resultsdir2 = 'rna_V+J_combos_allele'
import matplotlib.ticker as mtick
plt.style.use('grayscale')

try:
    os.mkdir(resultsdir2)
except OSError as e:
    print(e)
survs ={}
##### THis is the final graph format 
def comboMakeKM(v,j,df):# category ie VID/JID , the VID ID, df - alldata
    #logrank = comboCalculateLogRank(v, j , df)
    parts = comboClean(v, j, df)
    
    fig,ax = plt.subplots(figsize=(9.5,6),dpi=200)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    #ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0)) #formats the prop as a %
    
    
    kmf = KaplanMeierFitter()
    lab = v+'@'+j
    kmf.fit(parts[0],event_observed=parts[1],label=lab)
    kmf.plot(ax=ax,show_censors=True,ci_show=False,
             censor_styles={'ms': 8, 'marker': 2,'zorder':2},linewidth=3.2,zorder=2)
    
    st15 = qth_survival_time(.5,kmf)
    st17 = qth_survival_time(.70,kmf)
    medt = median_survival_times(kmf)
    
    kmf.fit(parts[2],event_observed=parts[3],label=f'All but {lab}')
    kmf.plot(ax=ax,show_censors=True,ci_show=False,
             linestyle=(0,(1,.3)),censor_styles={'ms': 8, 'marker': 2,'zorder':1},
             linewidth=3.2,zorder=1)
    st25  = qth_survival_time(.5,kmf)
    st27  = qth_survival_time(.70,kmf)
    medr = median_survival_times(kmf)
    
    survs[(v,j)] = {'Trg-50th':st15,'Rst-50th':st25,'Trg-70':st17,
                    'Rst-70':st27,'Pt_N':len(parts[0]),
                    'Median_surv_targ':medt, "Median_surv_rest":medr   }

     
    lab = lab.replace('*','.')
    lab = lab.replace('/','.')
    
    subdir = v[0:3] # IGL, TRA, etc...
    try:
        os.mkdir(f'{resultsdir2}/{subdir}')
    except OSError as e:
        print(e)
    
    #ax.set_ylabel('Survival')
    ax.set_xlabel(None)
    ax.set_xticks([1000,3000,5000])
    ax.set_yticks([0,.5,1])
    ax.set_yticklabels([0,50,100])
    
    ax.tick_params(axis='x',labelsize=45,width=5,
                   length=22,labelrotation=-35)
    
    ax.tick_params(axis='y',labelsize=45,width=5,length=22)
    ax.set_ylim(0,1.02)
    ax.set_xlim(1,6050)
    #ax.text( .05,.16 ,f'P: {logrank.p_value:.3f}'  , transform=ax.transAxes)
    #ax.text( .1,.28 ,f'Test stat: {logrank.test_statistic:.3f}'  , transform=ax.transAxes)
    #ax.text( .05, .1 , f'# pt: {len(parts[0])}' ,transform=ax.transAxes)
    #ax.text( .22, .08 , f'50|33th% {st13:n}|{st15:n} --- All: {st23:n}|{st25:n}' ,transform=ax.transAxes)
    ax.get_legend().remove()
    fig.savefig(f'{resultsdir2}/{subdir}/{lab}.png',bbox_inches='tight',dpi=450)

    
    
for idx, row in g_comboPvals.iterrows():
    #if ('A' not in row[Vver]) and ('B' not in row[Vver]): # quit plotting not tra/b
    #    continue
    print(idx,'--',row)
    comboMakeKM(row[Vver],row[Jver],data)
    

survdf = pd.DataFrame.from_dict(survs,orient='index')\
                    .reset_index().rename(columns={'level_0':Vver,'level_1':Jver})


c2 = comboPvals.loc[:,[Vver,Jver,'pval']] 
c3 = c2.merge(survdf, on=[Vver,Jver],how='left')
c3.to_excel(f'{resultsdir2}/{resultsdir2}_v3.xlsx',index=False)

#%% looking at HLA + a certain VJ
# In paper.

Vver = 'VIDa'
Jver = 'JIDa'
data = rna_alldatahla

#stacking HLAs of each type into 1 column
hlatypes = ['A', 'A1', 'B', 'B1', 'C', 'C1', 'DP', 'DP1', 'DQ', 'DQ1', 'DR', 'DR1']
df = data.loc[ data['Receptor'].isin(["TRA","TRB"]) ,:]
df = df.melt(id_vars='Case ID', value_vars=hlatypes,value_name='hlas').drop('variable',axis=1)
df.drop_duplicates(inplace=True)
df = df.merge(data, on = 'Case ID')

df = df[df['Receptor'].isin(['TRA','TRB'])]

cts = df.drop_duplicates(subset=[Vver,Jver,'Case ID','hlas']).groupby([Vver,Jver,'hlas'],as_index=False).count()

good_combos = cts.query('days_to_last_follow_up >= 15').reset_index()

#%%

def vjhlaClean(v,j,hla,df,mode='b'):
    
    if mode=='b':
        df2 = df.loc[(df[Vver]==v) & (df[Jver]==j)&(df['hlas']==hla)  ,]
    elif mode =='v':
        df2 = df.loc[(df[Vver]==v) & (df['hlas']==hla)  ,]
    elif mode=='j':
        df2 = df.loc[ (df[Jver]==j)&(df['hlas']==hla)  ,]
    df2 = df2.drop_duplicates(subset='Case ID')
    
    targetT = df2.loc[:,'days_to_last_follow_up']
    targetE  = df2.loc[:,'alive']

    if(len(targetT)<15):
        raise TooFewCases()
    
    dfdd = df.drop_duplicates(subset='Case ID')
    mask = dfdd['Case ID'].isin(df2['Case ID'])
    del df2
    df3 = dfdd.loc[~mask,:].loc[:,['days_to_last_follow_up','alive']]
    #restT = df.loc[~mask,:].loc[:,'days_to_last_follow_up']
    #restE = df.loc[~mask,:].loc[:,'alive']
    restT = df3['days_to_last_follow_up']
    restE = df3['alive']
    return (targetT,targetE,restT,restE)
   
def vjhlaLogRank(v,j,hla,df,mode):
    tT,tE,rT,rE = vjhlaClean(v,j,hla,df,mode)
    #print(tT,tE,rT,rE)
    return logrank_test(tT,rT,tE,rE)


def vjdoHLACat(combos,df,mode='b'):
    
    pvals = np.zeros(len(combos))
    for idx,row in combos.iterrows():
        #print(idx)
        v = row[Vver]
        j= row[Jver]
        hla = row['hlas']
        try:
            logrank = vjhlaLogRank(v,j,hla, df,mode)
        except TooFewCases:
            pvals[idx] = 1 # ie, it will get ignored
        else:
            pvals[idx] = logrank.p_value
        if idx %15 ==0:
            print(f'{idx/len(combos):.1%}')
    
    #pvalues = pd.DataFrame({ 'ID':IDs, 'pval': pvals,'type':cat.split('c')[0]}  )
    pvalues = combos.copy()
    pvalues['pval'] = pvals
    return pvalues
procdf = df[[Vver,Jver,'Case ID','alive','days_to_last_follow_up','hlas']].dropna(axis=0,subset=['days_to_last_follow_up'])
procdf2 = procdf.drop_duplicates()
procdf3 = procdf2.groupby([Vver,Jver,'hlas'],as_index=False).count()#.query('alive>8')
procdf4 = df[ df[Vver].isin(procdf3[Vver]) & df[Jver].isin(procdf3[Jver]) & df['hlas'].isin(procdf3['hlas']) ]
procdf5= procdf4[[Vver,Jver,'Case ID','alive','days_to_last_follow_up','hlas']].dropna(axis=0,subset=['days_to_last_follow_up']).drop_duplicates()
procdf5[['Case ID',Vver,Jver,'hlas']] = procdf5[['Case ID',Vver,Jver,'hlas']].astype('category')
good_combos_f = good_combos[['Case ID',Vver,Jver,'hlas']].astype('category')

#comboPvalsJ = vjdoHLACat(good_combos_f,procdf5,mode='j')
#comboPvalsV = vjdoHLACat(good_combos_f,procdf5,mode='v')
comboPvalsB = vjdoHLACat(good_combos_f,procdf5,mode='b')

#comboPvalsJ_f = comboPvalsJ.query('pval <.05')
#comboPvalsV_f = comboPvalsV.query('pval <.05')
comboPvalsB_f = comboPvalsB.query('pval <.05')

#%% plotting HLA + V + J combo

hlavjdir= 'rna_hlavj_NOallele'


plt.style.use('grayscale')

try:
    os.mkdir(hlavjdir)
except OSError as e:
    print(e)
survs = {}
def vjhlaMakeKM(v,j,hla,df,mode):# category ie VID/JID , the VID ID, df - alldata
    #logrank = vjhlaLogRank(v,j,hla, df,mode)
    parts = vjhlaClean(v,j,hla, df,mode)
    
    fig,ax = plt.subplots(figsize=(9.5,6),dpi=200)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    
    kmf = KaplanMeierFitter()
    lab = v +' '+  j +' '+ hla + ' mode:' + mode
    
    kmf.fit(parts[0],event_observed=parts[1],label=lab)
    kmf.plot(ax=ax,show_censors=True,ci_show=False,censor_styles={'ms': 8, 'marker': 2,'zorder':2},linewidth=3,zorder=2)
    st15 = qth_survival_time(.5,kmf)
    st13 = qth_survival_time(.33,kmf)    
    
    kmf.fit(parts[2],event_observed=parts[3],label=f'All but {lab}')
    kmf.plot(ax=ax,show_censors=True,ci_show=False,linestyle=(0,(1,.3)),censor_styles={'ms': 8, 'marker': 2,'zorder':1},linewidth=3,zorder=1)

    st25  = qth_survival_time(.5,kmf)
    st23  = qth_survival_time(.33,kmf)

    survs[(v,j,hla)] = {'Trg-50th':st15,'Rst-50th':st25,'Trg-33':st13,'Rst-33':st23,'Pt_N':len(parts[0])}    

    
    lab = lab.replace('*','+')
    lab = lab.replace('/','.')
    lab = lab.replace(':','$')
    
    subdir = v[0:3]
    
    try:
        os.mkdir(f'{hlavjdir}/{mode}')
    except OSError as e:
        print(e)

    try:
        os.mkdir(f'{hlavjdir}/{mode}/{subdir}')
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
    #ax.text( .05,.16 ,f'P: {logrank.p_value:.3f}'  , transform=ax.transAxes)
    #ax.text( .1,.28 ,f'Test stat: {logrank.test_statistic:.3f}'  , transform=ax.transAxes)
    #ax.text( .05, .1 , f'# pt: {len(parts[0])}' ,transform=ax.transAxes)
    #ax.text( .22, .08 , f'50|33th% {st13:n}|{st15:n} --- All: {st23:n}|{st25:n}' ,transform=ax.transAxes)
    #ax.legend(None)
    ax.get_legend().remove()
    fig.savefig(f'{hlavjdir}/{mode}/{subdir}/{lab}.png',bbox_inches='tight',dpi=450)

'''  
for idx, row in comboPvalsV_f.iterrows():
    print(idx,'--',row)
    vjhlaMakeKM(row[Vver],row[Jver],row['hlas'],procdf5,'v')
    
for idx, row in comboPvalsJ_f.iterrows():
    print(idx,'--',row)
    vjhlaMakeKM(row[Vver],row[Jver],row['hlas'],procdf5,'j')
'''    
for idx, row in comboPvalsB_f.iterrows():
    print(idx,'--',row)
    vjhlaMakeKM(row[Vver],row[Jver],row['hlas'],procdf5,'b')

survdf_hlaB = pd.DataFrame.from_dict(survs,orient='index').reset_index().rename(columns={'level_0':Vver,'level_1':Jver,'level_2':'hlas'})


Bdat = comboPvalsB.merge(survdf_hlaB, on=[Vver,Jver,'hlas'],how='left')

#comboPvalsV.to_excel(f'{hlavjdir}/{hlavjdir}V.xlsx')
#comboPvalsJ.to_excel(f'{hlavjdir}/{hlavjdir}J.xlsx')
Bdat.to_excel(f'{hlavjdir}/{hlavjdir}B.xlsx')
#%%


#%% Below is looking at HLA alone. Functional 

#%% Lookig at HLA type and its influence of survival.

data = rna_alldatahla


hlatypes = ['A', 'A1', 'B', 'B1', 'C', 'C1', 'DP', 'DP1', 'DQ', 'DQ1', 'DR', 'DR1']
df = data.copy()

groupedhla = []
for i in hlatypes:
    
    groupedhla.append(data.groupby(i,as_index=False).count())
    
combinedGropedHLA = []
for i in range(0,len(groupedhla),2):
    a = groupedhla[i]
    b = groupedhla[i+1]
    
    c = a.merge(b, left_on=hlatypes[i],right_on=hlatypes[i+1])
    c['n'] = c['Case ID_x'] + c['Case ID_y']
    c=c.query('n>10')
    
    hlas = list(c[hlatypes[i]+'_x'].values)
    combinedGropedHLA.extend(hlas)
    
    
hlaToCheck = set(combinedGropedHLA) # remove dupes
print('finished finding sufficent n combos')

def combineHLAs(x):
    vals =[]
    for hla in hlatypes:
        vals.append(x[hla])
    return vals

df['hla'] = df.apply(combineHLAs,axis=1)
print('combined hlas into a list')

def hlaClean(hla,df):
    df = df.dropna(axis=0,subset=['days_to_last_follow_up'])
    #print(df['hla'])
    print(hla)
    # this line selects the rows that contain the hla in its list
    df2 = df.loc[df.hla.apply(lambda x: bool(set(x) & set([hla]))) ,]
    #print(df2)
    df2 = df2.drop_duplicates(subset='Case ID')
    
    targetT = df2.loc[:,'days_to_last_follow_up']
    targetE  = df2.loc[:,'alive']
    
    restT = df.loc[~df['Case ID'].isin(df2['Case ID']),:].drop_duplicates(subset='Case ID')
    restT = restT.loc[:,'days_to_last_follow_up']
    restE = df.loc[~df['Case ID'].isin(df2['Case ID'])].drop_duplicates(subset='Case ID')
    restE =  restE.loc[:,'alive']
    
    if(len(targetT)<5):
        #print( ' failed.')
        raise TooFewCases()
        #return 'fail'
    #print(value, ' passed')
    return (targetT,targetE,restT,restE)

#x = hlaClean('DPB1*02:01',df)

def hlaLogRank(hla,df):
    tT,tE,rT,rE = hlaClean(hla,df)
    #print(tT,tE,rT,rE)
    return logrank_test(tT,rT,tE,rE)


def doHLACat(hlas,df):
    LEN = len(hlas)
    pvals = np.zeros(LEN)
    for idx,hla in enumerate(hlas):
        #print(row)
        try:
            logrank = hlaLogRank(hla, df)
        except TooFewCases:
            pvals[idx] = 1 # ie, it will get ignored
        else:
            pvals[idx] = logrank.p_value
            
        if (idx % LEN ==0):
            print(f'{idx/len:.1%}')
            
    
    #pvalues = pd.DataFrame({ 'ID':IDs, 'pval': pvals,'type':cat.split('c')[0]}  )
    pvalues = pd.DataFrame({'HLA':list(hlas)})
    pvalues['pval'] = pvals
    return pvalues
res = doHLACat(hlaToCheck,df)
comboPvals = res.query('pval <.05')
#%%
#%% ploting hlas
hladir= 'rna_hla'
try:
    os.mkdir(hladir)
except OSError as e:
    print(e)

def hlaMakeKM(hla,df):# category ie VID/JID , the VID ID, df - alldata
    logrank = hlaLogRank(hla, df)
    parts = hlaClean(hla, df)
    
    fig,ax = plt.subplots()
    
    kmf = KaplanMeierFitter()
    lab = hla
    kmf.fit(parts[0],event_observed=parts[1],label=lab)
    kmf.plot(ax=ax,show_censors=True)
    
    
    kmf.fit(parts[2],event_observed=parts[3],label=f'All but {lab}')
    kmf.plot(ax=ax)
    
    
    lab = lab.replace('*','+')
    lab = lab.replace('/','.')
    lab = lab.replace(':','$')
    
    ax.set_ylabel('Survival')
    ax.set_xlabel('Days Since Last Follow-Up')
    ax.text( .1,.2 ,f'P:{logrank.p_value:.3f}'  , transform=ax.transAxes)
    ax.text( .1,.28 ,f'Test stat: {logrank.test_statistic:.3f}'  , transform=ax.transAxes)
    
    fig.savefig(f'{hladir}/{lab}.png',bbox_inches='tight',dpi=400)

for idx, row in comboPvals.iterrows():
    print(idx,'--',row)
    hlaMakeKM(row['HLA'],df)
    
res.to_excel(f'{hladir}/hla.xlsx', index=False)
#%%