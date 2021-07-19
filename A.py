import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import qth_survival_time
import os
import time
from scipy.stats import f_oneway, kruskal, ttest_ind
#########

class TooFewCases(Exception):
    pass

##########
#%% loading all the data and formatting it. Note this one is modifed and containts
## code to make sure 0 read patients are not excluded.

start= time.time()
print('loading clinical data')

rna_cpd = pd.read_csv('RNA_clinical.tsv',sep='\t',usecols=['case_id','days_to_last_follow_up','vital_status','gender','age_at_diagnosis',"ann_arbor_clinical_stage","international_prognostic_index","site_of_resection_or_biopsy",'tissue_or_organ_of_origin',"tumor_stage"],na_values=('\'--'))
dna_cpd = pd.read_csv('DNA_clinical.tsv',sep='\t',usecols=['case_id','days_to_last_follow_up','vital_status','gender','age_at_diagnosis',"ann_arbor_clinical_stage","international_prognostic_index","site_of_resection_or_biopsy",'tissue_or_organ_of_origin',"tumor_stage"],na_values=('\'--'))
print(f'loaded clinical - {time.time()-start:.5f}s')

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
############## modified for proper 0 receptor vs >1 comparison
dna_alldata = dna_cpd.merge(dna_vdjdata, right_on='Case ID',left_on='case_id',how='left',indicator=True).loc[:,['Case ID','Receptor','CDR3','VIDc','JIDc','VIDa','JIDa','days_to_last_follow_up','gender','alive','File ID','age_at_diagnosis',"ann_arbor_clinical_stage",'case_id','_merge']]
rna_alldata = rna_cpd.merge(rna_vdjdata, right_on='Case ID',left_on='case_id').loc[:,['Case ID','Receptor','CDR3','VIDc','JIDc','VIDa','JIDa','days_to_last_follow_up','gender','alive','File ID','age_at_diagnosis',"ann_arbor_clinical_stage"]]


dna_alldata.loc[dna_alldata['_merge'] == "left_only",'Case ID'] = dna_alldata.query('_merge == "left_only"')['case_id']
######## end mmodified.
#loading HLA data.
print(f'loading HLA data - {time.time()-start:.4f}s')
hladata = pd.read_excel('hla.xlsx',engine='openpyxl',sheet_name=0).drop(columns=[1,2])
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


print(f'loading physiochem data - {time.time()-start:.3f}s')
#pkdata  = pd.read_excel('VDJ_recoveries_fix.xlsx',sheet_name=3,engine='openpyxl')
#pkdata.to_csv('csv_vdj_pkdata.tsv',sep='\t')
dna_pkdata  = pd.read_csv('DNA_csv_vdj_pkdata.tsv', sep = '\t').iloc[:,3:].drop(columns=["Receptor","length","Case ID"]) # drop useless columns
dna_pkdata = dna_pkdata.drop_duplicates(subset="CDR3")
dna_alldatapk = dna_alldata.merge(dna_pkdata, on="CDR3",how="left")

rna_pkdata = pd.read_csv('RNA_VDJ_Recoveries_pk.csv')
rna_pkdata = rna_pkdata.drop_duplicates(subset="CDR3")
rna_alldatapk = rna_alldata.merge(rna_pkdata,on="CDR3",how="left")

print(f'physiochem data loaded - {time.time()-start:.4f}s')

del dna_pkdata, rna_pkdata

print(f'done - {time.time()-start:.4f}s')  
#cleaning up
del start


#%%    

#%% counting receptors, survival association, staging with receptor cts
# This cell just laods the receptor data and fills in 0 forpeople with none.
# Table 1 is generated here.
data = dna_alldata

q1 = data.dropna(axis=0,subset=['days_to_last_follow_up']).groupby(['Case ID','Receptor'],as_index=False).count()
q2 = q1[["Case ID","Receptor","CDR3"]].rename(columns={'CDR3':'ct'})

filled_rtype_data = {} # 240 have surivival data
receptors = ['TRA', 'IGK', 'TRD', 'TRB', 'TRG', 'IGH', 'IGL']
for rtype in receptors:
    a1 = data.query('Receptor == @rtype')
    
    have0  = data[~data["Case ID"].isin(a1["Case ID"])].dropna(axis=0,subset=['days_to_last_follow_up'])
    have0['Receptor'] = rtype
    have0['ct'] = 0
    have0 = have0[["Case ID","Receptor","ct"]]

    a2 = pd.concat([q2,have0],ignore_index=True,axis=0).drop_duplicates()
    
    a3 = a2.query('Receptor == @rtype')
    
    filled_rtype_data[rtype]  = a3
#%%
results = []
stage_res = []
w1 = data.drop_duplicates(subset=['Case ID'])[['Case ID','days_to_last_follow_up','alive',"ann_arbor_clinical_stage"]].dropna(axis=0,subset=['days_to_last_follow_up'])
for k in filled_rtype_data.keys():
    df = filled_rtype_data[k]
    print(f'######################## {k} ####################')
    
    e1 = df.merge(w1, on = 'Case ID')
    ####### 0 vs non zero ########
    g1 = e1[e1.ct>0]
    g0 = e1[e1.ct==0]
    
    targetT = g1.loc[:,'days_to_last_follow_up']
    targetE = g1.loc[:,'alive']
    
    restT = g0.loc[:,'days_to_last_follow_up']
    restE = g0.loc[:,'alive']
    
    avg_pos = df[df['ct']>0]['ct'].mean()
    
    res = logrank_test(targetT, restT, targetE, restE)
    #print(res)

    fig,ax = plt.subplots()


    kmf = KaplanMeierFitter()
    lab = k
    
    kmf.fit(g1['days_to_last_follow_up'],event_observed=g1['alive'],label=lab)
    kmf.plot(ax=ax,show_censors=True)
    
    pos_med_surv = qth_survival_time(.5,kmf)
    neg_med_surv=0
    if len(g0) >1:
        kmf.fit(g0['days_to_last_follow_up'],event_observed=g0['alive'],label=f'0 {k}')
        kmf.plot(ax=ax,show_censors=True)
        neg_med_surv = qth_survival_time(.5,kmf)
        
    results.append(pd.DataFrame([[f'{k} vs 0 {k}', res.p_value,len(targetT),len(restT),avg_pos,
                                  pos_med_surv,neg_med_surv]] ,
                                columns=["class","pval",
                                         "Present",'absent','positive_average_count',
                                         "pos_median_survival","neg_median_survival"]) )

    ## top third vs bottom ##############
    frac = .33
    length= int(len(e1)*frac)
    top = e1.nlargest(length, columns="ct",keep='all')
    bot = e1.nsmallest(length, columns="ct",keep='all')
    
    targetT = top.loc[:,'days_to_last_follow_up']
    targetE = top.loc[:,'alive']
    
    restT = bot.loc[:,'days_to_last_follow_up']
    restE = bot.loc[:,'alive']
    
    res = logrank_test(targetT, restT, targetE, restE)
    print(res)
    results.append(pd.DataFrame([[f'top {frac:.0%} {k} vs bot {frac:.0%} {k}', res.p_value]] ,columns=["class","pval"]) )
    fig,ax = plt.subplots()


    kmf = KaplanMeierFitter()
    lab = f'top {frac:.0%} {k}'
    
    kmf.fit(top['days_to_last_follow_up'],event_observed=top['alive'],label=lab)
    kmf.plot(ax=ax,show_censors=True)
    
    kmf.fit(bot['days_to_last_follow_up'],event_observed=bot['alive'],label=f'bottom {frac:.0%}')
    kmf.plot(ax=ax,show_censors=True)
    '''
    ### Clinical stage using top bot method
    binned = e1.copy()
    binned['top'] =  1 * (e1['Case ID'].isin(top['Case ID'] )    )
    binned['bot'] =  1 * (e1['Case ID'].isin(bot['Case ID'] )    )
    
    clin = binned.groupby('ann_arbor_clinical_stage',as_index=False).mean()
    clin2 = binned.groupby('ann_arbor_clinical_stage',as_index=False).std()

    stage_res.append( ( clin,clin2, k ) )
    '''
    
    

resultsdf = pd.DataFrame(columns=['class','pval','ntarget','nrest'])
for x in results:
   
    resultsdf =  resultsdf.append([x])

resultsdf.to_excel('ReceptorCountSurvival.xlsx', index=False)

'''  
for x in stage_res:
    print(x[2])
    print(x[0])
    print(x[1])
    
'''
#%%
#%% below not in paper


igkmerged = w1.merge(filled_rtype_data['IGL'], on = 'Case ID')

num = len(igkmerged.query('ct > 20')) / 221

quants = igkmerged.quantile([.2,.4,.5,.6,.8,.9])

#%% my staging system

data  = dna_alldata

w1 = data.drop_duplicates(subset=['Case ID'])[['Case ID','days_to_last_follow_up','alive',"ann_arbor_clinical_stage","age_at_diagnosis"]].dropna(axis=0,subset=['days_to_last_follow_up'])

data = w1.merge(filled_rtype_data['IGL'], on = 'Case ID')
tranums = w1.merge(filled_rtype_data['TRA'], on = 'Case ID').rename(columns={'ct':'tract'})

data = data.merge(tranums[['Case ID','tract']])
data['ratio'] =  data.ct / data.tract

m = data.loc[data['ratio'] != np.inf, 'ratio'].max()
data['ratio'].replace(np.inf,m,inplace=True)
data['ratio'].replace(np.nan,m,inplace=True) # not sure if min or max here


def assignClass(x):
    x = x['ct']
    if x <=2:
        return "A"
    elif x>2 and x<16:
        return "B"
    
    elif x >=16:
        return "C"
    
def assignClass2(x):
    x = x['ratio']
    if x <=.05:
        return "A"
    elif x>.05: # and x<.70:
        return "B"
    
    #elif x >=.70:
    #    return "C"    
    
data['kstage'] = data.apply(assignClass,axis=1)
data['kstage2'] = data.apply(assignClass2,axis=1)

td = []
alive = []
for stage in ["A","B","C"]:
    td.append(data[data.kstage2==stage]['days_to_last_follow_up'].values)
    alive.append(data[data.kstage2==stage]['alive'].values)
    
r1 = logrank_test(td[0],td[1],alive[0],alive[1])    
r2 = logrank_test(td[1],td[2],alive[1],alive[2])   
r3 = logrank_test(td[0],td[2],alive[0],alive[2]) 

ttest = ttest_ind(td[0],td[1])

Asurv = data[data.kstage2 =='A'].median()
Bsurv = data[data.kstage2 =='B'].median()

resanova = f_oneway(td[0],td[1],td[2])
reskru = kruskal(td[0],td[1],td[2])

bystage = data.groupby('kstage2').agg([np.mean, np.std])

import statsmodels.formula.api as smf

data = data.dropna(subset=["age_at_diagnosis"]).reset_index()


mod = smf.ols('days_to_last_follow_up ~ ct * tract + age_at_diagnosis + ratio ', data=data)

res = mod.fit()

print(res.summary())

#%% VIRAL CDR EPITOPES

data = dna_alldata

epitopes  = pd.read_csv('viral_epitopes.tsv',sep='\t')

filepi = epitopes.groupby('CDR3', as_index=False).count().query('Epitope == 1')['CDR3'].to_frame().merge(epitopes, on='CDR3')


dna_matches = filepi[filepi['CDR3'].isin(data['CDR3'])]['CDR3'].to_frame().merge(epitopes, on='CDR3')

merged = filepi.merge(data, on='CDR3').drop_duplicates()

counts = merged['Epitope gene'].value_counts()

aw = merged.drop_duplicates(subset=["Case ID","Epitope gene"]).groupby('Epitope gene',as_index=False).count()

epi_to_check = aw.nlargest(10,columns="Case ID")['Epitope gene']
epi_results = []
for epi in epi_to_check:
    have = merged[merged['Epitope gene'] == epi].drop_duplicates(subset="Case ID")
    nothave = merged[~merged['Case ID'].isin(have["Case ID"])].drop_duplicates(subset="Case ID")
    
    targetT = have.loc[:,'days_to_last_follow_up']
    targetE = have.loc[:,'alive']
    
    restT = nothave.loc[:,'days_to_last_follow_up']
    restE = nothave.loc[:,'alive']
    
    targetM = targetT.median()
    restM = restT.median()
    print(f'{epi} == {targetM} -- rest == {restM}')
    
    res = logrank_test(targetT, restT, targetE, restE)
    
    print(f"{epi} --- {res.p_value} ")
    
    epi_results.append([epi,res.p_value, targetM-restM, abs(len(targetT)*(targetM-restM)/100/res.p_value)])

epiresdf = pd.DataFrame(epi_results)
#%% viral cdr epitopes pastor set

data = dna_alldata

epitopes  = pd.read_csv('pastor.csv')
epitopes = epitopes[epitopes['Species']=='Human']
a = epitopes
a['CDR3'] = a['CDR3.alpha.aa']
b = epitopes
b['CDR3'] =  b['CDR3.beta.aa']

epitopes = pd.concat([a,b])
GP_TAR = 'Antigen_protein' # Antigen_protein
filepi = epitopes.groupby('CDR3', as_index=False).count().query('{0} <1000'.format(GP_TAR))['CDR3'].to_frame().merge(epitopes, on='CDR3')


dna_matches = filepi[filepi['CDR3'].isin(data['CDR3'])]['CDR3'].to_frame().merge(epitopes, on='CDR3')

merged = filepi.merge(data, on='CDR3').drop_duplicates()

counts = merged[GP_TAR].value_counts()

epi_ct = merged.drop_duplicates(subset=["Case ID",GP_TAR]).groupby(GP_TAR,as_index=False).count()

epi_to_check = epi_ct.nlargest(15,columns="Case ID")[GP_TAR]
epi_results = []
for epi in epi_to_check:
    have = merged[merged[GP_TAR] == epi].drop_duplicates(subset="Case ID")
    nothave = merged[~merged['Case ID'].isin(have["Case ID"])].drop_duplicates(subset="Case ID")
    
    targetT = have.loc[:,'days_to_last_follow_up']
    targetE = have.loc[:,'alive']
    
    restT = nothave.loc[:,'days_to_last_follow_up']
    restE = nothave.loc[:,'alive']
    
    targetM = targetT.median()
    restM = restT.median()
    print(f'{epi} == {targetM} -- rest == {restM}')
    
    res = logrank_test(targetT, restT, targetE, restE)
    
    print(f"{epi} --- {res.p_value} ")
    
    epi_results.append([epi,res.p_value, targetM-restM, abs(len(targetT)**2*(targetM-restM)/(res.p_value**4))])

epiresdf = pd.DataFrame(epi_results)