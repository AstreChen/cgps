### test for svm function in run_mulmds.py
import sys,os 
import numpy as np
import pandas as pd

def read_gmt_file(gmt):
    gset_des={}    #description of gene sets
    gset_genes = {}  # genes of gene sets
    gmt_file = open(gmt,'r')
    line = gmt_file.readline()
    try :
        while line :
            if '#' in line:
                line = gmt_file.readline()
                continue
            words = line.rstrip('\n').split('\t')

            if len(words)<3 :
                raise exception.GmtformatErr, "Error Message:\nBad gene set database format. Please check your gene set database format, only gmt file, a tab delimited file are allowed. \nAbout file format, please click www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats "

            if words[0] in gset_des.keys():
                continue
                print "Warning Message from KOBAS: \n These are at leaset two gene sets share the same gene set ID. Please check. We only use the gene set first emerges."

            gset_des[words[0]] = words[1]
            gset_genes[words[0]] = words[2:]
            line = gmt_file.readline()
        gmt_file.close()
    except IOError,e:
        print "Error Message:\nBad gene set database format. Please check your gene set database format, only gmt file are allowed."
        gmt_file.close()
        sys.exit(e)
    return gset_des, gset_genes


    #datadir = #save all 9 methods results
    #gmtf = # gmt file
    #svmfile =  # svm file

def run_svm_18dims(datadir,gmtf,outdir,svmfile='./data/cgps_model.pkl'):
    ### read svm model
    from sklearn.externals import joblib
    svc = joblib.load(svmfile)
    svc_features = [
        'cepa',
        'gage',
        'ganpa',
        'globaltest',
        'gsa',
        'gsea',
        'padog',
        'plage',
        'safe']  


    ## read gene sets ###
    gmt_dict = {}
    gset_des, gset_genes = read_gmt_file(gmtf)
    gset_df = pd.DataFrame({'GENE.SET':gset_des.keys(), 'TITLE':gset_des.values()})

    ### read table ###
    tbl = {}
    for md in svc_features:
        tbl[md] = pd.read_csv(datadir + md +'.tsv', sep='\t')
        if md == 'gage':
            tbl[md] = tbl[md].rename(columns={'geneset':'GENE.SET'})

    ### process table to svc feature vector ###
    md_rank={}
    rank_tb = gset_df
    for md in svc_features:
        if md == 'gage':
            pval = 'pval.min'
        elif md == 'plage':
            pval = 'P.Value'
        else:
            pval = 'P.VALUE'
        tbl[md]['rank'] = (np.arange(tbl[md].shape[0]) + 1.0) / gset_df.shape[0]  ### divide the total of the gene sets invert to rank percent
        md_rank[md] = pd.DataFrame({'GENE_SET': tbl[md]['GENE.SET'], (md+'_pval'):tbl[md][pval],(md+'_rank'):tbl[md]['rank'] })        
        rank_tb = rank_tb.merge(md_rank[md], left_on = 'GENE.SET', right_on = 'GENE_SET', how = 'outer')       
        rank_tb = rank_tb.drop('GENE_SET',axis=1)
    rank_arr = rank_tb.iloc[:,2:]
    rank_arr = rank_arr.fillna(1.0)

    #### classify ######
    comb_cls = svc.predict(rank_arr).astype(int)
    comb_prb = svc.predict_proba(rank_arr)   # 2 class, so 2-d array : probability to predict as 0 and 1
    comb_prb = np.max(comb_prb,axis=1)   #  the probability to predict correctly
    comb_dis = svc.decision_function(rank_arr)

    ### output the results ###
    rank_out = rank_tb.iloc[:,:2]    ### first 2 columns: gene set id, gene set description
    tmp = pd.DataFrame(np.column_stack([comb_dis,comb_cls,comb_prb]))
    tmp.columns = ['distance','class','prob']  ### 3-5 column: class,distance, probability ,
    rank_out['ENT'] = tmp.loc[:,'class']==1
    rank_out = pd.concat([rank_out, tmp.loc[:,['distance','prob']]], axis=1, ignore_index=True)

    rank_out.columns = ['GENE_SET','NAME','ENRICH_CLASS','R_SCORE','PROBABILITY']
    rank_out = rank_out.sort_values('R_SCORE', ascending=False)
    rank_out.to_csv( outdir +'combination_results.tsv',  sep="\t", header=True, index=False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "python "+ sys.argv[0] + "datadir  outdir "
    datadir = sys.argv[1]
    #save all 9 methods results
    gmtf = './data/kegg.pathway.hsa.320.gmt'
    # gmt file
    outdir = sys.argv[2]
    run_svm_18dims(datadir,gmtf,outdir)