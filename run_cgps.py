import sys, os
import subprocess
from optparse import OptionParser
import numpy as np
import pandas as pd

def get_dir():
    srcf = os.environ['HOME']+'/.cgpsrc'
    if not os.path.isfile(srcf):
        sys.exit('please add .cgpsrc under HOME folder\n')
    with open(srcf,'r') as f:
        tmp = f.readline()
    cgps_home = tmp.split('=')[1].strip()
    return cgps_home

def config_option():
    usage = 'Usage: run_cgps.py -e expfile -p phefile -d datatype -s species -o outdir'
    p = OptionParser(usage)
    p.add_option(
        '-e','--expfile',dest='expfile',action='store',
        help='')
    p.add_option(
        '-p','--phefile',dest='phefile',action='store',
        help='')
    p.add_option(
        '-d','--datatype',dest='datatype',action='store',
        help='')
    p.add_option(
        '-s','--species',dest='spe',action='store',
        help='')
    p.add_option(
        '-o','--outdir',dest='outdir',action='store',
        help='absolute path')
    p.add_option('-g','--gmtf',dest='gmtf',action='store',
                 help='absolute path')
    opt, args = p.parse_args()
    return p,opt,args

#################################
##### run individual methods ####
#################################
cgps_home = get_dir()
opt_parser, opt, args = config_option()
print opt.expfile,opt.phefile,opt.datatype,opt.outdir
if len(sys.argv) == 1:
    opt_parser.print_help()
    sys.exit(1)

cmdline = 'Rscript '+cgps_home+'/scripts/combined_methods.R '+ ' '.join([cgps_home,opt.expfile, opt.phefile, opt.datatype, opt.spe, opt.outdir,opt.gmtf])
print cmdline
subprocess.call(cmdline,shell=True)

################################
###### Predict by CGPS #########
################################

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

def run_svm_18dims(datadir,gmtf,outdir,svmfile):
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
    #comb_dis = svc.decision_function(rank_arr)
    comb_dis = - np.log( 1.0 - comb_prb  )

    ### output the results ###
    rank_out = rank_tb.iloc[:,:2]    ### first 2 columns: gene set id, gene set description
    tmp = pd.DataFrame(np.column_stack([comb_dis,comb_cls,comb_prb]))
    tmp.columns = ['distance','class','prob']  ### 3-5 column: class,distance, probability ,
    rank_out['ENT'] = tmp.loc[:,'class']==1
    rank_out = pd.concat([rank_out, tmp.loc[:,['distance','prob']]], axis=1, ignore_index=True)

    rank_out.columns = ['GENE_SET','NAME','ENRICH_CLASS','R_SCORE','PROBABILITY']
    rank_out = rank_out.sort_values('R_SCORE', ascending=False)
    rank_out = rank_out.loc[:,['GENE_SET','NAME','R_SCORE']]
    rank_out.to_csv( outdir +'combination_results.tsv',  sep="\t", header=True, index=False)

gmtf = opt.gmtf
#gmtf = cgps_home + '/data/kegg.'+opt.spe+'.gmt'
svmfile = cgps_home + '/data/cgps_model.pkl'
#outdir = cgps_home + '/test/res/'
run_svm_18dims(opt.outdir+'/',gmtf,opt.outdir+'/',svmfile)

