import pdb
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import os
from multiprocessing import Pool,Manager
import numpy as np
from matplotlib import pyplot as plt

samplenum_THRD = 6
pvalue_THRD = 0.05
def extract_gene_pair(tissue,levelname):
    '''plot VENN and return the intersction tuple'''
    pathdir = r'/hwdata/home/wuyj2/eKPI/ArticleCancerStat/spearman_clean/'
    wpathdir = r'/hwdata/home/wuyj2/eKPI/ArticleCancerStat/VENN_spearman_Thrd%s/'%samplenum_THRD
    if not os.path.exists(wpathdir):
        os.makedirs(wpathdir)
    level_name_items = levelname.keys()
    mylength = len(level_name_items)
    level_name_items = sorted(level_name_items)
    level_sets = []
    set_labels = []
    set_colors_two = ('r','b')
    set_colors_three = ('r', 'b','g')
    fig,ax = plt.subplots()
    level_num = {}
    # pdb.set_trace()
    for level in level_name_items:
        file = levelname[level]
        tmp_set = set()
        with open(pathdir+file) as f:
            f.readline()
            if level=='Kinase.Phos':
                level = 'Phos'
                for line in f:
                    line = line.strip('\n').split('\t')
                    # print(line)
                    rnaname = line[0].split('_')[0]
                    genesite = line[1]
                    sample_num = int(line[4])
                    psdjustvalue = float(line[6])
                    if sample_num >= samplenum_THRD and psdjustvalue<pvalue_THRD:
                        pair = rnaname+'#'+genesite
                        tmp_set.add(pair)
                    else:
                        continue
            else:
                for line in f:
                    line = line.strip('\n').split('\t')
                    # print(line)
                    # pdb.set_trace()
                    rnaname = line[0]
                    genesite = line[1]
                    sample_num = int(line[4])
                    psdjustvalue = float(line[6])
                    if sample_num >= samplenum_THRD and psdjustvalue<pvalue_THRD:
                        pair = rnaname+'#'+genesite
                        tmp_set.add(pair)
                    else:
                        continue
        if len(tmp_set)==0:
            mylength -= 1
            continue
        level_sets.append(tmp_set)
        level_num[level] = len(tmp_set)
        set_labels.append(level)
        print(tissue,level,len(tmp_set))
    # pdb.set_trace()
    if mylength==2:
        v = venn2(subsets=level_sets,set_labels=tuple(set_labels),set_colors = set_colors_two,ax=ax,layout_algorithm=matplotlib_venn.layout.venn2.DefaultLayoutAlgorithm(fixed_subset_sizes=(1,1,1)))
        a = level_sets[0]
        b = level_sets[1]
        subsets = (len(a - b), len(b - a), len(a & b))
    elif mylength==3:
        v = venn3(subsets=level_sets, set_labels=tuple(set_labels), set_colors=set_colors_three,ax=ax,layout_algorithm=matplotlib_venn.layout.venn3.DefaultLayoutAlgorithm(fixed_subset_sizes=(1,1,1,1,1,1,1)))
        a = level_sets[0]
        b = level_sets[1]
        c = level_sets[2]
        subsets = (
            len(
                a - (b | c)
            ),  # TODO: This is certainly not the most efficient way to compute.
            len(b - (a | c)),
            len((a & b) - c),
            len(c - (a | b)),
            len((a & c) - b),
            len((b & c) - a),
            len(a & b & c),
        )
    # plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70, -70),
    #              ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #              arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='gray'))
    # v3.get_label_by_id('100').set_text(fontsize=15)
    elif len(level_name_items)==1:
        subsets = ()
        return tissue, level_name_items, subsets, level_num
    ax.set_title(tissue,fontsize=20)
    plt.savefig(wpathdir+tissue+'.pdf')
    # print(subsets)
    # pdb.set_trace()
    return tissue,level_name_items,subsets,level_num

#obtain the file_names of different level for each cancer
tissue_TN_level = {}
pathdir = r'/hwdata/home/wuyj2/eKPI/ArticleCancerStat/spearman_clean/'
for file in os.listdir(pathdir):
    if not file.endswith("txt"):
        continue
    fileN = file.split('_')
    tumorsize = fileN[0]
    TN = fileN[2]
    thelevel = fileN[3]
    name = tumorsize+'_'+TN
    # pdb.set_trace()
    # if tumorsize not in tumors:
    #     continue
    if name not in tissue_TN_level:
        tissue_TN_level[name] = {}
    tissue_TN_level[name][thelevel] = file

#obtain the intersection tuple for each cancer
tissue_TN_level_item =  tissue_TN_level.items()
pool = Pool(processes=40)
return_list = []
for tissue,level_name in tissue_TN_level_item:
    # df.loc[len(df)] = [tissue,level_name.get('Phos'),level_name.get('Pro'),level_name.get('RNA')]
    # if tissue=='GBM_Normal':
    # tissue,level_name_items,subsets,level_num = extract_gene_pair(tissue,level_name)
    # pdb.set_trace()
    rt = pool.apply_async(extract_gene_pair, args=(tissue,level_name,))
    return_list.append(rt.get())
pool.close()
pool.join()


#statistic the multiple ratios for each cancer
df = pd.DataFrame(columns=['Cancertype','Phos','Pro','Phos:Pro','RNA','Phos:RNA','Pro:RNA','Phos:Pro:RNA','Total_Phos','Total_Pro','Total_RNA','(Phos:Pro)/Pro','(Phos:Pro)/Phos','(Pro:RNA)/Pro','(Pro:RNA)/RNA'])
levels = ['Phos','Pro','RNA']
for tissue,levelname,values,level_num in return_list:
    myindex = len(df)
    if len(levelname)==3:
        tmp_value = [tissue]
        values = list(values)
        tmp_value.extend(list(values))
        tmp_value.extend([0,0,0,0,0,0,0])
        df.loc[myindex] = tmp_value
    elif len(levelname)==2:
        tmp_value = [tissue]
        values = list(values)
        tmp_value.extend(list(values))
        tmp_value.extend([0,0,0,0,0,0,0,0,0,0,0])
        df.loc[myindex] = tmp_value
        # for i in range(len(levelname)):
        #     df.loc[myindex,levelname[i]] = values[i]
        # df.loc[myindex,"%s:%s"%(levelname[0],levelname[1])] = values[-1]
    else:
        df.loc[myindex] = [tissue,0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0]
        # df.loc[myindex, levelname[0]] = values[0]
    for level in levels:
        if level_num.get(level):
            df.loc[myindex,'Total_%s'%level] = level_num.get(level)
            if level=='Pro':
                df.loc[myindex, '(Phos:Pro)/Pro'] = '{:.2%}'.format(df.loc[myindex, 'Phos:Pro']/level_num.get(level))
                df.loc[myindex, '(Pro:RNA)/Pro'] = '{:.2%}'.format(df.loc[myindex, 'Pro:RNA']/level_num.get(level))
            elif level=='RNA':
                df.loc[myindex, '(Pro:RNA)/RNA'] = '{:.2%}'.format(df.loc[myindex, 'Pro:RNA']/level_num.get(level))
            else:
                df.loc[myindex, '(Phos:Pro)/Phos'] = '{:.2%}'.format(df.loc[myindex, 'Phos:Pro']/level_num.get(level))
df.to_csv(r'/hwdata/home/wuyj2/eKPI/ArticleCancerStat/VENN_spearman_Thrd6/Cancertype_VENN_stat_cor_THRD%s.csv'%(samplenum_THRD))
