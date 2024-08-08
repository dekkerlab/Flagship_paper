
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.style.use('seaborn-white')
import multiprocess as mp
import numpy as np
import pandas as pd
import cooltools
import cooler
#import bbi
import bioframe

conditions  = [
"H1-hESC_MicroC",
"H1-hESC_FA_DSG_DpnII_Hi-C",
"H1-hESC_PLAC-Seq",
"H1-hESC_ChIA_PET_PolII",
"H1-hESC_ChIA_PET_CTCF",
"H1-hESC_DNA_clusters_R1_R2_2-100",
"H1-hESC-GAM"]



path="mcool/"
out="insulation/"


path2="insulation/"

ins_files=[
"4DNFI9GMP2J8_H1-hESC_MicroC.25kb.cis.bed",
"4DNFI82R42AD_H1-hESC_FA_DSG_DpnII_Hi-C.25kb.cis.bed",
"4DNFICOGAKW2_H1-hESC_PLAC-Seq.25kb.cis.bed",
"4DNFIF1J6GC5_H1-hESC_ChIA_PET_PolII.25kb.cis.bed",
"4DNFIK3276U7_H1-hESC_ChIA_PET_CTCF.25kb.cis.bed",
"4DNFIV3PDS5F_4DNFIIY1TXUZ_H1-hESC_DNA_clusters_R1_R2_2-100.25kb.cis.bed",
"GAM_H1_new.IS.scores.at25Kb_formatted.bed"]


names = [
"MicroC",
"Hi-C",
"PLAC-Seq",
"ChIA-PET PolII",
"ChIA-PET CTCF",
"SPRITE",
"GAM"]


filename1="4DNFI9GMP2J8_H1-hESC_MicroC.mcool"
filename2="4DNFI82R42AD_H1-hESC_FA_DSG_DpnII_Hi-C.mcool"
filename3="4DNFICOGAKW2_H1-hESC_PLAC-Seq.mcool"
filename4="4DNFIF1J6GC5_H1-hESC_ChIA_PET_PolII.mcool"
filename5="4DNFIK3276U7_H1-hESC_ChIA_PET_CTCF.mcool"
filename6="H1-hESC_DNA_clusters_R1_R2_2-100.mcool"
filename7="H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at25Kb.mcool"




cooler_paths = {

    
    "H1-hESC_MicroC" : path+filename1+'::/resolutions/25000',
    "H1-hESC_FA_DSG_DpnII_Hi-C": path+filename2+'::/resolutions/25000',
    "H1-hESC_PLAC-Seq": path+filename3+'::/resolutions/25000',
    "H1-hESC_ChIA_PET_PolII": path+filename4+'::/resolutions/25000',
    "H1-hESC_ChIA_PET_CTCF": path+filename5+'::/resolutions/25000',
    "H1-hESC_DNA_clusters_R1_R2_2-100": path+filename6+'::/resolutions/25000',
    "H1-hESC-GAM": 'GAM_2024/'+filename7+'::/resolutions/25000'
    
}


long_names = {
    
    "H1-hESC_MicroC" : "H1-hESC_MicroC",
    "H1-hESC_FA_DSG_DpnII_Hi-C" : "H1-hESC_FA_DSG_DpnII_Hi-C",
    "H1-hESC_PLAC-Seq":  "H1-hESC_PLAC-Seq",
    "H1-hESC_ChIA_PET_PolII": "H1-hESC_ChIA_PET_PolII",
    "H1-hESC_ChIA_PET_CTCF": "H1-hESC_ChIA_PET_CTCF",
    "H1-hESC_DNA_clusters_R1_R2_2-100" : "H1-hESC_DNA_clusters_R1_R2_2-100",
    "H1-hESC-GAM" : "H1-hESC-GAM"
    
}

#exp_name_list=[file.split(".")[0]+"25000.cis.expected" for file in cool]   




#ins_name_list=[file.split(".")[0]+".txt" for file in cool]   

path_gam="GAM_2024/"
insulation_path={
    "H1-hESC_MicroC" : path2+ins_files[0],
    "H1-hESC_FA_DSG_DpnII_Hi-C": path2+ins_files[1],
    "H1-hESC_PLAC-Seq": path2+ins_files[2],
    "H1-hESC_ChIA_PET_PolII": path2+ins_files[3],
    "H1-hESC_ChIA_PET_CTCF": path2+ins_files[4],
    "H1-hESC_DNA_clusters_R1_R2_2-100": path2+ins_files[5],
    "H1-hESC-GAM" : path2+ins_files[6] 
}


all_m=[]
new_data=pd.DataFrame()

index_list_dict={}
for i,cond in enumerate(conditions):
    insulations_new_all=pd.DataFrame()
    insulations = pd.read_table(insulation_path[cond])
    print(insulations.head(5))
    insulations = insulations.dropna()
    if insulations.shape[1]>4:
        insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    else:
        insulations.columns=["chrom","start","end","log2_insulation_score_100000"]

    strong_weak=insulations[insulations["log2_insulation_score_100000"] !=0]
    print(strong_weak.shape)
    m=np.mean(strong_weak["log2_insulation_score_100000"])
    all_m.append(m)
    #name=insulation_path[cond].split("/")[-1].split(".")[0]
    #print(name)
    x = np.log10(strong_weak['log2_insulation_score_100000'].values)
    bins = np.linspace(x.min(), x.max(), num=100)
    boundary_list=strong_weak[strong_weak["log2_insulation_score_100000"] <=  -0.226]
    if insulations.shape[1]>4:
        insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    else:
        insulations.columns=["chrom","start","end","log2_insulation_score_100000"]
    #boundary_list.to_csv("ESC_"+cond+"_strong_boundaries.bed",sep="\t",index=None)
    for chrom in np.unique(insulations[["chrom"]]):
        insulations_new_chr=pd.DataFrame()
        boundary_list_chr=boundary_list[boundary_list["chrom"]==chrom]
        mid_point_new=boundary_list_chr[["start","end"]].mean(axis=1)+1
        insulations_chr=insulations[insulations["chrom"]==chrom]
        for i in mid_point_new:
            temp=insulations_chr[insulations_chr.eval('(start < {}) & ({} < end)'.format(i,i))]
            insulations_new_chr=pd.concat([insulations_new_chr,pd.DataFrame(temp)],ignore_index=False)
        insulations_new_all=pd.concat([insulations_new_all,insulations_new_chr])
    print(insulations_new_all.shape)
    index_list_dict[cond]=insulations_new_all.index


def pile_up(df,index_list):
    start = df.iloc[index_list].index-10
    end  = df.iloc[index_list].index+11
    #print(end)
    return (np.nanmean(np.array(list(map(lambda x: df["log2_insulation_score_100000"].values[x[0]:x[1]],  zip(start,end)))), axis=0))



i=0
ins_data=pd.DataFrame()
file_list=[]
for file_name in ins_files:
    print(file_name)
    file_list.append(file_name.split('.')[0])
    file=pd.read_csv(path2+file_name,sep="\t")
    #print(file.head(5))
    ins=file["log2_insulation_score_100000"]
    ins_data=pd.concat([ins_data,ins],axis=1)

    



hg38 = bioframe.fetch_chromsizes('hg38')
chromsizes = bioframe.fetch_chromsizes('hg38')
chromosomes = list(chromsizes.index)

print(filename1,filename2,filename3)


path="mcool/"
out="insulation/"



hg38 = bioframe.fetch_chromsizes('hg38')
chromsizes = bioframe.fetch_chromsizes('hg38')
chromosomes = list(chromsizes.index)

print(filename1,filename2,filename3,filename4,filename5,filename6,filename7)


cool = [filename1,filename2,filename3,filename4,filename5,filename6,filename7]



clrs = {
    cond: cooler.Cooler(cooler_paths[cond]) for cond in conditions
}


gs = GridSpec(nrows=1, ncols=1)
plt.figure(figsize=(8,8))

k=0


names=conditions



names = [
"MicroC",
"Hi-C",
"PLAC-Seq",
"ChIA-PET PolII",
"ChIA-PET CTCF",
"SPRITE",
"GAM"]

c=['g','dodgerblue','orange','cyan','pink','r','k']


#adds=[0.04, 0, 0.08,0,0,0.12,0.13 ]

#adds=[0.04, 0, 0.08,0.07,0,-0.12,-0.03 ]

for i, cond in enumerate(conditions):
    insulations = pd.read_table(insulation_path[cond])
    #insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    if insulations.shape[1]>4:
        insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    else:
        insulations.columns=["chrom","start","end","log2_insulation_score_100000"]
    ax = plt.subplot(gs[k])
    print(cond)
    #print(np.max(pile_up(insulations,index_list_dict[cond].values)))
    sub=pile_up(insulations,index_list_dict[cond].values)[0]

    #print(np.min(pile_up(insulations,index_list_dict[cond].values)))
    ins=pile_up(insulations,index_list_dict[cond].values)
    ins_new=ins-sub
    img=ax.plot(ins_new,label=names[i],color=c[i])
    ax.xaxis.tick_bottom()
    ax.set_ylim(-0.65,0.4)
ax.legend()
#plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.ylabel("log2 insulation score")


plt.tight_layout()
plt.savefig("Figure_2e_Hi_hESC_insulation_pileup_25kb.pdf")


##### ##################



path="dcic_mcool/mcool/"
out="insulation/"


path2="insulation/"


ins_files=[
"4DNFI9FVHJZQ_HFFc6_MicroC.25kb.cis.bed",
"4DNFIAVXXO55_HFFc6_FA_DSG_DpnII_Hi-C.25kb.cis.bed",
"4DNFI9REIU8H_HFFc6_PLAC-Seq.25kb.cis.bed",
"4DNFIOOSHTQV_HFFc6_ChIA_PET_PolII.25kb.cis.bed",
"4DNFIELOAD41_HFFc6_ChIA_PET_CTCF.25kb.cis.bed",
"4DNFIRXON7Z2_HFFc6_DNA_clusters_R1_R2_2-100.25kb.cis.bed"]

names = [
"MicroC",
"Hi-C",
"PLAC-Seq",
"ChIA-PET PolII",
"ChIA-PET CTCF",
"SPRITE"]

filename1="4DNFI9FVHJZQ_HFFc6_MicroC.mcool"
filename2="4DNFI9REIU8H_HFFc6_PLAC-Seq.mcool"
filename3="4DNFIAVXXO55_HFFc6_FA_DSG_DpnII_Hi-C.mcool"
filename4="4DNFIELOAD41_HFFc6_ChIA_PET_CTCF.mcool"
filename5="4DNFIOOSHTQV_HFFc6_ChIA_PET_PolII.mcool"
filename6="4DNFIRXON7Z2_HFFc6_DNA_clusters_R1_R2_2-100.mcool"


path="mcool/"
out="insulation/"


cooler_paths = {

    
    "HFFc6_MicroC" : path+filename1+'::/resolutions/25000',
    "HFFc6_FA_DSG_DpnII_Hi-C": path+filename2+'::/resolutions/25000',
    "HFFc6_PLAC-Seq": path+filename3+'::/resolutions/25000',
    "HFFc6_ChIA_PET_PolII": path+filename4+'::/resolutions/25000',
    "HFFc6_ChIA_PET_CTCF": path+filename5+'::/resolutions/25000',
    "HFFc6_DNA_clusters_R1_R2_2-100": path+filename6+'::/resolutions/25000'    
}


long_names = {
    
    "HFFc6_MicroC" : "HFFc6_MicroC",
    "HFFc6_FA_DSG_DpnII_Hi-C" : "HFFc6_FA_DSG_DpnII_Hi-C",
    "HFFc6_PLAC-Seq":  "HFFc6_PLAC-Seq",
    "HFFc6_ChIA_PET_PolII": "HFFc6_ChIA_PET_PolII",
    "HFFc6_ChIA_PET_CTCF": "HFFc6_ChIA_PET_CTCF",
    "HFFc6_DNA_clusters_R1_R2_2-100" : "HFFc6_DNA_clusters_R1_R2_2-100",
    
}

#exp_name_list=[file.split(".")[0]+"25000.cis.expected" for file in cool]   

conditions  = [
"HFFc6_MicroC",
"HFFc6_FA_DSG_DpnII_Hi-C",
"HFFc6_PLAC-Seq",
"HFFc6_ChIA_PET_PolII",
"HFFc6_ChIA_PET_CTCF",
"HFFc6_DNA_clusters_R1_R2_2-100"]


#ins_name_list=[file.split(".")[0]+".txt" for file in cool]   

path_gam="GAM_2024/"
insulation_path={
    "HFFc6_MicroC" : path2+ins_files[0],
    "HFFc6_FA_DSG_DpnII_Hi-C": path2+ins_files[1],
    "HFFc6_PLAC-Seq": path2+ins_files[2],
    "HFFc6_ChIA_PET_PolII": path2+ins_files[3],
    "HFFc6_ChIA_PET_CTCF": path2+ins_files[4],
    "HFFc6_DNA_clusters_R1_R2_2-100": path2+ins_files[5]
}


all_m=[]
new_data=pd.DataFrame()

index_list_dict={}
for i,cond in enumerate(conditions):
    insulations_new_all=pd.DataFrame()
    insulations = pd.read_table(insulation_path[cond])
    print(insulations.head(5))
    insulations = insulations.dropna()
    if insulations.shape[1]>4:
        insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    else:
        insulations.columns=["chrom","start","end","log2_insulation_score_100000"]

    strong_weak=insulations[insulations["log2_insulation_score_100000"] !=0]
    print(strong_weak.shape)
    m=np.mean(strong_weak["log2_insulation_score_100000"])
    all_m.append(m)
    #name=insulation_path[cond].split("/")[-1].split(".")[0]
    #print(name)
    x = np.log10(strong_weak['log2_insulation_score_100000'].values)
    bins = np.linspace(x.min(), x.max(), num=100)
    boundary_list=strong_weak[strong_weak["log2_insulation_score_100000"] <=  -0.25]
    if insulations.shape[1]>4:
        insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    else:
        insulations.columns=["chrom","start","end","log2_insulation_score_100000"]
    #boundary_list.to_csv("ESC_"+cond+"_strong_boundaries.bed",sep="\t",index=None)
    for chrom in np.unique(insulations[["chrom"]]):
        insulations_new_chr=pd.DataFrame()
        boundary_list_chr=boundary_list[boundary_list["chrom"]==chrom]
        mid_point_new=boundary_list_chr[["start","end"]].mean(axis=1)+1
        insulations_chr=insulations[insulations["chrom"]==chrom]
        for i in mid_point_new:
            temp=insulations_chr[insulations_chr.eval('(start < {}) & ({} < end)'.format(i,i))]
            insulations_new_chr=pd.concat([insulations_new_chr,pd.DataFrame(temp)],ignore_index=False)
        insulations_new_all=pd.concat([insulations_new_all,insulations_new_chr])
    print(insulations_new_all.shape)
    index_list_dict[cond]=insulations_new_all.index


def pile_up(df,index_list):
    start = df.iloc[index_list].index-10
    end  = df.iloc[index_list].index+11
    #print(end)
    return (np.nanmean(np.array(list(map(lambda x: df["log2_insulation_score_100000"].values[x[0]:x[1]],  zip(start,end)))), axis=0))


c=['g','dodgerblue','orange','cyan','pink','r']

names = [
"MicroC",
"Hi-C",
"PLAC-Seq",
"ChIA-PET PolII",
"ChIA-PET CTCF",
"SPRITE"]
i=0
ins_data=pd.DataFrame()
file_list=[]
for file_name in ins_files:
    print(file_name)
    file_list.append(file_name.split('.')[0])
    file=pd.read_csv(path2+file_name,sep="\t")
    #print(file.head(5))
    ins=file["log2_insulation_score_100000"]
    ins_data=pd.concat([ins_data,ins],axis=1)

    



hg38 = bioframe.fetch_chromsizes('hg38')
chromsizes = bioframe.fetch_chromsizes('hg38')
chromosomes = list(chromsizes.index)

print(filename1,filename2,filename3)



hg38 = bioframe.fetch_chromsizes('hg38')
chromsizes = bioframe.fetch_chromsizes('hg38')
chromosomes = list(chromsizes.index)

print(filename1,filename2,filename3,filename4,filename5,filename6,filename7)


cool = [filename1,filename2,filename3,filename4,filename5,filename6,filename7]



clrs = {
    cond: cooler.Cooler(cooler_paths[cond]) for cond in conditions
}


gs = GridSpec(nrows=1, ncols=1)
plt.figure(figsize=(8,8))
k=0





#adds=[0.04, 0, 0.08,0,0,0.12,0.13 ]

#adds=[0.04, 0, 0.08,0.07,0,-0.12,-0.03 ]

for i, cond in enumerate(conditions):
    insulations = pd.read_table(insulation_path[cond])
    #insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    if insulations.shape[1]>4:
        insulations.columns=["chrom","start","end","is_bad_bin","log2_insulation_score_100000","n_valid_pixels_100000","boundary_strength_100000"]
    else:
        insulations.columns=["chrom","start","end","log2_insulation_score_100000"]
    ax = plt.subplot(gs[k])
    print(cond)
    #print(np.max(pile_up(insulations,index_list_dict[cond].values)))
    sub=pile_up(insulations,index_list_dict[cond].values)[0]

    #print(np.min(pile_up(insulations,index_list_dict[cond].values)))
    ins=pile_up(insulations,index_list_dict[cond].values)
    ins_new=ins-sub
    img=ax.plot(ins_new,label=names[i],color=c[i])
    ax.xaxis.tick_bottom()
    ax.set_ylim(-0.65,0.4)
ax.legend()
#plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.ylabel("log2 insulation score")


plt.tight_layout()
plt.savefig("Figure_2e_HFFc6_insulation_pileup_25kb.pdf")
