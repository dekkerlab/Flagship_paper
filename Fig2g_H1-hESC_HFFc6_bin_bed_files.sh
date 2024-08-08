bigWigToBedGraph H1-hESC_ChIA_PET_CTCF_4DNFIKVLJ93N.bw H1-hESC_ChIA_PET_CTCF_4DNFIKVLJ93N.bedGraph
bigWigToBedGraph H1-hESC_ChIA_PET_PolII_4DNFIAAHLWON.bw H1-hESC_ChIA_PET_PolII_4DNFIAAHLWON.bedGraph
bigWigToBedGraph H1-hESC_HiC_4DNFI4RWBB4U.bw H1-hESC_HiC_4DNFI4RWBB4U.bedGraph
bigWigToBedGraph H1-hESC_MicroC_4DNFI475YIT8.bw H1-hESC_MicroC_4DNFI475YIT8.bedGraph
bigWigToBedGraph H1-hESC_PLAC_Seq_4DNFINB2O3WR.bw H1-hESC_PLAC_Seq_4DNFINB2O3WR.bedGraph
bigWigToBedGraph HFFc6_ChIA_PET_CTCF_4DNFIKLR76CJ.bw HFFc6_ChIA_PET_CTCF_4DNFIKLR76CJ.bedGraph
bigWigToBedGraph HFFc6_ChIA_PET_PolII_4DNFIC2HQQIU.bw HFFc6_ChIA_PET_PolII_4DNFIC2HQQIU.bedGraph
bigWigToBedGraph HFFc6_Hi-C_4DNFINQZ5JHV.bw HFFc6_Hi-C_4DNFINQZ5JHV.bedGraph
bigWigToBedGraph HFFc6_MicroC_4DNFIFIPWIKW.bw HFFc6_MicroC_4DNFIFIPWIKW.bedGraph
bigWigToBedGraph HFFc6_PLAC_Seq_4DNFIW537ID7.bw HFFc6_PLAC_Seq_4DNFIW537ID7.bedGraph
bigWigToBedGraph HFFc6_SPRITE_4DNFIPBBLIPO.bw HFFc6_SPRITE_4DNFIPBBLIPO.bedGraph


fetchChromSizes hg38 | awk -v OFS="\t" '($1!~/_/){ print $1, "0", $2 }' | sort -k1,1 -k2,2 - > hg38.bed
binSize=50000
for i in $(ls | grep "bedGraph" | grep "E_L_counts");do
filename=$(echo $i | cut -f1 -d".");
bedtools sort -i $i | awk -v OFS="\t" '{ print $1, $2, $3, ".", $4 }' - > $filename.bed;
bedops --chop ${binSize} hg38.bed | bedmap --echo --max --prec 3 - $filename.bed > $filename.temp.bed
cat $filename.temp.bed | awk -F "|" '{print $1"\t"$2}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $filename.binned.50kb.bed;done


fetchChromSizes hg38 | awk -v OFS="\t" '($1!~/_/){ print $1, "0", $2 }' | sort -k1,1 -k2,2 - > hg38.bed
binSize=50000
for i in $(ls | grep "_2.bedGraph" | grep "ESC" );do
filename=$(echo $i | cut -f1 -d".");
bedtools sort -i $i | awk -v OFS="\t" '{ print $1, $2, $3, ".", $4 }' - > $filename.bed;
bedops --chop ${binSize} hg38.bed | bedmap --echo --max --prec 3 - $filename.bed > $filename.temp.bed
cat $filename.temp.bed | awk -F "|" '{print $1"\t"$2}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4}' > $filename.binned.50kb.bed;done


fetchChromSizes hg38 | awk -v OFS="\t" '($1!~/_/){ print $1, "0", $2 }' | sort -k1,1 -k2,2 - > hg38.bed
binSize=50000
for i in $(ls | grep "_2.bedGraph" | grep "HFFc6" );do
filename=$(echo $i | cut -f1 -d".");
bedtools sort -i $i | awk -v OFS="\t" '{ print $1, $2, $3, ".", $4 }' - > $filename.bed;
bedops --chop ${binSize} hg38.bed | bedmap --echo --max --prec 3 - $filename.bed > $filename.temp.bed
cat $filename.temp.bed | awk -F "|" '{print $1"\t"$2}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4}' > $filename.binned.50kb.bed;done
