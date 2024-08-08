
################ Testing the new saddle plot #############


input_path=/Users/betulakgoloksuz/Desktop/gam_new
bed_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD 
comp_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD
out_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD/npz_July2024_final_normalized_mod/


comp_file=H1-hESC_DamID-seq_normalized_counts_4DNFIXNBG8L1_2.binned.bed

for i in $(ls $input_path | grep "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mcool") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
filename_final=$(echo "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mod.cis.expected"
echo $input
echo $comp_file
echo $expected

cooltools compute-saddle --strength --vmin 0.1 --vmax 0.4 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -2.24 3.93 --contact-type cis --out-prefix $out_path/$filename'_GAM_ESC' --fig pdf $input $comp_path/$comp_file $expected;

done

##

comp_file=H1-hESC_E_L_counts_4DNFI4BSJRMF_4DNFIISI1ZA8.binned.50kb.bed


for i in $(ls $input_path | grep "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mcool") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
filename_final=$(echo "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mod.cis.expected"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength ---vmin 0.1 --vmax 0.4  --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 28.30 --contact-type cis --out-prefix $out_path/$filename'_GAM_ESC' --fig pdf $input $comp_path/$comp_file $expected;
done

###


##

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFI625PP2A_2_binned.bed


for i in $(ls $input_path | grep "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mcool") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
filename_final=$(echo "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mod.cis.expected"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.1 --vmax 0.4  --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -4.76 6.94 --contact-type cis --out-prefix $out_path/$filename'_GAM_ESC' --fig pdf $input $comp_path/$comp_file $expected;
done

#####


comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFICLC27JC_2_binned.bed


for i in $(ls $input_path | grep "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mcool") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
filename_final=$(echo "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mod.cis.expected"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.1 --vmax 0.4  --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -3.42 5.8 --contact-type cis --out-prefix $out_path/$filename'_GAM_ESC' --fig pdf $input $comp_path/$comp_file $expected;
done

###

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFIFNIK4HD_2_binned.bed

for i in $(ls $input_path | grep "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mcool") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
filename_final=$(echo "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mod.cis.expected"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.1 --vmax 0.4  --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -4.12 8.07 --contact-type cis --out-prefix $out_path/$filename'_GAM_ESC' --fig pdf $input $comp_path/$comp_file $expected;
done

####

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFILHYV7Q7_2_binned.bed


for i in $(ls $input_path | grep "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mcool") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
filename_final=$(echo "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mod.cis.expected"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.1 --vmax 0.4  --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -3.76 4.41 --contact-type cis --out-prefix $out_path/$filename'_GAM_ESC' --fig pdf $input $comp_path/$comp_file $expected;
done

##

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFILYDJU8T_2_binned.bed


for i in $(ls $input_path | grep "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mcool") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
filename_final=$(echo "H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/H1hESC.GAM.NPMI.normalized.pairwise.curated.matrices.at50Kb.mod.cis.expected"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.1 --vmax 0.4  --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -10.29 3.98 --contact-type cis --out-prefix $out_path/$filename'_GAM_ESC' --fig pdf $input $comp_path/$comp_file $expected;
done


##################################


input_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD/h1ESC_mcool
bed_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD 
comp_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD
out_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD/npz_July2024_final_normalized_mod/


comp_file=H1-hESC_DamID-seq_normalized_counts_4DNFIXNBG8L1_2.binned.bed

for i in $(ls $input_path | grep ".mcool$"| grep "H1-hESC") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $bed_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -2.24 3.93 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

##

comp_file=H1-hESC_E_L_counts_4DNFI4BSJRMF_4DNFIISI1ZA8.binned.50kb.bed


for i in $(ls $input_path | grep ".mcool$"| grep "H1-hESC") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $bed_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 28.30 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

###


##

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFI625PP2A_2_binned.bed


for i in $(ls $input_path | grep ".mcool$"| grep "H1-hESC") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $bed_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -4.76 6.94 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

#####


comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFICLC27JC_2_binned.bed


for i in $(ls $input_path | grep ".mcool$"| grep "H1-hESC") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $bed_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -3.42 5.8 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

###

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFIFNIK4HD_2_binned.bed

for i in $(ls $input_path | grep ".mcool$"| grep "H1-hESC") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $bed_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -4.12 8.07 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

####

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFILHYV7Q7_2_binned.bed


for i in $(ls $input_path | grep ".mcool$"| grep "H1-hESC") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $bed_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -3.76 4.41 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

##

comp_file=H1-hESC_TSA-seq_normalized_counts_4DNFILYDJU8T_2_binned.bed

for i in $(ls $input_path | grep ".mcool$"| grep "H1-hESC") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $bed_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range -10.29 3.98 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done


####################################################################



############  HFFc6   ####################################



##################################

input_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD
bed_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD 
comp_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD/npz_July2024_final/
out_path=/Users/betulakgoloksuz/Desktop/gam_new/oneD/npz_July2024_final_normalized_mod/


comp_file=HFFc6_DamID-seq_normalized_counts_4DNFI7724Y7Q_2_50kb.bed

for i in $(ls $input_path | grep ".mcool$"| grep "HFFc6") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 3.16 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

##

comp_file=HFFc6_E_L_counts_4DNFI7EOU166_4DNFI9FTS684_50kb.bed


for i in $(ls $input_path | grep ".mcool$"| grep "HFFc6") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 37.37 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

###


##

comp_file=HFFc6_TSA-seq_normalized_counts_4DNFI6FTPH5V_2_50kb.bed


for i in $(ls $input_path | grep ".mcool$"| grep "HFFc6") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 2.89 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

#####


comp_file=HFFc6_TSA-seq_normalized_counts_4DNFI8CJWFUT_2_50kb.bed


for i in $(ls $input_path | grep ".mcool$"| grep "HFFc6") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 4.71 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

###

comp_file=HFFc6_TSA-seq_normalized_counts_4DNFIEWPL92Y_2_50kb.bed

for i in $(ls $input_path | grep ".mcool$"| grep "HFFc6") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 2.61 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

####

comp_file=HFFc6_TSA-seq_normalized_counts_4DNFIHYF6H13_2_50kb.bed


for i in $(ls $input_path | grep ".mcool$"| grep "HFFc6") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 2.94 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done

##

comp_file=HFFc6_TSA-seq_normalized_counts_4DNFIMTKDNJW_2_50kb.bed

for i in $(ls $input_path | grep ".mcool$"| grep "HFFc6") ; do \
filename=$(echo $comp_file | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i::resolutions/50000
cool_name=$(echo $i | cut -f1 -d".");
#compartment=$comp_path/$filename".vecs.tsv"
expected=$input_path"/$cool_name.50000.cis.expected.tsv"
echo $input
echo $comp_file
echo $expected
cooltools compute-saddle --strength --vmin 0.5 --vmax 2 --regions $input_path/hg38_chromsizes.bed --qrange 0.02 0.98 --range 0 2.52 --contact-type cis --out-prefix $out_path/$filename'_'$cool_name --fig pdf $input $comp_path/$comp_file $expected;
done