#!/bin/bash
###################INSTRUCTIONS########################

#YOU WILL NEED 6 FILES in order to run this program: (1) reference cds sequence file (NCBI): input= "$name"2.txt 
#(2) The converted genomic positions to cds positions file (from UCSC site): input= "$name"_converted_data.txt 
#(3) Rscript merge file: input= Rscript.R 
#(4) Rscript merge file for odd columns: input= Rscript_o.R 
#(5) Rscript merge file for even columns: input= Rscript_e.R 
#(6) "$name" vcf file 
#be sure to get transcript name (transcript="insert name here with quotes")
#be sure of the orientation of the bases (for example, do you need to change nucleotide base to its complement?)
#be sure to change the paths and lists in this program BUT make sure quotations are vertical as opposed to being at an angle

#######################################################

echo "--------------------BEGIN Here By Defining Variables -------------------------"
read -p "Enter the name of your gene: " name
read -p "What is the name of your transcript? " transcript
read -p "Are sure this transcript name is correct? [y/n]" response

if [[ $response = "n" ]];
then
	exit
elif [[ $response = "y" ]]; 
then
	echo "Your gene name is $name!" 
	echo "The name of your transcript is $transcript!"
	echo "Proceed young Jedi"
fi 
####CHANGE RScript sed command below#######################
mypath=/home/wgriffith/scripts
until [ -f "$mypath"/"$name"/"$name"_MACPRF.fas ]
do

echo "--------------------Create R scripts Based on template ----------------"
touch "$mypath"/RScript_"$name".R
touch "$mypath"/RScript_"$name"_o.R
touch "$mypath"/RScript_"$name"_e.R

cp "$mypath"/RScript_CDH2.R "$mypath"/RScript_"$name".R
cp "$mypath"/RScript_CDH2_o.R "$mypath"/RScript_"$name"_o.R
cp "$mypath"/RScript_CDH2_e.R "$mypath"/RScript_"$name"_e.R

var="CDH2"
sed -i "s/$var/$name/g" "$mypath"/RScript_"$name".R #replace name CDH2 with new name
sed -i "s/$var/$name/g" "$mypath"/RScript_"$name"_o.R #replace name CDH2 with new name
sed -i "s/$var/$name/g" "$mypath"/RScript_"$name"_e.R #replace name CDH2 with new name

##NO MORE CHANGES BEYOND THIS POINT#########################
echo "------------------- Format mRNA reference sequence for Homo from UCSC -----------------------"

sed 's/[0-9]//g' <"$mypath"/"$name"2.txt | tr '\n' ' ' | sed 's/ //g' | awk '{ print toupper($0) }'| sed 's/./& /g' | tr ' ' '\n' | grep 'A\|T\|G\|C' > "$mypath"/"$name"3.txt
sed = "$mypath"/"$name"3.txt | sed 'N;s/\n/ /g' > "$mypath"/"$name"5.txt 
#sed '$d' "$mypath"/"$name"4.txt > "$mypath"/"$name"5.txt

#sed 's/[0-9]//g' <"$mypath"/"$name"2.txt | tr '\n' ' ' | sed 's/ //g' | awk '{ print toupper($0) }'| sed 's/./& /g' | tr ' ' '\n' > "$mypath"/"$name"3.txt
#sed = "$mypath"/"$name"3.txt | sed 'N;s/\n/ /g' > "$mypath"/"$name"4.txt 
#sed '$d' "$mypath"/"$name"4.txt > "$mypath"/"$name"5.txt

echo "------------------Isolate correct transcript from variant annotator at UCSC ($8 or $9 contain data) -------------------"

grep -w "$transcript" <"$mypath"/"$name"_converted_data.txt | awk '{ print $1 " " $2 " " $3 " " $9 }' | awk '{ print $4 " " $2 " " $3 " " }' | sed 's/ch.*://g'> "$mypath"/"$name"_converted_data1.txt #arrange the columns

if [ -s "$mypath"/"$name"_converted_data1.txt ]; then
echo "Moving on to R";
else 
echo "Your file is empty, girl!"
exit;
fi; 

#############################
#############################
#Merge file1 and file2 in R #
#CONTROL D to exit from R   #
#############################
#############################
#/ycga-ba/ba_sequencers1/scratch/mcc77/PRM1_merged_data2.txt
#R CMD BATCH Rscript.R

#mypath=/ycga-ba/ba_sequencers1/scratch/mcc77

#name[1]="ASPM"
#transcript="uc001gtu.3"

#for (( i=1; i<=1; i++ ))
#do
R CMD BATCH /home/wgriffith/scripts/RScript_"$name".R

echo "--------------------Clean up R output file ------------------------"

awk '{ print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' <"$mypath"/"$name"_merged_data.txt | sed 's/"//g;s/_//g;s/ch.*://g' | grep -v 'REJECT\|V2.x\|V2.y\|V3\|V4\|-' > "$mypath"/"$name"_merged_data1.txt

echo "---------------------Rearrange the columns in the previously merged file---------------------------"

awk '{ print $1 "\t" $2 "\t" $4 "\t" $3 }' "$mypath"/"$name"_merged_data1.txt > "$mypath"/"$name"_merged_data2.txt #now chromosomal position is in the same column as the impute.legend file

echo "----------------- Convert to phased genotypes --------------------------------"

/usr/java/jdk1.8.0_65/bin/java -jar "$mypath"/beagle.07Jan16.9b4.jar gt="$mypath"/1000Genomes_P3_"$name".recode.vcf out="$mypath"/"$name"_phased

#ref="$mypath"/ALL."$chr".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz 

echo "----------------- Convert vcf of interest to an IMPUTE file ----------------------------" 

"$mypath"/vcftools_0.1.13/bin/./vcftools --gzvcf "$mypath"/"$name"_phased.vcf.gz --IMPUTE --phased --out "$mypath"/"$name"

grep -v 'REJECT\|ID pos allele0 allele1' "$mypath"/"$name".impute.legend | awk '{ print $2 " " $3 " " $4 }' > "$mypath"/"$name".impute1.legend
 
sed 's/ //g;s/^ //g' "$mypath"/"$name".impute.hap > "$mypath"/"$name".impute1.hap

paste "$mypath"/"$name".impute1.legend "$mypath"/"$name".impute1.hap > "$mypath"/"$name".impute2.hap

awk '{ gsub("0",$2,$4); print $0}' "$mypath"/"$name".impute2.hap > "$mypath"/"$name".impute3.hap

awk '{ gsub("1",$3,$4); print $0}' "$mypath"/"$name".impute3.hap > "$mypath"/"$name".impute4.hap

awk '{ print $4 }' "$mypath"/"$name".impute4.hap | sed 's/./& /g;s/^ //g' > "$mypath"/"$name".impute5.hap

echo "----------------- Print odd columns ------------------------"

awk  '{for (i=1;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' "$mypath"/"$name".impute5.hap > "$mypath"/"$name".phased_odd.hap

echo "----------------- Print even columns ------------------------"

awk  '{for (i=2;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' "$mypath"/"$name".impute5.hap > "$mypath"/"$name".phased_even.hap

echo "---------------- Paste the map file to the haplotype file --------------------------" 

awk '{ print $1 " " $1 " " $1 " " $1 }' "$mypath"/"$name".impute1.legend > "$mypath"/"$name".phased_ids.tped

paste "$mypath"/"$name".phased_ids.tped "$mypath"/"$name".phased_odd.hap > "$mypath"/"$name"_step1.txt #odd file

paste "$mypath"/"$name".phased_ids.tped "$mypath"/"$name".phased_even.hap > "$mypath"/"$name"_step2.txt #even file


#############################
#############################
#Paste file1 and file2 in R #
#############################
#############################

R CMD BATCH /home/wgriffith/scripts/RScript_"$name"_o.R
R CMD BATCH /home/wgriffith/scripts/RScript_"$name"_e.R

cd /home/wgriffith/scripts

	read -r -p "Do you need the reverse complementary sequence? [y/n]" response
if [[ $response = "y" ]];
then
	echo "--------------------Clean up R output files (odd and even) ------------------------"
	sed 's/"//g;s/_//g' <"$mypath"/"$name"_merged_data_even.txt | grep -v 'REJECT\|V2.x\|V2.y\|V3\|V4\|-' | awk '{ gsub("chr","")}1' | awk '{ gsub(":","")}1' > "$mypath"/"$name"_merged_data_even1.txt
	awk '{ $1=$2=$3=$4=$5=$6=$7=$8=""; print $0 }' "$mypath"/"$name"_merged_data_even1.txt | sed 's/./& /g' | sed 's/ //g;s/^ //g' > "$mypath"/"$name"_merged_data_even2.txt #insert spaces between characters, excluding the beginning of line
	sed 's/"//g' <"$mypath"/"$name"_merged_data_odd.txt  | grep -v 'REJECT\|V2.x\|V2.y\|V3\|V4\|-' > "$mypath"/"$name"_merged_data_odd1.txt
	awk '{ $1=$2=$3=$4=$5=$6=$7=$8=""; print $0 }' "$mypath"/"$name"_merged_data_odd1.txt | sed 's/./& /g' | sed 's/ //g;s/^ //g' > "$mypath"/"$name"_merged_data_odd2.txt #remove all spaces between characters, including the beginning of line

	echo "--------------------Complement of bases ------------------------"
	sed 's/T/W/g;s/C/X/g;s/G/Y/g;s/A/Z/g' <"$mypath"/"$name"_merged_data_odd2.txt | sed 's/W/A/g;s/X/G/g;s/Y/C/g;s/Z/T/g' > "$mypath"/"$name"_merged_data_odd2a.txt
	sed 's/T/W/g;s/C/X/g;s/G/Y/g;s/A/Z/g' <"$mypath"/"$name"_merged_data_even2.txt | sed 's/W/A/g;s/X/G/g;s/Y/C/g;s/Z/T/g' > "$mypath"/"$name"_merged_data_even2a.txt

	echo "--------------------Continue cleaning up R output files (odd and even files) ------------------------"
	awk '{ print $2 " " $3 " " $4 }' "$mypath"/"$name"_merged_data_even1.txt > "$mypath"/"$name"_merged_data_id.txt
	paste "$mypath"/"$name"_merged_data_id.txt "$mypath"/"$name"_merged_data_even2a.txt > "$mypath"/"$name"_merged_data_even3.txt
	awk '{ print $1 " " $2 " " $3 " " $4 }' "$mypath"/"$name"_merged_data_even3.txt > "$mypath"/"$name"_merged_data_even4.txt
	paste "$mypath"/"$name"_merged_data_id.txt "$mypath"/"$name"_merged_data_odd2a.txt > "$mypath"/"$name"_merged_data_odd3.txt
	awk '{ print $1 " " $2 " " $3 " " $4 }' "$mypath"/"$name"_merged_data_odd3.txt > "$mypath"/"$name"_merged_data_odd4.txt

	echo "-----------------Sort by mRNA position --------------------"

	sort -n -k 2 "$mypath"/"$name"_merged_data_even4.txt > "$mypath"/"$name"_merged_data_even5.txt
	sort -n -k 2 "$mypath"/"$name"_merged_data_odd4.txt > "$mypath"/"$name"_merged_data_odd5.txt

	echo "--------------------Populate the rest of sequence with reference-------------------"
	awk '{ gsub("NT",$3,$4); print $0}' "$mypath"/"$name"_merged_data_odd5.txt  | awk '{ gsub("NA","?",$3); print $0}' | grep -v "?" > "$mypath"/"$name"_merged_data_odd6.txt
	awk '{ gsub("NT",$3,$4); print $0}' "$mypath"/"$name"_merged_data_even5.txt | awk '{ gsub("NA","?",$3); print $0}' | grep -v "?" > "$mypath"/"$name"_merged_data_even6.txt

	echo "--------------------Create spaces between bases --------------------"
	awk '{ $1=$2=$3=""; print $0 }' "$mypath"/"$name"_merged_data_odd6.txt | sed 's/./& /g' > "$mypath"/"$name"_merged_data_odd7.txt #insert spaces between characters, excluding the beginning of line
	awk '{ $1=$2=$3=""; print $0 }' "$mypath"/"$name"_merged_data_even6.txt | sed 's/./& /g' > "$mypath"/"$name"_merged_data_even7.txt #insert spaces between characters, excluding the beginning of line
else
if [[ $response = "n" ]];
then 
	echo "--------------------Clean up R output files (odd and even) ------------------------"
	sed 's/"//g;s/_//g' <"$mypath"/"$name"_merged_data_even.txt | grep -v 'REJECT\|V2.x\|V2.y\|V3\|V4\|-' | awk '{ gsub("chr","")}1' | awk '{ gsub(":","")}1' > "$mypath"/"$name"_merged_data_even1.txt
	awk '{ $1=$2=$3=$4=$5=$6=$7=$8=""; print $0 }' "$mypath"/"$name"_merged_data_even1.txt | sed 's/./& /g' | sed 's/ //g;s/^ //g' > "$mypath"/"$name"_merged_data_even2.txt #insert spaces between characters, excluding the beginning of line
	sed 's/"//g' <"$mypath"/"$name"_merged_data_odd.txt  | grep -v 'REJECT\|V2.x\|V2.y\|V3\|V4\|-' > "$mypath"/"$name"_merged_data_odd1.txt
	awk '{ $1=$2=$3=$4=$5=$6=$7=$8=""; print $0 }' "$mypath"/"$name"_merged_data_odd1.txt | sed 's/./& /g' | sed 's/ //g;s/^ //g' > "$mypath"/"$name"_merged_data_odd2.txt #remove all spaces between characters, including the beginning of line

	echo "--------------------Continue cleaning up R output files (odd and even files) ------------------------"
	awk '{ print $2 " " $3 " " $4 }' "$mypath"/"$name"_merged_data_even1.txt > "$mypath"/"$name"_merged_data_id.txt
	paste "$mypath"/"$name"_merged_data_id.txt "$mypath"/"$name"_merged_data_even2.txt > "$mypath"/"$name"_merged_data_even3.txt
	awk '{ print $1 " " $2 " " $3 " " $4 }' "$mypath"/"$name"_merged_data_even3.txt > "$mypath"/"$name"_merged_data_even4.txt
	paste "$mypath"/"$name"_merged_data_id.txt "$mypath"/"$name"_merged_data_odd2.txt > "$mypath"/"$name"_merged_data_odd3.txt
	awk '{ print $1 " " $2 " " $3 " " $4 }' "$mypath"/"$name"_merged_data_odd3.txt > "$mypath"/"$name"_merged_data_odd4.txt

	echo "-----------------Sort by mRNA position --------------------"

	sort -n -k 2 "$mypath"/"$name"_merged_data_even4.txt > "$mypath"/"$name"_merged_data_even5.txt
	sort -n -k 2 "$mypath"/"$name"_merged_data_odd4.txt > "$mypath"/"$name"_merged_data_odd5.txt

	echo "--------------------Populate the rest of sequence with reference-------------------"
	awk '{ gsub("NA",$3,$4); print $0}' "$mypath"/"$name"_merged_data_odd5.txt  | awk '{ gsub("NA","?",$3); print $0}' | grep -v "?" > "$mypath"/"$name"_merged_data_odd6.txt
	awk '{ gsub("NA",$3,$4); print $0}' "$mypath"/"$name"_merged_data_even5.txt | awk '{ gsub("NA","?",$3); print $0}' | grep -v "?" > "$mypath"/"$name"_merged_data_even6.txt

	echo "--------------------Create spaces between bases --------------------"
	awk '{ $1=$2=$3=""; print $0 }' "$mypath"/"$name"_merged_data_odd6.txt | sed 's/./& /g' > "$mypath"/"$name"_merged_data_odd7.txt #insert spaces between characters, excluding the beginning of line
	awk '{ $1=$2=$3=""; print $0 }' "$mypath"/"$name"_merged_data_even6.txt | sed 's/./& /g' > "$mypath"/"$name"_merged_data_even7.txt #insert spaces between characters, excluding the beginning of line
else 
      echo "ERROR!"
fi
	fi

echo "--------------------Transpose nucleotides to sequence data -----------------------"

perl -anle '
    $l[$_] .= $F[$_] for 0..$#F;
    END {
        print join " ", split // for @l;
    }' "$mypath"/"$name"_merged_data_even7.txt > "$mypath"/"$name"_merged_data_even8.txt

perl -anle '
    $l[$_] .= $F[$_] for 0..$#F;
    END {
        print join " ", split // for @l;
    }' "$mypath"/"$name"_merged_data_odd7.txt > "$mypath"/"$name"_merged_data_odd8.txt

echo "----------------- Combine chromosome 1 with chromosome 2 for same individual -------------------"

sed = "$mypath"/"$name"_merged_data_even8.txt | sed 'N;s/\n/ /g' > "$mypath"/"$name"_merged_data_even9.txt
sed = "$mypath"/"$name"_merged_data_odd8.txt | sed 'N;s/\n/ /g' > "$mypath"/"$name"_merged_data_odd9.txt

cat "$mypath"/"$name"_merged_data_odd9.txt "$mypath"/"$name"_merged_data_even9.txt > "$mypath"/"$name"_whole.txt
sort -n -k 1 "$mypath"/"$name"_whole.txt > "$mypath"/"$name"_whole1.txt 
awk '{ $1=""; print $0 }' "$mypath"/"$name"_whole1.txt | sed 's/0/?/g' | sed 's/ //g;s/^ //g' > "$mypath"/"$name"_whole2.txt
#sed -i "s/[0-9]//g" "$mypath"/"$name"_whole2.txt

echo "----------------- Make directory for gene  -------------------"

mkdir -p "$mypath"/"$name"/
mv "$mypath"/"$name"_whole2.txt "$mypath"/"$name"/
mv "$mypath"/"$name"_merged_data_odd6.txt "$mypath"/"$name"/
mv "$mypath"/"$name"_merged_data_even6.txt "$mypath"/"$name"/
mv "$mypath"/"$name"_phased.vcf.gz "$mypath"/"$name"/

echo "------------------- duplicate each line of sample labels -----------"

awk '{ print $1 " " $1 }' "$mypath"/"$name".impute.hap.indv | sed 's/ /\n/g' > "$mypath"/"$name".tfam_dup
sed -e 's/ /\n/g;1~2s/$/_1/g' "$mypath"/"$name".tfam_dup > "$mypath"/"$name".tfam_dup2
paste "$mypath"/"$name".tfam_dup2 "$mypath"/"$name"/"$name"_whole2.txt > "$mypath"/"$name"/"$name"_whole3.txt
paste "$mypath"/new_merged1.txt "$mypath"/"$name"/"$name"_whole3.txt > "$mypath"/"$name"/"$name"_whole4.txt
#sed 's/_1//g' < "$mypath"/"$name"/"$name"_whole4.txt | awk '{ print $1 "_" $2 " " $3 }' > "$mypath"/"$name"/"$name"_whole5.txt
awk '{ print $1 "_" $2 " " $3 }' "$mypath"/"$name"/"$name"_whole4.txt > "$mypath"/"$name"/"$name"_whole5.txt
awk '{ print $1 "\t" $2 }' "$mypath"/"$name"/"$name"_whole5.txt | sed 's/\t/\n/g' | sed 's/^NA/>NA/g;s/^HG/>HG/g' | grep -A1 '>' > "$mypath"/"$name"/"$name".fasta #we are taking the first 2802 haplotypes for 1401 individuals

echo "-----------------Remove excess files--------------------------------"

rm "$mypath"/"$name"_converted_data1.txt
rm "$mypath"/"$name"_merged*
rm "$mypath"/"$name"_step*
rm "$mypath"/"$name"_whole*
rm "$mypath"/"$name"/"$name"_phased*
rm "$mypath"/"$name"3.*
rm "$mypath"/"$name"5.*
rm "$mypath"/"$name"_phased.*
rm "$mypath"/"$name".impute*
rm "$mypath"/"$name".phased*
rm "$mypath"/"$name".tfam*
rm "$mypath"/"$name".log
rm "$mypath"/RScript_"$name"*

echo "-----------------Create input for MACPRF--------------------------------"

touch "$mypath"/"$name"/"$name"_MACPRF.fas
cp "$mypath"/"$name"/"$name".fasta "$mypath"/"$name"/"$name"_MACPRF.fas
echo "Program Complete. Good-bye!"
done