#!/bin/bash
	read -p "Enter your chromosome of choice: " j
	echo "You have selected chromosome $j!"
while [[ "$j" -gt 0 && "$j" -lt 24 ]];
do
	echo "You are now downloading chromosome $j:$k-$l"

wget --continue ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr"$j".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget --continue ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr"$j".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

if [ -f ALL.chr"$j".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ];
then
	echo "Your download is complete!"
	break
	continue
fi
done
	read -r -p "Do you want to extract your polymorphisms of interest? [y/n]" response
if [[ $response = "n" ]]
then
	echo "Thank you. Good-bye"
	exit
elif [[ $response = "y" ]];
then 
	read -p "Enter the name of your gene: " i
	read -p "Enter the start position of your gene: " k
	read -p "Enter the end position of your gene: " l
fi
mypath=/home/wgriffith/scripts
while [[ "$k" -lt "$l" ]];
do
	echo "Your start and end positions are correct"
	"$mypath"/vcftools_0.1.13/bin/./vcftools --gzvcf "$mypath"/ALL.chr"$j".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --min-alleles 2 --max-alleles 2 --remove-filtered-all --remove-indels --max-missing 1 --chr "$j" --from-bp "$k" --to-bp "$l" --keep "$mypath"/292_illumina_samples1.txt --recode --out "$mypath"/1000Genomes_P3_"$i"

if [ -f 1000Genomes_P3_"$i".recode.vcf ];
then
	echo "You have successfully extracted your polymorphisms of interest!"
	break
	continue
fi
done
if [ ! -f 1000Genomes_P3_"$i".recode.vcf ];
then
	echo "You have not successfully extracted your polymorphisms of interest! Check your start and end positions."
	exit	
fi