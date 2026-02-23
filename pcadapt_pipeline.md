conda activate test_env

#view vcf file header
bcftools view -h .vcf.gz

##change headers
#extract header
bcftools view -h chr17.GFM.withId.vcf.gz > headerGFM.txt   
bcftools view -h chr17.WMA.withID.vcf.gz > headerWMA.txt 

#move txt file
mv headerGFM.txt /home/las80898/Mallard/
mv headerWMA.txt /home/las80898/Mallard/

#manually edit header
##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Average mapping quality">  
note: first attempt, i didn't change format, just type=float

#move txt file back
mv /home/las80898/Mallard/headerGFM.txt /home/las80898/mallard_data 
mv /home/las80898/Mallard/headerWMA.txt /home/las80898/mallard_data 

#reheader the file
bcftools reheader -h headerGFM.txt -o chr17.GFMIDreheader.vcf.gz chr17.GFM.withId.vcf.gz   
bcftools reheader -h headerWMA.txt -o chr17.WMAIDreheader.vcf.gz chr17.WMA.withID.vcf.gz   
##

#create index files
bcftools index chr17.GFMIDreheader.vcf.gz 
bcftools index chr17.WMAIDreheader.vcf.gz 

#merge vcf files
bcftools merge chr17.GFMIDreheader.vcf.gz chr17.WMAIDreheader.vcf.gz -Oz -o chr17merge.vcf.gz

#convert vcf to bed
plink2 --vcf chr17merge.vcf.gz --make-bed --out chr17
