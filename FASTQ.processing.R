### Script to process RAW amplicon data ###

# First, we need to extract the R1 and R2 fastq files from each folder
# and then we need to place these in a single folder. This new folder will be 
# placed in our working directory for R, which will allow us to process the 
# reads using the dada2 workflow in R. 

find ~/Documents/Plastic_Expt/DefaultProject/ -type f -name '*16S*' -exec mv {} ~/Documents/Plastic_Expt/16S_Amps/ \;
find ~/Documents/Plastic_Expt/DefaultProject/ -type f -name '*ITS*' -exec mv {} ~/Documents/Plastic_Expt/ITS_Amps/ \;
