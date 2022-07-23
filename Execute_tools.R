#######################################################################################
# Use the tools for initial process of 2bRAD sequence data ############################
#######################################################################################
## Edit only the following chunk of code ##
prefix<- "PlateX_test"  # a prefix to write in the name of output files
path.to.RawData<- "/home/eduardo/Documents/Test/RawData" # the directory where are located the raw sequences as provided by sequencer (fastq.gz format)
path.to.Scripts<- "/home/eduardo/Documents/Test/scripts" # a directory where we will save the shell scripts but also the perl script
path.to.Data<- "/home/eduardo/Documents/Test/data" # a directory to save some output data including the md5sum and a list of proceced samples
path.to.md5sum<- "/home/eduardo/Documents/Test/RawData/md5sum.txt" # the path to the md5sum file provided by the sequencer it has two collumns 1. md5sum and 2. Name_of_files
path.to.plate.map<-"/home/eduardo/Documents/Test/dataplate10_map.csv" # the path to the file with a plate map (8x12) in csv format and with no header and no row.names
path.to.demultiplexed <-"/home/eduardo/Documents/Test/demultiplexed" # this path will be created after running the demultiplex.sh script it will contain the demultiplexed sequences in .fq format
########################################################################################
## load functions
source("~/Documents/Github/2bRAD-tools/2bRAD_tools_functions.R")

## If you define  all the object above the functions should run without providing any argument
## to see what arguments are being passed to each function execute:
args(md5sum.check)
args(gunzip.fastq.gz)
args(write.demultiplex.sh)
args(write.fastq.sh)

# This function will generate a md5sum of the files provided by the sequences in fast.gz format, and compare to the md5sum provided  by the sequencer in a txt file
md5sum.check() 
# This function will decompress fastq.qz files to fastq files
gunzip.fastq.gz()
# This function will write  shell script that can be executed for demultiplexing samples, for running the shell script one must have the ... file in in the directory specified in path.to.Scripts directory
write.demultiplex.sh()
# This function will write a shell script for executing fastx_tools to rename fastq files based on the positions provided in a plate_map.csv, this will also perform a quality filter
write.fastq.sh()



