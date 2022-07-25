# 2bRAD-tools
A set of tools to process 2bRAD sequencing, compatible with the workflow in https://github.com/Kenizzer/2bRAD-Edwards-Lab
## md5sum.check() 
This function will generate a md5sum of the files provided by the sequences in fast.gz format, and compare to the md5sum provided  by the sequencer in a txt file

## gunzip.fastq.gz()
This function will decompress fastq.qz files to fastq files

## write.demultiplex.sh()
This function will write  shell script that can be executed for demultiplexing samples, for running the shell script one must have the ... file in in the directory specified in path.to.Scripts directory

## write.fastq.sh()
This function will write a shell script for executing fastx_tools. Such a bash script will rename fastq files with names provided in a plate_map.csv depending on the row/colum possition (this asummes that adaptors used are the same as in https://github.com/Kenizzer/2bRAD-Edwards-Lab), the bash script will also perform a quality filter.
