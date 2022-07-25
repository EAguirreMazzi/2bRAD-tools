
#### The following script has some functions to automate some task when working with 2b-RAD sequences
#### Step1 : Md5summary check ####
md5sum.check<-function(name=prefix, RawData=path.to.RawData,
                       Data=path.to.Data, md5sum=path.to.md5sum){
  library(tools)
  list.of.files<-list.files(RawData)
  md5.ori<-read.table(md5sum)
  names(md5.ori)<-c("md5.ori","files" )
  md5.new<-md5sum(list.of.files[grep(".fastq",list.of.files)])
  md5.new<-data.frame(md5.new)
  md5.new$files<-rownames(md5.new)
  md5.check <-merge(md5.new, md5.ori, by= c("files"))
  write.csv(md5.check, 
          paste(Data,"/", name ,"_md5sum_report.txt",sep = ""))
}

#### Step2: decompress fastq.gz ####
gunzip.fastq.gz<-function(RawData=path.to.RawData){
  library(R.utils)
  setwd(RawData)
  list.of.files<-list.files(RawData)
  files<-paste(RawData,list.of.files[grep(".fastq.gz",list.of.files)],sep = "/")
  lapply(files,gunzip)
}

#Step3: prepare shell scripts for de-multiplexing (use decompressed file names)
write.demultiplex.sh<-function(name=prefix,RawData=path.to.RawData, 
                               Scripts=path.to.Scripts){
  list.of.files<-list.files(RawData)
  decompressed_files<-list.of.files[grep(".fastq",list.of.files)]
  demultiplex.sh<-c("#!/bin/bash",
                    paste("cd ", RawData,sep = ""),
                    paste("cp", Scripts, "/trim2bRAD_2barcodes_noAdap.pl ." ,sep=""),
                    "chmod +x trim2bRAD_2barcodes_noAdap.pl",
                    paste("perl trim2bRAD_2barcodes_noAdap.pl fastq=",
                          decompressed_files,
                          " site=",dQuote('.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}', q=getOption("double")),
                          " barcode2=",dQuote("[ATGC]{4}",q=getOption("double")),
                          sep=""),
                    "mkdir demultiplexed",
                    "mv *TGGT.fq demultiplexed",
                    "mv *AGAC.fq demultiplexed",
                    "mv *ACCA.fq demultiplexed",
                    "mv *AGTG.fq demultiplexed",
                    "mv *CATC.fq demultiplexed",
                    "mv *GTGA.fq demultiplexed",
                    "mv *TCAG.fq demultiplexed",
                    "mv *GCTT.fq demultiplexed",
                    "mv *CTAC.fq demultiplexed",
                    "mv *TGTC.fq demultiplexed",
                    "mv *TCAC.fq demultiplexed",
                    "mv *GACT.fq demultiplexed",
                    "rm *.fq",
                    "rm trim2bRAD_2barcodes_noAdap.pl")
  write.table(demultiplex.sh,
              paste(Scripts, "/", name, "_demultiplex.sh",sep = ""),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

############################################################################################
#### Step4 : after running the bash script above (_demultiplex.sh) the following script will create a bash script for renaming files and applying some quality control filters
############################################################################################
write.fastq.sh <-function(plate.map=path.to.plate.map,
                          demultiplexed=path.to.demultiplexed,
                          Scripts=path.to.Scripts,
                          Data=path.to.Data,
                          name=prefix,
                          q=30,
                          p=100){
  library(reshape2)
  plate<-as.matrix(read.csv(plate.map, header = F))
  plate.list <-melt(plate)
  names(plate.list)<-c("Rows","Cols","sample")
  plate.list$Rows<-paste("row",plate.list$Rows,sep = "")
  plate.list$Cols<-gsub("V","col",plate.list$Cols)
  plate.list$Ill_adaptor<-c(rep("TGGT",8),rep("AGAC",8),rep("ACCA",8),rep("AGTG",8),rep("CATC",8),rep("GTGA",8),rep("TCAG",8),rep("GCTT",8),rep("CTAC",8),rep("TGTC",8),rep("TCAC",8),rep("GACT",8))  
  plate.list$row_adap <- paste(plate.list$Rows,plate.list$Ill_adaptor, sep = "_")
  
  demultiplexed.files<- list.files(demultiplexed)
  Ill_adaptor<- rep("",length(demultiplexed.files))
  for (i in 1:length(demultiplexed.files)) {
    adapList<-c("TGGT", "AGAC", "ACCA", "AGTG", "CATC", "GTGA",
               "TCAG", "GCTT", "CTAC", "TGTC", "TCAC", "GACT")
    for(adap in 1:12){
      index<-grep(adapList[adap],demultiplexed.files, ignore.case = T)
      Ill_adaptor[index]<-adapList[adap]
    }
  }
  
  row_position<-rep("",length(demultiplexed.files))
    for (i in 1:length(demultiplexed.files)) {
    rowList<- paste("row", 1:8, sep = "")
    for(row in 1:8){
      index<-grep(rowList[row],demultiplexed.files, ignore.case = T)
      row_position[index]<-rowList[row]
    }
  }
  
  demultiplexed.files<- data.frame(demultiplexed.files=demultiplexed.files,
                                   row_adap= paste(row_position,Ill_adaptor, sep = "_"))
  List<-merge(demultiplexed.files, plate.list, by= "row_adap")
  core.code<-paste("fastq_quality_filter -i ", List$demultiplexed.files," -o ",List$sample,".fastq -q ", q ," -p ", p , sep="")
  ####PREPARE SHELL SCRIPT FOR RENAMING AND FILTERING DEMULTIPLEXED FILES WITH FASTQC####
  code<-c("#!/bin/bash",
          paste("cd ", demultiplexed, sep = ""),
          core.code,
          paste("mkdir ",demultiplexed, "/named_filtered",sep = ""),
          paste("mv ", demultiplexed, "*.fastq ", demultiplexed ,"/named_filtered",sep = ""),
          paste("echo ",prefix," ready",sep = ""))
  write.table(code,
              file = paste(Scripts,"/", name,"_fastx_qcfilter_naming.sh",sep = ""),
              quote = FALSE, row.names = FALSE,col.names = FALSE)
  write.csv(List, paste(Data,"/",name,"_list_of_samples.csv",sep = ""))
}
