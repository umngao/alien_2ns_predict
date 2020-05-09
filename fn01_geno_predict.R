
geno_predict = function(atag.file, wtag.file, tbt, log, reads.cut=10^5){
    library(data.table)
    ## This function requires data.table package
    ## It takes as input the 
    ## 1. alien tag file: 'data/tags.alien.txt'
    ## 2. wheat.tag file: 'data/tags.wheat.txt'
    ## 3. TBT (tags by taxa file for selected tags): 'data/All.seltag.TBT.txt'
    ## 4. per sample read depth log file: 'data/All_ReadsPerSample.log'
    
    ######## step (1) read in input tags, tagsByTaxa (TBT, tbt) file and raw.cnts file
    alien.tag = fread(atag.file,head=T)
    colnames(alien.tag)[1] ='Tags'
    wheat.tag = fread(wtag.file,head=T)
    colnames(wheat.tag)[1] ='Tags'
    aTag.num = nrow(alien.tag)
    wTag.num = nrow(wheat.tag)
    
    TBT = fread(tbt, head=T)  # TagsByTaxa file
    colnames(TBT)[1] = 'Tags'
    reads.cnt.raw = data.table(read.table(log,header=T,stringsAsFactors=F))
    reads.cnt.raw = reads.cnt.raw[goodBarcodedReads >= reads.cut]
    
    
    #### step (2) first normalization based on raw reads sequencing depth (in million)
    alien.tag[,specific:='alien']; wheat.tag[,specific:='wheat']; 
    tags.all = rbind(alien.tag,wheat.tag)
    colnames(tags.all)[1] = 'Tags'
    DT = tags.all[TBT, on='Tags', nomatch=0]
    DT.sub = DT[,c(reads.cnt.raw$FullSampleName),with=F]
    DT.sub.norm1 = sweep(DT.sub, 2, reads.cnt.raw$goodBarcodedReads/10^6, `/`) 
    d.tmp = data.table(DT[,.(Tags, specific)],DT.sub.norm1)
    
    
    #### step (3) then calculate per type (alien or wheat) sum
    calc_sum_based_on_specific = function(dt, specific = d.tmp$specific ){
        cnt1 = sum(dt[specific =='alien'])
        cnt2 = sum(dt[specific =='wheat'])
        return(c(cnt1,cnt2))
    }
    results = apply(d.tmp[,-c(1:2)],2,calc_sum_based_on_specific)
    cnt.alien = results[1,]
    cnt.wheat = results[2,]
    DT = data.table(line=names(cnt.alien),cnt.alien=cnt.alien,cnt.wheat=cnt.wheat)
    
    
    ##### step (4) Further normalization by relative alien Tag or wheat Tag numbers
    DT[,norm1:=cnt.alien/aTag.num][,norm2:=cnt.wheat/wTag.num]
    DT=DT[!(norm1+norm2)==0]
    DT[,rf:=norm1/(norm1+norm2)]
    
    ######  Step (5) picking colors and predict, and adjust prediction results
    DT[,predict:=norm1>norm2]
    DT[,COL:=ifelse(predict==T, 'red','blue')]
    DT[,adjusted.predict:=predict]
    DT[(cnt.alien+cnt.wheat)<20, adjusted.predict:=NA]
    DT[rf>0.2 & rf<0.7, adjusted.predict:=NA]
    DT[,genotype:=ifelse(predict==T, 'alien','wheat')]
    DT[,zygosity:=ifelse(rf>0.2 & rf<0.7, "heterozygous", "homozygous")]
    DT[]
}
