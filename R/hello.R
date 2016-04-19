plotAxtl<- function(gene,axtl=annotated_tT) {
  #copy the name of the gene for the plot title
  genet<-gene
  #Add carot and dollar sign to gene to ensure nothing else matches gene name in search
  gene<-paste0("^",gene,"$")
  #All data columns have "dpi" in the name.  This searches the axtl data for "gene" 
  #and assigns the columns that match "dpi" to the temporary data frame temp.
  temp<-as.data.frame(axtl[grep(gene,axtl$Gene),grep("dpi",colnames(axtl))])
  #Create a new column in temp with the probe.ids which are currently the rownames.
  temp$probe <-rownames(temp)
  #Add a column for this dataset so they can be split out later in the plot
  temp$species<-"Axolotl"
  #This converts the data from "wide" to "long" format, for easier plotting.
  temp<-gather(temp,variable,value,-probe,-species)
  #The new column names will become the axis and legend labels
  colnames(temp)<-c("probe","species","DPI","logFC")
  
  #Now do humanGSE28914
  temp2<-as.data.frame(humanGSE28914[grep(gene,humanGSE28914$symbol),grep("dpi",colnames(humanGSE28914))])
  temp2$probe <-rownames(temp2)
  temp2$species<-"GSE28914"
  temp2<-gather(temp2,variable,value,-probe,-species)
  colnames(temp2)<-c("probe","species","DPI","logFC")
  
  #Now do humanGSE28914
  temp3<-as.data.frame(humanGSE50425[grep(gene,humanGSE50425$symbol),grep("dpi",colnames(humanGSE50425))])
  temp3$probe <-rownames(temp3)
  temp3$species<-"GSE50425"
  temp3<-gather(temp3,variable,value,-probe,-species)
  colnames(temp3)<-c("probe","species","DPI","logFC")
  
  temp4 <-rbind(temp,temp2,temp3)
  p<-ggplot(temp4,aes(x=DPI,y=logFC,colour=probe,group=probe))
  p+geom_line()+geom_point()+theme_bw()+ggtitle(genet)+facet_grid(. ~ species, scales="free_x")
  
}

#Define a new function that will accept Genes or Probes or long lists of genes/probes.
plotAxtlv3<- function(genelist,v=T) {
  #create a list to hold all the data as we loop through the gene list
  data<-list()
  
  #start of loop
  for (gene in genelist) {
  #check to see if gene is in the axolotl list
  if (gene %in% axtl$Gene) {
    #copy the name of the gene for the plot title
    genet<-gene
    #Add carot and dollar sign to gene to ensure nothing else matches gene name in search
    gene<-paste0("^",gene,"$")
    #All data columns have "dpi" in the name.  This searches the axtl data for "gene" 
    #and assigns the columns that match "dpi" to the temporary data frame temp.
    temp<-as.data.frame(axtl[grep(gene,axtl$Gene),grep("dpi",colnames(axtl))])
    #This will use use probe name if it is not
    gene<-paste0("^",gene,"$")
    is_a_probe<-F
  } else if (gene %in% rownames(axtl)) {
    temp<-axtl[gene,2:4]
    genet<-axtl[gene,"symbol"]
    gene<-paste0("^",genet,"$")
    is_a_probe<-T
  } else {return(paste0(gene," not found as a gene name or probeid"))}
  
  #Look up human gene values only if there is a human homolog for the axolotl gene
if (nchar(gene) > 2) {
  #Now do humanGSE28914
  temp2<-as.data.frame(humanGSE28914[grep(gene,humanGSE28914$symbol),grep("dpi",colnames(humanGSE28914))])
  
  #Now do humanGSE28914
  temp3<-as.data.frame(humanGSE50425[grep(gene,humanGSE50425$symbol),grep("dpi",colnames(humanGSE50425))])
  
  #combine two microarrays to get pseudohuman
  print("Warning:  Combining GSE28914 and GSE50425 data.")
  human<-matrix(rep(0,4))
  human[1]<-median(temp2$acute_dpi,na.rm=T)
  human[2]<-median(temp2$day3_dpi,na.rm=T)
  human[3]<-median(temp3$day14_dpi,na.rm=T)
  human[4]<-median(temp3$day21_dpi,na.rm=T)
  
  humanDF<-as.data.frame(human)
  humanDF$species<-"Human"
  humanDF$DPI<-c("0_dpi","2_dpi","14_dpi","21_dpi")
  humanDF$DPI<-factor(humanDF$DPI,levels=c("0_dpi","2_dpi","14_dpi","21_dpi"))
  colnames(humanDF)<-c("logFC","Species","DPI")
  
  } else { 
    human<-matrix(rep(NA,4))
    humanDF<-NULL 
    }

#  } else { humanDF<-as.data.frame(matrix(rep(0,4)))}
  
#print(paste0(nrow(temp)," probe(s) found for gene ",genet,":  ",rownames(temp)))
  for (probe in rownames(temp)) {
    #print(temp[i,])
    #probe<-rownames(temp[i,])
    
    axolotl<-as.numeric(c(0,temp[probe,]))
    if (nchar(gene) > 2) { y<-as.numeric(cor(human,axolotl)) } else {  y<-NA  }
#   data[[probe]]<-c(y,axolotl,human)
    data[[probe]]<-c(y,axtl[probe,"Gene"],axolotl[2:4],human,axtl[probe,"Name"],axtl[probe,"Source.Seq"])
    
    if (v) {
    axolotlDF<-as.data.frame(axolotl)
    axolotlDF$species<-"Axolotl"
    axolotlDF$DPI<-c("0_dpi","2_dpi","14_dpi","21_dpi")
    axolotlDF$DPI<-factor(axolotlDF$DPI,levels=c("0_dpi","2_dpi","14_dpi","21_dpi"))
    colnames(axolotlDF)<-c("logFC","Species","DPI")
    temp5<-rbind(axolotlDF,humanDF)
    #temp5<-melt(temp5,id.vars="Species","DPI")
    #title<-paste0(genet," (",probe,") Pearson=",y)
    title1<-axtl[probe,"description"]
    if (is_a_probe) { title2<-paste0(probe, " Pearson's Coefficient=",round(y,4)) } else {
      title2<-paste0(genet,"  (",probe,")"," Pearson's Coefficient=",round(y,4))}
    #title2<-paste0("Pearson=  ",round(y,4))
    print(ggplot(temp5,aes(x=DPI,y=logFC,colour=Species,group=Species))+geom_line()+geom_point()+theme_bw()+ggtitle(bquote(atop(.(title1), atop(italic(.(title2)), "")))))
    }
  }
  }
  data2<-do.call("rbind",data)
#  colnames(data2)<-c("corr","axo_0_dpi","axo_2_dpi","axo_14_dpi","axo_21_dpi","human_0_dpi","human_3_dpi","human_14_dpi","human_21_dpi")
  colnames(data2)<-c("corr","symbol","axo_2_dpi","axo_14_dpi","axo_21_dpi","human_0_dpi","human_3_dpi","human_14_dpi","human_21_dpi","name","description")
  if (v) {return(paste0(nrow(temp)," probe(s) found for gene ",genet,":  ",paste(rownames(temp),collapse=",")))} else {return(data2) }
}
