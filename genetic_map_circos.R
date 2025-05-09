# before starting, i create a file analogous to the one in Carolina's paper for spruce. that file is copied here as Consensus_maps_SPRUCE.txt, downloaded on 2025-24-24 from https://zenodo.org/records/1472158

# for n in {1..12}; do awk -v lg="$n" '{print lg "\t" $0}' consensus_LG"$n".tsv; done > consensus_LGs.tsv
## i change the first "1" in the first line to "LG", and then delete the other 11 header lines before the start of a new LG

#grep -v "#" AxiomGT1.calls_Y3088_HetMarkedMissing_newCoords_newHeader_sorted_variable.vcf | cut -f1,3 > tmp.tsv
#grep -v "#" AxiomGT1.calls_AC1017_strict_HetMarkedMissing_newCoords_newHeader_sorted_variable.vcf | cut -f1,3 >> tmp.tsv 
#sort -u tmp.tsv > tmp_unique.tsv
#awk 'NR==FNR {chr[$2]=$1; next} {if ($2 in chr) print $0 "\t" chr[$2]; else print $0 "\tNA"}' tmp_unique.tsv consensus_LGs.tsv > consensus_LGs_chr.tsv
#change last column name to "chr"

## source code downloaded on 2025-04-24 from https://github.com/parkingvarsson/HaploidSpruceMap/blob/master/map_analysis/omicCircos.R
### heavily modified, the original code is kept here as omicCircos_SPRUCE.R

read.table("./consensus_LGs_chr.tsv",head=T)->maps

library (OmicCircos)
options (stringsAsFactors=FALSE) ;

##Segment frame
cons.name <- paste("LG",maps$LG,sep=" ")
cons.start <- as.character(round(maps$consensus),digit=3)
cons.end <- as.character(round(maps$consensus),digit=3)
cons.f <- data.frame(seg.name=cons.name, seg.start=cons.start,seg.end=cons.end,the.v=runif(length(cons.end)),Note=as.character(maps$marker))
#####
#> head(cons.f)
#seg.name seg.start seg.end       the.v         Note
#1     LG 1         0       0 0.002544849 AX-601671694
#2     LG 1         0       0 0.306732609 AX-601671709
#3     LG 1         0       0 0.127236336 AX-389938569
#4     LG 1         0       0 0.003909409 AX-601671720
#5     LG 1         0       0 0.853212124 AX-602442804
#6     LG 1         0       0 0.576358559 AX-601115205
#####
## i don't really understand what's the column "the.v"

db <- segAnglePo(cons.f,seg=unique(cons.name),angle.start = 3,angle.end = 357);
####
#> head(db)
#seg.name angle.start        angle.end          seg.sum.start seg.sum.end seg.start seg.end
#[1,] "LG 1"   "273"              "302.692186266772" "0"           "114"       "0"       "114"  
#[2,] "LG 2"   "304.692186266772" "332.300710339384" "114"         "220"       "0"       "106"  
#[3,] "LG 3"   "334.300710339384" "360.867403314917" "220"         "322"       "0"       "102"  
#[4,] "LG 4"   "362.867403314917" "392.820047355959" "322"         "437"       "0"       "115"  
#[5,] "LG 5"   "394.820047355959" "422.949486977111" "437"         "545"       "0"       "108"  
#[6,] "LG 6"   "424.949486977111" "449.172059984215" "545"         "638"       "0"       "93"        
####

##Segment mapping
cons.v<-maps[,c(1,3)]
names(cons.v)<-c("seg.name","seg.pos")
cons.v$seg.pos<-as.character(cons.v$seg.pos)
cons.v$seg.name<-paste("LG",cons.v$seg.name,sep=" ")
####
#> head(cons.v)
#seg.name            seg.pos
#1     LG 1                  0
#2     LG 1              0.092
#3     LG 1 0.0930000000000035
#4     LG 1 0.0930000000000035
#5     LG 1 0.0930000000000035
#6     LG 1  0.367999999999998
####

###Segment links
cons.links <- NULL
cons.links$chr1 <-NULL
cons.links$pos1 <-NULL
cons.links$marker1 <-NULL
cons.links$chr2 <-NULL
cons.links$pos2 <-NULL
cons.links$marker2 <-NULL
cons.links$ann1 <- NULL
cons.links$ann2 <- NULL
for (i in 1:length(unique(maps$chr))){
  tmp <- maps[maps$chr == unique(maps$chr)[i],]
  if (dim(tmp)[1] > 1){
    for (j in 2:dim(tmp)[1]){
      cons.links$chr1 <- c(cons.links$chr1,as.character(paste("LG",tmp$LG[j-1],sep=" ")))
      cons.links$chr2 <- c(cons.links$chr2,as.character(paste( "LG", tmp$LG[j],sep=" ")))
      cons.links$pos1 <- c(cons.links$pos1,as.character(tmp$consensus[j-1]))
      cons.links$pos2 <- c(cons.links$pos2,as.character(tmp$consensus[j]))
      cons.links$marker1 <- c(cons.links$marker1,as.character(tmp$marker[j-1]))
      cons.links$marker2 <- c(cons.links$marker2,as.character(tmp$marker[j]))
      cons.links$ann1 <- c(cons.links$ann1,as.character(tmp$chr[j-1]))
      cons.links$ann2 <- c(cons.links$ann2,as.character(tmp$chr[j]))
    }
  }
}

cons.links <- as.data.frame(cons.links)
cons.links <- cons.links[,c(1,3,5,7,2,4,6,8)]
####
# head(cons.links)
#chr1               pos1      marker1     ann1 chr2               pos2      marker2     ann2
#1 LG 1                  0 AX-601671694 PS_chr01 LG 1              0.092 AX-601671709 PS_chr01
#2 LG 1              0.092 AX-601671709 PS_chr01 LG 1 0.0930000000000035 AX-389938569 PS_chr01
#3 LG 1 0.0930000000000035 AX-389938569 PS_chr01 LG 1 0.0930000000000035 AX-601671720 PS_chr01
#4 LG 1 0.0930000000000035 AX-601671720 PS_chr01 LG 1 0.0930000000000035 AX-602442804 PS_chr01
#5 LG 1 0.0930000000000035 AX-602442804 PS_chr01 LG 1  0.367999999999998 AX-601115205 PS_chr01
#6 LG 1  0.367999999999998 AX-601115205 PS_chr01 LG 1  0.367999999999998 AX-600017573 PS_chr01
####

intra_cons.links <- cons.links[cons.links$chr1 == cons.links$chr2,]
####
# head(intra_cons.links)
#chr1               pos1      marker1     ann1 chr2               pos2      marker2     ann2
# LG 1                  0 AX-601671694 PS_chr01 LG 1              0.092 AX-601671709 PS_chr01
#2 LG 1              0.092 AX-601671709 PS_chr01 LG 1 0.0930000000000035 AX-389938569 PS_chr01
# LG 1 0.0930000000000035 AX-389938569 PS_chr01 LG 1 0.0930000000000035 AX-601671720 PS_chr01
#4 LG 1 0.0930000000000035 AX-601671720 PS_chr01 LG 1 0.0930000000000035 AX-602442804 PS_chr01
#5 LG 1 0.0930000000000035 AX-602442804 PS_chr01 LG 1  0.367999999999998 AX-601115205 PS_chr01
#6 LG 1  0.367999999999998 AX-601115205 PS_chr01 LG 1  0.367999999999998 AX-600017573 PS_chr01
####

inter_cons.links <- cons.links[cons.links$chr1 != cons.links$chr2,]
####
#> head(inter_cons.links)
#hr1    pos1      marker1     ann1 chr2   pos2      marker2     ann2
#1265 LG 1 114.073 AX-613846573 PS_chr01 LG 4 21.145 AX-606636740 PS_chr01
#2620 LG 2 105.592 AX-606664754 PS_chr02 LG 4 21.145 AX-599615339 PS_chr02
#626 LG 2  59.572 AX-600441685 PS_chr06 LG 4 21.145 AX-599460881 PS_chr06
#2627 LG 4  21.145 AX-599460881 PS_chr06 LG 5 18.898 AX-601581620 PS_chr06
#2628 LG 5  18.898 AX-601581620 PS_chr06 LG 6      0 AX-606708523 PS_chr06
#4180 LG 2  74.625 AX-599724810 PS_chr09 LG 4 21.145 AX-601649940 PS_chr09
####

col.intra.links <- rep("grey90", nrow(intra_cons.links))  ## default to grey90

for (i in 1:nrow(intra_cons.links)) {
  if (abs(as.numeric(intra_cons.links$pos2[i]) - as.numeric(intra_cons.links$pos1[i])) > 5 &&
      as.numeric(gsub("LG ", "", intra_cons.links$chr1[i])) == as.numeric(gsub("PS_chr", "", intra_cons.links$ann1[i])) &&
      as.numeric(gsub("LG ", "", intra_cons.links$chr2[i])) == as.numeric(gsub("PS_chr", "", intra_cons.links$ann2[i]))) {
    col.intra.links[i] <- "grey35"
  }
  if (abs(as.numeric(intra_cons.links$pos2[i]) - as.numeric(intra_cons.links$pos1[i])) > 5 &&
      as.numeric(gsub("LG ", "", intra_cons.links$chr1[i])) != as.numeric(gsub("PS_chr", "", intra_cons.links$ann1[i])) &&
      intra_cons.links$ann1[i] == intra_cons.links$ann2[i]) {
    col.intra.links[i] <- "red"
  }
}

####
#> unique(col.intra.links)
#[1] "grey90" "grey35" "red"  
####

col.inter.links <- rep("black", nrow(inter_cons.links))  ## default to black

for (i in 1:nrow(inter_cons.links)) {
  
  ann1_is_chr <- grepl("^PS_chr", inter_cons.links$ann1[i])
  ann2_is_chr <- grepl("^PS_chr", inter_cons.links$ann2[i])
  
  if (ann1_is_chr & ann2_is_chr) {
    chr1_num <- as.numeric(gsub("LG ", "", inter_cons.links$chr1[i]))
    ann1_num <- as.numeric(gsub("PS_chr", "", inter_cons.links$ann1[i]))
    chr2_num <- as.numeric(gsub("LG ", "", inter_cons.links$chr2[i]))
    ann2_num <- as.numeric(gsub("PS_chr", "", inter_cons.links$ann2[i]))
    
    if (chr1_num == ann1_num | chr2_num == ann2_num) {
      col.inter.links[i] <- "orange"
    } else {
      col.inter.links[i] <- "darkblue"
    }
  }
}

####
#> unique(col.inter.links)
#[1] "orange"   "darkblue" "black"
####

intra_cons.links.fixed <- intra_cons.links[,c("chr1","pos1","marker1","chr2","pos2","marker2")]
inter_cons.links.fixed <- inter_cons.links[,c("chr1","pos1","marker1","chr2","pos2","marker2")]
# the two lines above are for formatting purposes; it's what the circos function expects

inter_cons.links2<-inter_cons.links.fixed[which(col.inter.links == "darkblue"),]
intra_cons.links2<-intra_cons.links.fixed[which(col.intra.links== "grey35"),]
intra_cons.links3<-intra_cons.links.fixed[which(col.intra.links== "red"),]

#pdf("circos_pine.pdf")
png("circos_pine.png", width=4000, height=4000, res=600)#, bg = "transparent")
par(mar=c(2,2,2,2),mfrow=c(1,1))
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
circos(R=350, cir=db, type="chr",col="grey50", print.chr.lab=FALSE, W=4, scale=FALSE)
circos(R=310, cir=db, mapping=cons.v, type="b3",col="black", W=40 ,lwd=1,B=TRUE,scale=FALSE)
circos(R=305, cir=db, mapping=intra_cons.links.fixed, type="link2",B=TRUE, col=col.intra.links)
circos(R=305, cir=db, mapping=intra_cons.links2, type="link2", col="grey35")
circos(R=305, cir=db, mapping=intra_cons.links3, type="link2", col="red")
circos(R=185, cir=db, mapping=inter_cons.links.fixed, type="link",B=TRUE, col=col.inter.links)
circos(R=185, cir=db, mapping=inter_cons.links2, type="link", col="dark blue")
text(395,730,"A")
text(395,650,"B")
text(395,550,"C")
text(540,750,"LG I")
text(700,640,"LG II")
text(785,480,"LG III")
text(765,255,"LG IV")
text(645,100,"LG V")
text(470,30,"LG VI")
text(280,40,"LG VII")
text(100,150,"LG VIII")
text(20,330,"LG IX")
text(30,500,"LG X")
text(110,650,"LG XI")
text(250,750,"LG XII")
dev.off()
