library('rioja')
library('data.table')
library('rioja')
options(scipen=20)
bin_ord_fh<-"Normal_rep2/raw/100000/Normal_rep2_100000_abs.bed"
re1<-100000
raw_fh <- "Normal_rep2/iced/100000/Normal_rep2_100000_iced.matrix"
all_valid <- "../data/Normal_rep2/Normal_rep2.allValidPairs"
out_prefix <- "Normal_rep2"
csoretool <- 'CscoreTool1.1'
out_folder <- 'local_compartment'
cmd <- paste0('awk  \'BEGIN{OFS="\\t"}($2=="AA" && $3 >BB && $5=="AA" && $6<=CC) ',
              '{print $1,"chr1",$3-BB,$4,"chr1",$6-BB,$7}\' DD > EE')
get_chpos <- function(bin_ord_fh){
  dat <- read.table(bin_ord_fh,
                    colClasses = rep(c('character','integer'),c(1,3)))
  tmp <- tapply(dat[,4],dat[,1],max)
  chpos <- data.frame(ch=names(tmp),len=tmp,stringsAsFactors = F)
  chpos <- chpos[order(chpos$len),]
  return(chpos)
}
read_inter <- function(fh,chpos1){
  # read interaction matrix (HiC-Pro format)
  # bin1 bin2 IF
  dat <- fread(paste(fh))
  colnames(dat) <- c('st','ed','IF')
  dat$ch1 <- as.integer(cut(dat$st,c(0,chpos1[,2])))
  dat$ch2 <- as.integer(cut(dat$ed,c(0,chpos1[,2])))
  dat <- dat[st!=ed,]
  return(dat)
}
onnorm <- function(dat){
  #normalize cis interaction by genomic distance
  data <- dat[,.(M=mean(IF)),keyby=.(d=ed-st)]
  datb <- data.table(M=rep(0,max(data$d)))
  datb$M[data$d] <- data$M
  dat$IF <- dat$IF/datb$M[dat[,ed-st]]
  return(dat)
}
cuthl <- function(vc,h=F,l=F){
  # change the range of a matrix or vector
  # the value greater than h will be assigned as h
  # and value less than l will be assigned as l
  if(h){
    vc[vc>h] <- h
  }
  if(l){
    vc[vc<l] <- l
  }
  return(vc)
}
cuthl_q <- function(vc,h_q=F,l_q=F){
  #
  if(h_q){
    h <- quantile(vc,h_q,na.rm=T)
    vc <- cuthl(vc,h=h)
  }
  if(l_q){
    l <- quantile(vc,l_q,na.rm=T)
    vc <- cuthl(vc,l=l)
  }
  return(vc)
}
int2matf <- function(dat,pos1,pos2,h_q=0.975,l_q=F){
  dat <- dat[st>pos1 & ed<=pos1+pos2,]
  dat$IF <- cuthl_q(dat$IF,h_q,l_q)
  mat <- array(0,dim=rep(pos2,2))
  mat[as.matrix(dat[,1:2])-pos1] <- dat$IF
  mat[as.matrix(dat[,2:1])-pos1] <- dat$IF
  return(mat)
}
partitionf <- function(mat,mnum=250,extend=50){
  d <- dist(mat)
  ch1 <- chclust(d)
  tmp <- bstick(ch1,20,plot=F)
  md <- which(tmp[,2]-tmp[,3]< max(c(tmp[,3]),tmp[,2]/10))[1]
  part=c(0,cumsum(table(cutree(ch1,md))))
  return(part)
}
cmdf <- function(cmd,ch,st,ed,all_valid,out,nch){
  cmd1 <- gsub('AA',ch,cmd)
  cmd1 <- gsub('BB',st,cmd1)
  cmd1 <- gsub('CC',ed,cmd1)
  cmd1 <- gsub('DD',all_valid,cmd1)
  cmd1 <- gsub('EE',out,cmd1)
  return(cmd1)
}
if(!dir.exists(out_folder)){
  dir.create(out_folder)
}
chpos1 <- get_chpos(bin_ord_fh)
dat <- read_inter(raw_fh,chpos1)
dat <- dat[ch1==ch2]
datn <- onnorm(dat)
pos1 <- c(0,chpos1[,2])
parts <- numeric()
res <- numeric()
for(n in 1:nrow(chpos1)){
  ch <- chpos1[n,1]
  clen=pos1[n+1]-pos1[n]
  matn <- int2matf(datn,pos1[n],clen)
  part <- partitionf(matn)
  
  for(i in 2:length(part)){
    st <- part[i-1]
    ed <- part[i]
    bed1 <- data.frame('chr1',seq(1,ed-st)*re1-re1+1,seq(ed-st)*re1)
    out1 <- paste0(out_folder,'/',out_prefix,'_',ch,'_',i-1,'.bed')
    out2 <- paste0(out_folder,'/',out_prefix,'_',ch,'_',i-1,'.summary')
    out3 <- paste0(out_folder,'/',out_prefix,'_',ch,'_',i-1)
    write.table(bed1,out1,sep='\t',quote=F,col.names=F,row.names=F)
    cmd1 <- cmdf(cmd,ch,st*re1,ed*re1,all_valid,out2)
    system(cmd1)
    cmd2 <- paste(csoretool,out1,out2,out3,'4',re1*2,'chr1')
    system(cmd2)
    l_cpt <- read.table(paste0(out3,'_cscore.txt'))
    res <- rbind(res,data.frame(ch,bed1[,2:3]+st*re1,i-1,l_cpt[,2]))
  }
  parts <- rbind(parts,data.frame(ch,1:length(part[-1]),
                                  part[-length(part)],part[-1]))
  
}
colnames(parts) <- c('chr','blocks','st','ed')
colnames(res) <- c('chr','st','ed','block','local_compartment')
write.table(parts,paste0(out_prefix,'_',re1,'_parts.txt'),
            row.names=F,sep='\t',quote=F)
write.table(res,paste0(out_prefix,'_',re1,'_local_compartment.txt'),
            row.names=F,sep='\t',quote=F)
