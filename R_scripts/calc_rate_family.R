#!/bin/env Rscript
library(dplyr)

specie <- commandArgs(TRUE)[1]
child <- commandArgs(TRUE)[2]
out_dir <- commandArgs(TRUE)[3]
dat_dir <- "/project/MutationRates/NewDeNovoCalling/dat_files/gatk/"
dat_dir2 <- file.path(dat_dir,specie)
filename.homref_parents <- file.path(dat_dir2,paste("homoref_test_SNV_w_DPS20_w_lowMQ_w_pc_w_repeat_w_sd.dat", sep=''))
filename.homref_vs_alt <- file.path(dat_dir2,paste("het_test_SNV_w_DPS20_w_lowMQ_w_pc_w_repeat_w_sd.dat", sep=''))
filename.denovo <- file.path(dat_dir2,paste("denovo_raw_SNV_",child,"_w_DPS20_w_known_w_lowMQ_w_pc_w_repeat_w_sd.dat", sep=''))
callabilityDir <- "/project/MutationRates/NewDeNovoCalling/new_family_coverage_wr/"

Homo.minPL01 <- as.integer(commandArgs(TRUE)[4])
Het.minPL00 <- as.integer(commandArgs(TRUE)[5])
Homo.maxAD2 <- as.integer(commandArgs(TRUE)[6])
Homo.maxAR <- as.numeric(commandArgs(TRUE)[7])
minDP <- as.numeric(commandArgs(TRUE)[8])
maxDP <- as.numeric(commandArgs(TRUE)[9])
max_lowMQ_AD2 <- as.numeric(commandArgs(TRUE)[10])
min_relcov <- as.numeric(commandArgs(TRUE)[11])
max_relcov <- as.numeric(commandArgs(TRUE)[12])
SiteQC.minReadPosRankSum <-  as.numeric(commandArgs(TRUE)[13])
SiteQC.maxReadPosRankSum <- as.numeric(commandArgs(TRUE)[14])

MAX.OAD <- 1000000
SiteQC.minMQ <- 0
SiteQC.minFS <- 30
SiteQC.minBaseQRankSum <- -13
SiteQC.maxBaseQRankSum <- 13
SiteQC.minMQRankSum <-  -10 
SiteQC.maxMQRankSum <- 10
Het.minStrandCount <- 1
minMQ=20

add_LGQ <- function(d){
  d$CHILD.rc <- d$CHILD.DP/d$CHILD.COV
  d$FATHER.rc <- d$FATHER.DP/d$FATHER.COV
  d$MOTHER.rc <- d$MOTHER.DP/d$MOTHER.COV
  d$CHILD.minAltSAC <- pmin(d$CHILD.SAC3,d$CHILD.SAC4)
  d$FATHER.lowMQ.AD2 <- d$FATHER.AD2.MQ0 + d$FATHER.AD2.MQ1t20 + d$FATHER.AD2.MQ20t39 + d$FATHER.AD2.MQ40p
  d$MOTHER.lowMQ.AD2 <- d$MOTHER.AD2.MQ0 + d$MOTHER.AD2.MQ1t20 + d$MOTHER.AD2.MQ20t39 + d$MOTHER.AD2.MQ40p
  d$CHILD.lowMQ.AD2 <- d$CHILD.AD2.MQ0 + d$CHILD.AD2.MQ1t20 + d$CHILD.AD2.MQ20t39 + d$CHILD.AD2.MQ40p
  d$CHILD.LGQ <- pmax(pmin(d$CHILD.PL00,d$CHILD.PL01),pmin(d$CHILD.PL01,d$CHILD.PL11), pmin(d$CHILD.PL00, d$CHILD.PL11))
  d$FATHER.LGQ <- pmax(pmin(d$FATHER.PL00,d$FATHER.PL01),pmin(d$FATHER.PL01,d$FATHER.PL11), pmin(d$FATHER.PL00, d$FATHER.PL11))
  d$MOTHER.LGQ <- pmax(pmin(d$MOTHER.PL00,d$MOTHER.PL01),pmin(d$MOTHER.PL01,d$MOTHER.PL11), pmin(d$MOTHER.PL00, d$MOTHER.PL11))
  d
}

alpha.ReadPos <- pnorm(SiteQC.minReadPosRankSum)+pnorm(-SiteQC.maxReadPosRankSum)
alpha.MQRS <- pnorm(SiteQC.minMQRankSum)+pnorm(-SiteQC.maxMQRankSum)
alpha.BQRS <- pnorm(SiteQC.minBaseQRankSum)+pnorm(-SiteQC.maxBaseQRankSum)
alpha.FS <- 10^(-SiteQC.minFS/10)
# alpha.site is the expected number of true mutations that are filtered away
# because of the site filters on ReadPosRankSum, MQRankSum, BaseQRankSum and FS.
alpha.site <- 1- ((1-alpha.FS)*(1-alpha.ReadPos)*(1-alpha.MQRS)*(1-alpha.BQRS))

print(alpha.site)

filter_sites <- function(d) {
  d %>% filter(
    CHROM != "chrX",
    FS < SiteQC.minFS,
    ReadPosRankSum >= SiteQC.minReadPosRankSum,
    ReadPosRankSum <= SiteQC.maxReadPosRankSum,
    MQRankSum >= SiteQC.minMQRankSum,
    MQRankSum <= SiteQC.maxMQRankSum,
    MQ >= SiteQC.minMQ,
    !is.na(ReadPosRankSum),
    !is.na(MQRankSum),
    !is.na(BaseQRankSum),
    SD == 'OUT')
}

#### Make callability ####
d <- read.table(filename.homref_parents, head=T)
d <- filter_sites(d)
d <- add_LGQ(d)

print(nrow(d))

callability.HomRef <- d %>%
  filter(FATHER.DPS20 >= minDP,
         MOTHER.DPS20 >= minDP,
         FATHER.DPS20 <= maxDP,
         MOTHER.DPS20 <= maxDP,
         FATHER.rc > min_relcov,
         FATHER.rc < max_relcov,
         MOTHER.rc > min_relcov,
         MOTHER.rc < max_relcov,
         FATHER.PL01 > max(55, Homo.minPL01),
         MOTHER.PL01 > max(55, Homo.minPL01),
         FATHER.AD2 == 0,
         MOTHER.AD2 == 0,
         !is.na(ReadPosRankSum),
         !is.na(MQ),
         !is.na(BaseQRankSum),
         CHILD.AD1+CHILD.AD2 > 0) %>%
  mutate(CHILD.AR=CHILD.AD2/CHILD.DP, DP=CHILD.DPS20) %>%
  group_by(DP) %>%
  summarise(callability.Hom = mean(CHILD.GT == 0 &
                                   CHILD.PL01 >= Homo.minPL01 &
                                   CHILD.AD2 <= Homo.maxAD2 &
                                   CHILD.lowMQ.AD2 <= max_lowMQ_AD2 &
                                   CHILD.AR <= Homo.maxAR &
                                   CHILD.DPS20 <= maxDP &
                                   CHILD.DPS20 >= minDP))



d <- read.table(filename.homref_vs_alt, head=T)
d <- filter_sites(d)
d <- add_LGQ(d)

print(nrow(d))

callability.Het <- d %>%
  filter(FATHER.DPS20 >= minDP,
         MOTHER.DPS20 >= minDP,
         FATHER.DPS20 <= maxDP,
         MOTHER.DPS20 <= maxDP,
         FATHER.rc > min_relcov,
         FATHER.rc < max_relcov,
         MOTHER.rc > min_relcov,
         MOTHER.rc < max_relcov,
         FATHER.PL01 > max(55, Homo.minPL01),
         MOTHER.PL01 > max(55, Homo.minPL01),
         !is.na(ReadPosRankSum),
         !is.na(MQ),
         !is.na(BaseQRankSum),
         CHILD.AD1+CHILD.AD2>0) %>%
  mutate(CHILD.allele.balance=CHILD.AD2/(CHILD.AD1+CHILD.AD2),
         DP = CHILD.DPS20,
         CHILD.minhom=pmin(CHILD.PL00,CHILD.PL11)) %>%
  group_by(DP) %>%
  summarise(callability.Het = mean(CHILD.GT == 1 &
                                   CHILD.PL00 >= Het.minPL00 &
                                   CHILD.minAltSAC > Het.minStrandCount &
                                   CHILD.allele.balance >= 0.3 &
                                   CHILD.DPS20 <= maxDP &
                                   CHILD.DPS20 >= minDP
                                   ))

print(nrow(d))

pdf(file.path(out_dir,"callability.pdf"), height=4, width=8)
par(mfrow=c(1,2))
plot(callability.HomRef,xlim=c(0,100), type="o")
plot(callability.Het,xlim=c(0,100), type="o")
dev.off()

d <- read.table(filename.denovo, head=T)
d <- filter_sites(d)
d <- add_LGQ(d)

callability.Father <- callability.HomRef
callability.Father[callability.Father$DP>d$FATHER.COV[1]*max_relcov,]$callability.Hom <- 0
if (min_relcov>0){
  callability.Father[callability.Father$DP<d$FATHER.COV[1]*min_relcov,]$callability.Hom <- 0
}
callability.Mother <- callability.HomRef
callability.Mother[callability.Mother$DP>d$MOTHER.COV[1]*max_relcov,]$callability.Hom <- 0
if (min_relcov>0){
  callability.Mother[callability.Mother$DP<d$MOTHER.COV[1]*min_relcov,]$callability.Hom <- 0
}
callability.Het[callability.Het$DP>d$CHILD.COV[1]*max_relcov,]$callability.Het <- 0
if (min_relcov>0){
  callability.Het[callability.Het$DP<d$CHILD.COV[1]*min_relcov,]$callability.Het <- 0
}

get_nsites <- function(family, CHILD.callability=callability.Het) {
  filename <- paste(callabilityDir,'/', family ,"_minQ",minMQ,"_minq10_type_autosomes.txt", sep="")
  d.sites <- read.table(filename, head=F)
  names(d.sites) <- c("FATHER.DP","MOTHER.DP","CHILD.DP", "SITETYPE", "REPEAT", "COUNT")
  d.sites <- merge(d.sites, callability.Father, by.x="FATHER.DP", by.y="DP")
  names(d.sites)[ncol(d.sites)] <- "callability.FATHER"
  d.sites <- merge(d.sites, callability.Mother, by.x="MOTHER.DP", by.y="DP")
  names(d.sites)[ncol(d.sites)] <- "callability.MOTHER"
  d.sites <- merge(d.sites, CHILD.callability, by.x="CHILD.DP", by.y="DP")
  names(d.sites)[ncol(d.sites)] <- "callability.CHILD"
  d.sites %>% group_by(SITETYPE,REPEAT) %>% summarise(SNV.SUM=sum(COUNT*callability.FATHER*callability.MOTHER*callability.CHILD))
}

d.tmp <- get_nsites(child)
d.nsites <- data.frame(CHILD=child,
                       CSUM.CpG=sum(d.tmp[d.tmp$SITETYPE=="CpG",]$SNV.SUM),
                       CSUM.nonCpGstrong=sum(d.tmp[d.tmp$SITETYPE=="nonCpG-Strong",]$SNV.SUM),
                       CSUM.weak=sum(d.tmp[d.tmp$SITETYPE=="Weak",]$SNV.SUM),
                       CSUM.repeat=sum(d.tmp[d.tmp$REPEAT==TRUE,]$SNV.SUM),
                       CSUM.nonrepeat=sum(d.tmp[d.tmp$REPEAT==FALSE,]$SNV.SUM))


print(d.nsites)

#### Analyse mutations ####

d <- d %>%
  mutate(
    trio.min.DPS = pmin(FATHER.DPS20, MOTHER.DPS20, CHILD.DPS20),
    trio.max.DPS = pmax(FATHER.DPS20, MOTHER.DPS20, CHILD.DPS20),
    trio.min.rc = pmin(FATHER.rc, MOTHER.rc, CHILD.rc),
    trio.max.rc = pmax(FATHER.rc, MOTHER.rc, CHILD.rc),
    PARENTS.min.GQ = pmin(FATHER.GQ, MOTHER.GQ),
    PARENTS.max.AD2 = pmax(d$FATHER.AD2, d$MOTHER.AD2),
    PARENTS.max.lowMQ.AD2 = pmax(FATHER.lowMQ.AD2, MOTHER.lowMQ.AD2),
    PARENTS.max.AR = pmax(d$FATHER.AD2/d$FATHER.DP, d$MOTHER.AD2/d$MOTHER.DP),
    CHILD.allele.balance = CHILD.AD2/(CHILD.AD1 + CHILD.AD2)
  ) %>%
  filter(
    trio.max.DPS <= maxDP,
    trio.min.DPS >= minDP,
    trio.max.rc < max_relcov,
    trio.min.rc > min_relcov,
    PARENTS.min.GQ >= Homo.minPL01,
    CHILD.PL00 >= Het.minPL00,
    PARENTS.max.AD2 <= Homo.maxAD2,
    PARENTS.max.lowMQ.AD2 <= max_lowMQ_AD2,
    PARENTS.max.AR <= Homo.maxAR,
    MAX.OTHER.AD <= MAX.OAD,
    CHILD.minAltSAC >= Het.minStrandCount
  )

add_context_type <- function(d){
  d$context_type <-  "Weak"
  d[substr(d$CONTEXT,2,2)=='C',]$context_type <- "nonCpG-Strong"
  if (length(grep('CG',d$CONTEXT))>0) {
    d[grep('CG',d$CONTEXT),]$context_type <- "CpG"
  }
  d$cpg_type <- "nonCpG"
  if (length(grep('CG',d$CONTEXT))>0) {
    d[grep('CG',d$CONTEXT),]$cpg_type <- "CpG"
  }
  d
}
d <- add_context_type(d)

# Write a table with all mutations including likely somatic mutations 
write.table(d,file=file.path(out_dir,"all_mutations.dat"),row.names=F, col.names=T, quote=F)

# Remove likely somatic mutations (VAF < 0.3)
d <- d %>%
  filter(CHILD.allele.balance >= 0.3)

# Write a table with the found germline de novo mutations
write.table(d,file=file.path(out_dir,"germline_mutations.dat"),row.names=F, col.names=T, quote=F)

# Calculate summary statistics such as rate, fraction inherited etc.
rate.data <- d %>%
  group_by(CHILD) %>%
  summarise(nmut=sum(!is.na(ts)),
            nCpG=sum(context_type=="CpG"),
            nCpGtrans=sum(context_type=="CpG" & ts==1),
            nNonCpGstrong=sum(context_type=="nonCpG-Strong"),
            nWeak=sum(context_type=="Weak"),
            n.ts=sum(ts==1),
            n.tv=sum(ts==0),
            known.frac=mean(known_variant),
            GC.N=sum(GC.N),
            frac.inherited=sum(GC.AC)/sum(GC.N),
            GC.AC=sum(GC.AC),
            nPat = sum(POO=='P'),
            nMat = sum(POO=='M'),
            nRepeat = sum(repeat.),
            nNonRepeat = sum(repeat.==FALSE)) %>%
  merge(d.nsites) %>%
  as.data.frame()
  
par.info <- cbind(HOMGQ=Homo.minPL01, HETGQ=Het.minPL00, AD2=Homo.maxAD2, AR=Homo.maxAR, SiteQC.minFS,
                  minRPRS=SiteQC.minReadPosRankSum, maxRPRS=SiteQC.maxReadPosRankSum,
                  minBQRS=SiteQC.minBaseQRankSum, maxBQRS=SiteQC.maxBaseQRankSum,
                  minMQRS=SiteQC.minMQRankSum, maxMQRS=SiteQC.maxMQRankSum,
                  minDP=minDP, maxDP=maxDP, alpha.site=alpha.site, minMQ=SiteQC.minMQ,
                  lowMQAD2=max_lowMQ_AD2, minrc=min_relcov, maxrc=max_relcov)

# Write the summary statistics and the used filtering parameters.
write.table(
  cbind(rate.data, par.info),
  file=file.path(out_dir,"avg_rate.txt"),
  quote=FALSE,
  col.names=TRUE,
  row.names=FALSE)
