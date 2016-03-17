sp<-c("Agambiae", "Amelifera", "Athaliana", "Bdistachyon", "Bmori", "Cbriggsae", "Ccornix", "Clanatus", "Clupus","Crubella", "Csativus", "Dmelanogaster", "Dsimulans", "Dyakuba", "Falbicolis", "Gaculeatus", "Gmax", "Hmelpomene_2", "Mgallopavo", "Mtruncaluta","Ncrassa", "Olatipes", "Orufipogon", "Ptrichocarpa", "Sitalica", "Sscrofa", "Zmays")

filters<-list()

for (species in sp) {
  path<-paste0("filter_stats/POLYMORPHISM/", species, "/filter.stats", collapse="")
  filters[[species]]<-read.table(path, sep="\t", header=F, stringsAsFactors=F)
  filters[[species]]$sp=species
}

all.filt<-do.call("rbind", filters)
names(all.filt)=c("set", "site", "filter", "count", "sp")

library(plyr)
totals<-ddply(all.filt, .(set,site,sp), summarize, ct=sum(as.numeric(count)))

filt.stats<-merge(all.filt, totals)
filt.stats$prop=filt.stats$count/filt.stats$ct

subset(filt.stats, filter=="noCoverage", select=c("set", "site", "sp", "prop"))

write.table(filt.stats, file="filterings_stats.tsv", sep="\t", quote=F)
