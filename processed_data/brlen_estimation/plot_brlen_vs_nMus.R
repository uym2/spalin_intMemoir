setwd("/Users/uym2/my_gits/spalin_intMemoir/processed_data/brlen_estimation")
require(ggplot2)


d = read.table("brlen_vs_nMus.txt",header=T)

ggplot(d,aes(x=trueTime,y=dvalue)) + geom_point(alpha=0.3) + stat_smooth() + facet_wrap(~dtype,scale="free")
