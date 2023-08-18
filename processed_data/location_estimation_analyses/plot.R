setwd("/Users/uym2/my_gits/spalin_intMemoir/processed_data/location_estimation_analyses")
require(ggplot2)

d = read.table("estimated_location_combined.txt",header=T)
d$error = sqrt((d$true_x-d$estimated_x)**2+(d$true_y-d$estimated_y)**2)

ggplot(d,aes(y=error,x=model,color=model)) + 
  #geom_boxplot() + 
  stat_summary(position=position_dodge(width=1)) + 
  theme_classic()
ggsave("gauss_vs_t.pdf")

############ permutation test #########
d = read.table("permutation_test.txt",header=T)
ggplot(d,aes(x=type,y=error)) + geom_boxplot()

ggplot(d,aes(y=error,x=as.factor(b2root),fill=type)) + 
  geom_boxplot() + xlab("branch distance to root") + ylab("location error") + 
  theme_classic()
ggsave("error_vs_b2root.pdf")

ggplot(d,aes(y=error,x=b2root,color=type)) + 
  stat_summary() + geom_line(stat = "summary") + 
  xlab("branch distance to root") + ylab("location error") + 
  theme_classic()
ggsave("error_vs_b2root_lineplot.pdf")

############ coestimation with branch lengths #########
d = read.table("inferred_locations.txt",header=T)
d$type = factor(d$type,levels=c("estimated","brlen_coestimated","permutated"),
                labels=c("branch lengths known","branch lengths coestimated","random"))
ggplot(d,aes(x=type,y=error)) + 
  geom_boxplot() + stat_summary() +
  theme_classic() + theme(axis.title.x=element_blank()) 
ggsave("location_estimation.pdf")


