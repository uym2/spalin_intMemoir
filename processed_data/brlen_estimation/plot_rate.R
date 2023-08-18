setwd("/Users/uym2/my_gits/spalin_intMemoir/processed_data/brlen_estimation")

d = read.table("brlen.txt", header=T)
require(ggplot2)
require(chngpt)

d1 = d[d$nodeName != "virtual" & d$nodeAge >0 & d$nodeAge < 216 & d$trueTime > 20 & d$trueTime < 50,] 
       #& d$rep %in% c("s13c1","s2c1","s21c4","s5c7","s8c6","s15c7","s7c3","s5c12","s9c3","s5c3","s5c10","s2c6","s14c5","s8c5"),]

d1$mu = d1$estBrlen/d1$trueTime
#fit=chngptm(formula.1=estBrlen/trueTime~1, formula.2=~nodeAge, data=d1, type="stegmented",family="gaussian")

d1$method = factor(d1$method,levels = c("without_location","with_all_location","with_leaf_location"))
ggplot(d1[d1$method != "with_leaf_location",],
       aes(x=nodeAge,y=mu)) + 
  #geom_point(alpha=0.3) + 
  stat_summary() + 
  #geom_line(stat="summary") + 
  geom_smooth(method="lm") +
  #geom_smooth(aes(group=group),color="red")+
  facet_wrap(~method)+ 
  xlab("Time frame") + ylab("Mutation rate") +
  theme_classic()
ggsave("brlen.pdf")

ggplot(d1[d1$method != "with_leaf_location",],aes(x=trueTime,y=estBrlen)) + 
  #stat_summary() + 
  geom_point(alpha=0.3) + 
  #geom_smooth(method="lm") +
  #geom_smooth(aes(group=group),color="red")+
  facet_wrap(~method)+ 
  xlab("Time unit") + ylab("Mutation unit") +
  theme_classic()
ggsave("brlen_vs_time.pdf")

s10c3 #leaves: 20
s12c4 #leaves: 20
s15c8 #leaves: 20
s16c6 #leaves: 20
s17c1 #leaves: 20
s16c1 #leaves: 21
s18c2 #leaves: 21
s6c3 #leaves: 22
s9c4 #leaves: 22
s8c5 #leaves: 23
s14c5 #leaves: 25
s2c6 #leaves: 26
s5c10 #leaves: 26
s5c3 #leaves: 26
s9c3 #leaves: 26
s5c12 #leaves: 27
s7c3 #leaves: 27
s15c7 #leaves: 29
s8c6 #leaves: 29
s5c7 #leaves: 31
s21c4 #leaves: 32
s2c1 #leaves: 34
s13c1 #leaves: 39