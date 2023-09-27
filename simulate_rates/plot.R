setwd("/Users/uym2/my_gits/spalin_intMemoir/simulate_rates")
require(ggplot2)
require(reshape2)

a_start = 0.008
a_end = 0.004
n = 215

la_start = log(a_start)
la_end = log(a_end)
la = seq(la_start,la_end,(la_end-la_start)/(n-1))
a = exp(la)

a2_start = 0.01
a2_end = 0.002
n = 215

la2_start = log(a2_start)
la2_end = log(a2_end)
la2 = seq(la2_start,la2_end,(la2_end-la2_start)/(n-1))
a2 = exp(la2)

d = data.frame(frame=seq(1:n)-1,exp=a, exp2=a2,
               linear=seq(0.008,0.004,(0.004-0.008)/(n-1)),
               constant = rep(0.006,n))
d_true = melt(d,id.vars = c("frame"))
d_true$type = "true"
colnames(d_true) = c("frame","rateModel","rate","type")

d_est = read.table("problin_k10_inferred_rates_smooth_2500.txt",header=T)
d_est$type = "est"
d_est$rate = d_est$rate/1000
d = rbind(d_true,d_est)

d$type = factor(d$type,levels = c("est","true"),labels = c("estimated","true"))
d$rateModel = factor(d$rateModel,levels=c("constant","linear","exp2","exp"),labels=c("constant","linear","exponential","exp1"))

ggplot(d[d$rateModel != "exp1",],aes(x=frame,y=rate,linetype=type)) + geom_line() +
  facet_wrap(~rateModel)  + 
  theme_classic()
ggsave("problin_k20_inferred_rates_smooth_500.pdf")
