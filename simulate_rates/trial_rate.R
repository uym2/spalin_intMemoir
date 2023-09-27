setwd("/Users/uym2/my_gits/spalin_intMemoir/simulated_data")

d = read.table("k10_r1_brlens.txt",header=T)

ggplot(d[d$trueBrlen>0,],aes(x=method,y=abs(estBrlen-trueBrlen)/max(d$trueBrlen),color=method)) + geom_boxplot()


d = read.table("temp")
ggplot(d,aes(x=V1,y=V2)) + geom_line() + ylim(0.004,0.008)


a_start = 0.01
a_end = 0.002
n = 215

la_start = log(a_start)
la_end = log(a_end)
la = seq(la_start,la_end,(la_end-la_start)/(n-1))
a = exp(la)
d = data.frame(frame=seq(1:n)-1,rate_exp=a,rate_lin=seq(0.008,0.004,(0.004-0.008)/(n-1)))
dm = melt(d,id.vars = c("frame"))

ggplot(dm,aes(x=frame,y=value,color=variable)) + geom_line() + theme_classic()
