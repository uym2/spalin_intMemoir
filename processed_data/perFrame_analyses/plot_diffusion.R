setwd("/Users/uym2/my_gits/spalin_intMemoir/processed_data/perFrame_analyses")
require(ggplot2)

d = read.table("all_diffussion.txt",header=T)

# x-coordinate histogram => zero inflated 
ggplot(d,aes(x=par_x-x)) + 
  geom_histogram(color="black") + theme_classic() + 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  facet_wrap(~ID,scale="free")
ggsave("zero_inflated_x.pdf",width = 15,height=15)

# y-coordinate histogram => zero inflated 
ggplot(d,aes(x=par_y-y)) + 
  geom_histogram(color="black") + theme_classic() + 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  facet_wrap(~ID,scale="free")
ggsave("zero_inflated_y.pdf",width = 15,height=15)

# x coordinate: qq-plot without zeros
# include "division frames"
ggplot(d[(d$par_x-d$x) != 0,],aes(sample=par_x-x)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic() + facet_wrap(~ID,scale="free")
ggsave("qqplot_x_all.pdf",width=15,height=15)

# exclude "division frames"
ggplot(d[d$div_frame == "False" & (d$par_x-d$x) != 0,],aes(sample=par_x-x)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic() + facet_wrap(~ID,scale="free")
ggsave("qqplot_x_exclude_divframes.pdf",width=15,height=15)

# y coordinate: qq-plot without zeros
# include "division frames"
ggplot(d[(d$par_y-d$y) != 0,],aes(sample=par_y-y)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic() + facet_wrap(~ID,scale="free")
ggsave("qqplot_y_all.pdf",width=15,height=15)

# exclude "division frames"
ggplot(d[d$div_frame == "False" & (d$par_y-d$y) != 0,],aes(sample=par_y-y)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic() + facet_wrap(~ID,scale="free")
ggsave("qqplot_y_exclude_divframes.pdf",width=15,height=15)

#########################
d1 = d[d$ID == "s2c2",]
b = d1[d1$div_frame == "False" & d1$par_y-d1$y != 0 & d1$par_x-d1$x != 0,] # the "background" distribution
sigma_squared = with(b,(var(par_x-x)+var(par_y-y))/2)
mu_prime = with(d1[d1$div_frame == "True",],mean((par_x-x)**2+(par_y-y)**2))
delta = sqrt(mu_prime - 2*sigma_squared)

b$eps_x = b$x-b$par_x + runif(length(b$x))
b$eps_y = b$y-b$par_y + runif(length(b$y))
c = d1[d1$div_frame == "True",] # the dividing frames
c$eps_x = (1-delta/sqrt((c$x-c$par_x)**2+(c$y-c$par_y)**2))*(c$x-c$par_x)
c$eps_y = (1-delta/sqrt((c$x-c$par_x)**2+(c$y-c$par_y)**2))*(c$y-c$par_y)

d2 = rbind(b,c)

ggplot(d2,aes(sample=eps_x)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic() 

ggplot(d2,aes(sample=eps_y)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic() 

###############################
ggplot(d2,aes(sample=y-par_y)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic()

ggplot(d1[d1$par_x-d1$x != 0 & d1$div_frame == "False",],aes(sample=par_x-x)) + 
  geom_qq(size=0.1) + geom_qq_line() + 
  theme_classic()



