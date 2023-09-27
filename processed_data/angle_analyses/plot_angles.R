setwd("/Users/uym2/my_gits/spalin_intMemoir/processed_data/angle_analyses")

d = read.table("slope_and_theta.csv",header=T)

ggplot(d,aes(x=slope)) + geom_histogram(color="black",bins = 180) + facet_wrap(~sample)

ggplot(d,aes(x=slope)) + geom_histogram(color="black",bins = 180)


ggplot(d,aes(x=theta)) + geom_histogram(color="black",bins = 180) + 
  geom_vline(xintercept = 90,color="red") + 
  facet_wrap(~sample,nrow = 3) + theme_classic()
ggsave("angle_dividing_frame.pdf")

ggplot(d,aes(x=theta)) + geom_density() #+ facet_wrap(~sample)

# Gaussian model
x1 = rnorm(1000000,mean=0,sd=1.5)
y1 = rnorm(1000000,mean=0,sd=0.05)

x2 = rnorm(1000000,mean=0,sd=1.5)
y2 = rnorm(1000000,mean=0,sd=0.05)

cos_theta_1 = (x1*x2+y1*y2)/sqrt(x1*x1+y1*y1)/sqrt(x2*x2+y2*y2)
theta_1 = acos(cos_theta_1)/pi*180

x1 = rnorm(1000000,mean=0,sd=5.5)
y1 = rnorm(1000000,mean=0,sd=5.5)

x2 = rnorm(1000000,mean=0,sd=5.5)
y2 = rnorm(1000000,mean=0,sd=5.5)

cos_theta_2 = (x1*x2+y1*y2)/sqrt(x1*x1+y1*y1)/sqrt(x2*x2+y2*y2)
theta_2 = acos(cos_theta_2)/pi*180


d_gauss_1 = data.frame(theta = theta_1)
d_gauss_1$name = "theta_1"
d_gauss_2 = data.frame(theta = theta_2)
d_gauss_2$name = "theta_2"
d_gauss = rbind(d_gauss_1,d_gauss_2)

ggplot(d_gauss,aes(x=theta,color=name)) + 
  #geom_histogram(color="black",bins = 180) + 
  geom_density() + 
  #geom_vline(xintercept = 90,color="red") + 
  #facet_wrap(~sample,nrow = 3) + 
  theme_classic()

mean(theta)
ggsave("angle_gaussian.pdf")
