setwd("/Users/uym2/my_gits/spalin/perBranch_analyses/")
require(ggplot2)

d = read.table("tree_diffusion.txt",header=T)

ggplot(d,aes(x=dx_norm)) + geom_histogram(color="black") + theme_classic()
ggsave("histogram_dx.pdf")
ggplot(d,aes(sample=dx_norm)) + geom_qq() + geom_qq_line() + theme_classic()
ggsave("qqplot_dx.pdf")

ggplot(d,aes(x=dy_norm)) + geom_histogram(color="black") + theme_classic()
ggsave("histogram_dy.pdf")
ggplot(d,aes(sample=dy_norm)) + geom_qq() + geom_qq_line() + theme_classic()
ggsave("qqplot_dy.pdf")
