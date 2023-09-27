setwd("/Users/uym2/my_gits/spalin_intMemoir/data_for_slides/Inferring_Ancestral_Locations")

d = read.table("combined_data.txt",header=T)

d$method = factor(d$method,levels=c("true_brlen","problin","spalin","random_brlen"),
                  labels=c("true","problin","spalin","random"))
ggplot(d[d$method %in% c("problin","spalin","true","random"),],aes(x=method,y=error)) + 
  #geom_boxplot() + 
  stat_summary() + theme_classic()
ggsave("acestral_locations.pdf",width = 4,height=4)
