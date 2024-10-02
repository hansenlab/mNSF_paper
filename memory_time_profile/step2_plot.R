
library(ggplot2)
library(dplyr)

######
dir_processedData = "//dcs04/hansen/data/ywang/ST/memory_time_profile/"
setwd(dir_processedData)

# nfactor - run 3 times
df_nfactor = read.csv("memoryUsage_mouseSagittal_n3_iterL_runSeperately.csv")
df_nfactor$memory_peak = df_nfactor$peak_after_training
df_nfactor = df_nfactor %>%
  group_by(L) %>%
  summarise(std_memory = sd(memory_peak),memory_peak = mean(memory_peak), 
            std_runtime = sd(runtime),
            runtime = mean(runtime),n = n())

#ngene - run 3 times
df_ngene = read.csv("memoryUsage_mouseSagittal_iterG_runSeperately.csv")
df_ngene$memory_peak = df_ngene$peak_after_training
df_ngene = df_ngene %>%
  group_by(ngene) %>%
  summarise(std_memory = sd(memory_peak),memory_peak = mean(memory_peak), std_runtime = sd(runtime),
            runtime = mean(runtime),n = n())# df_ngene$memory_peak = df_ngene$peak_after_training

#nsample - run 3 times
df_nsample = read.csv("memoryUsage_mouseSagittal_L3_iterM_runSeperately.csv")
df_nsample$memory_peak = df_nsample$peak_after_training
df_nsample = df_nsample %>%
  group_by(M) %>%
  summarise(std_memory = sd(memory_peak),memory_peak = mean(memory_peak), std_runtime = sd(runtime),
            runtime = mean(runtime),n = n())
#nspot_perSample - run 3 times
df_nspot = read.csv("memoryUsage_mouseSagittal_iterNspot_runSeperately.csv")
df_nspot$memory_peak = df_nspot$peak_after_training
df_nspot = df_nspot %>%
  group_by(nspot) %>%
  summarise(std_memory = sd(memory_peak),memory_peak = mean(memory_peak), std_runtime = sd(runtime),
            runtime = mean(runtime),n = n())
#######  make plot
# nfactor
pdf("profile_nfactor_memory.pdf", height = 2, width = 5)
ggplot(df_nfactor, aes(y=memory_peak, x= L))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("memory ~ number of factors")+
  xlab("number of factors")+
  ylab("memory (GB)")+
 
  geom_ribbon(aes(ymin = memory_peak - 1.96*std_memory,
                  ymax = memory_peak + 1.96*std_memory),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)     +geom_point(size=2)  
dev.off()

pdf("profile_nfactor_time.pdf", height = 2, width = 5)
ggplot(df_nfactor, aes(y=runtime, x= L))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
   # +geom_point(size=2)    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("runtime ~ number of factors")+
  xlab("number of factors")+
  ylab("time (second)")+
  
  geom_ribbon(aes(ymin = runtime - 1.96*std_runtime,
                  ymax = runtime + 1.96*std_runtime),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)     +geom_point(size=2)  
dev.off()


#ngene
pdf("profile_ngene_memory.pdf", height = 2, width = 5)
ggplot(df_ngene, aes(y=memory_peak, x= ngene))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
   # +geom_point(size=2)   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("memory ~ number of genes")+
  xlab("number of genes")+
  ylab("memory (GB)")+
  
  geom_ribbon(aes(ymin = memory_peak - 1.96*std_memory,
                  ymax = memory_peak + 1.96*std_memory),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)   +geom_point(size=2)  
dev.off()

pdf("profile_ngene_time.pdf", height = 2, width = 5)
ggplot(df_ngene, aes(y=runtime, x= ngene))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
   # +geom_point(size=2)    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("runtime ~ number of genes")+
  xlab("number of genes")+
  ylab("time (second)")+
  
  geom_ribbon(aes(ymin = runtime - 1.96*std_runtime,
                  ymax = runtime + 1.96*std_runtime),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)    +geom_point(size=2)  
dev.off()

#nsample
pdf("profile_nsample_memory.pdf", height = 2, width = 5)
ggplot(df_nsample, aes(y=memory_peak, x= M))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
   # +geom_point(size=2)   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("memory ~ number of samples")+
  xlab("number of samples")+
  ylab("memory (GB)")+
  
  geom_ribbon(aes(ymin = memory_peak - 1.96*std_memory,
                  ymax = memory_peak + 1.96*std_memory),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)   +geom_point(size=2)  
dev.off()

pdf("profile_nsample_time.pdf", height = 2, width = 5)
ggplot(df_nsample, aes(y=runtime, x= M))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
   # +geom_point(size=2)    
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("runtime ~ number of samples")+
  xlab("number of samples")+
  ylab("time (second)")+
  
  geom_ribbon(aes(ymin = runtime - 1.96*std_runtime,
                  ymax = runtime + 1.96*std_runtime),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)    +geom_point(size=2)  
dev.off()

# nspot_perSample
pdf("profile_nspot_memory.pdf", height = 2, width = 5)
ggplot(df_nspot, aes(y=memory_peak, x= nspot))+ #+ facet_wrap(~donor,ncol=3)+
  # geom_violin(position=position_dodge(1))+
   # +geom_point(size=2)   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("memory ~ number of spots within each sample")+
  xlab("number of spots within each sample")+
  ylab("memory (GB)")+
  
  geom_ribbon(aes(ymin = memory_peak - 1.96*std_memory,
                  ymax = memory_peak + 1.96*std_memory),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)   +geom_point(size=2)  
dev.off()

pdf("profile_nspot_time.pdf", height = 2, width = 5)
ggplot(df_nspot, aes(y=runtime, x= nspot))+ #+ facet_wrap(~donor,ncol=3)+
  # df_nsample(position=position_dodge(1))+
   # +geom_point(size=2)   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("runtime ~ number of spots within each sample")+
  xlab("number of spots within each sample")+
  ylab("time (second)")   +
  
  geom_ribbon(aes(ymin = runtime - 1.96*std_runtime,
                  ymax = runtime + 1.96*std_runtime),    # shadowing cnf intervals
              fill = "steelblue2") + 
  geom_line(color = "firebrick",
            size = 1)              +geom_point(size=2)     #
dev.off()