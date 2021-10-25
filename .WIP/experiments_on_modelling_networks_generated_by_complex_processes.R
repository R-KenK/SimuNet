# Load packages and custom functions ----
source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")

# Case study ----
## Generating main variables ----
set.seed(42)

n.rep   <- 7L    # number of weighted networks produced
n.group <- 10L    # number of different groups of individuals generated
n.each  <- 50L    # number of repetition for a given group

## Simulations ----
### Exporting objects to parallel workers ----
cl <- snow::makeCluster(7)
snow::clusterCall(
  cl,function() {
    source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")
  }
)

### Running the simulations ----
param.n <- c(5,8,10,15,20)
param.samp.eff <- c(10,25,50,75,100,150,250,500)

param.list <-
  CJ(n = param.n,
     samp.eff = param.samp.eff,
     n.rep = n.rep,
     group.number = 1:n.group,
     group.rep = 1:n.each)

param.list <-
  param.list[,by = .(n,samp.eff,n.rep,group.number),
             .(group =
                 list(
                   data.table(
                     name = as.character(1:n),
                     value = runif(n),
                     rank = sample(1:n,n)
                   )
                 )
             )
  ] |>
  merge.data.table(x = param.list,by = c("n","samp.eff","n.rep","group.number"))
param.list

start.time <- Sys.time()
start.time
results.list <-
  param.list |>
  subset(n %in% c(5) & samp.eff %in% c(50) & group.number <= 3 & group.rep <= 4) |>
  run_simulations(verbose = FALSE)
end.time <- Sys.time()
end.time
end.time - start.time

snow::stopCluster(cl)

results <- aggregate_results(param.list[c(1,3,5,10),],results.list)
results

# results.multi <- readRDS(".WIP/simulation.data/results.n_10.samp_10.to.500.rds")
# results.multi

## Plotting results ----
x.labs <- c(
  "Expected null distance",
  paste0("... ",c("the real network",
                  "another network (similar generation)",
                  "an Erdős–Rényi graph (p = 0.5)",
                  "a Random network (fixed prob.)",
                  "a Random network (variable prob.)"))
)

x.labs.face <- c(
  "bold",
  rep("plain",length(x.labs) - 1)
)

x.fill <- c("tomato","royalblue",gg_color_hue(5)[2:5])

KS.stat <-
  results |>
  subset(samp.eff == 50) |>
  pull(edge.comp.dt) |>
  rbindlist() |>
  ggplot(aes(y,KS.stat,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  xlab("")+ylab("Kolmogorov-Smirnoff statistic")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))

KS.p <-
  results |>
  subset(samp.eff == 50) |>
  pull(edge.comp.dt) |>
  rbindlist() |>
  ggplot(aes(y,KS.p,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  xlab("")+ylab("Kolmogorov-Smirnoff p-value")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))

KL <-
  results |>
  subset(samp.eff == 50) |>
  pull(edge.comp.dt) |>
  rbindlist() |>
  ggplot(aes(y,KL,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  xlab("")+ylab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))

JS <-
  results |>
  subset(samp.eff == 50) |>
  pull(edge.comp.dt) |>
  rbindlist() |>
  ggplot(aes(y,JS,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  xlab("")+ylab("Jensen-Shannon distance")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))

EMD.e <-
  results |>
  subset(samp.eff == 50) |>
  pull(edge.comp.dt) |>
  rbindlist() |>
  ggplot(aes(y,EMD.e,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  xlab("")+ylab("Earth Mover Distance (Euclidean dist.)")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))

EMD.m <-
  results |>
  subset(samp.eff == 50) |>
  pull(edge.comp.dt) |>
  rbindlist() |>
  ggplot(aes(y,EMD.m,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  xlab("")+ylab("Earth Mover Distance (Manhattan dist.)")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...", gp=grid::gpar(fontsize = 16,fontface = "bold")))

# Graphic exploration ---------------------------------------------------------------------------------------
## Data handling and wrangling ----
type.str <- c("SimuNet network",
              "the real network",
              "another network (similar generation)",
              "an Erdős–Rényi graph (p = 0.5)",
              "a Random network (fixed prob.)",
              "a Random network (variable prob.)")
Adj.obs.dt <-
  expand.grid(i = 1:n,j = 1:n) |>
  subset(i < j) |>
  as.matrix() %>%
  {
    cbind(
      .,
      weight = results$real.dist[,,631][.]
    )
  } |>
  data.table()

real.dt <-
  results$edge.dt |>
  subset(samp.eff == 50 & group.rep == 1) |>
  subset(type == "real") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

SimuNet.dt <-
  results$edge.dt |>
  subset(samp.eff == 50 & group.rep == 1) |>
  subset(type == "SimuNet") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

## plots ----
results$edge.dt |>
  subset(samp.eff == 50 & group.rep == 1) |>
  ggplot(aes(weight, colour = type, fill = type))+
  facet_wrap(. ~ paste0("edge ",i,"-",j),nrow = n,scales = "free_y")+
  geom_histogram(binwidth = 1,alpha = 0.15,position = "identity",colour = NA)+
  geom_histogram(data = results$edge.dt[type %in% c("SimuNet","real") & samp.eff == 50 & group.rep == 1],
                 binwidth = 1,alpha = 0.15,position = "identity")+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(5)[2:5]))+
  scale_fill_manual(values = c("tomato","royalblue",gg_color_hue(5)[2:5]),labels = type.str)+
  labs(fill = "Network type",)+
  guides(colour = "none")+
  xlab("Edge weight")+ylab("Count")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results$edge.dt |>
  subset(samp.eff == 50 & group.rep == 1) |>
  ggplot(aes(weight, colour = type, fill = type))+
  facet_grid(i ~ j,scales = "free_y")+
  geom_density(bw = 1,alpha = 0.25,position = "identity",colour = NA)+
  geom_density(data = results$edge.dt[type %in% c("SimuNet","real") & samp.eff == 50 & group.rep == 1],
               bw = 1,alpha = 0.15,position = "identity")+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(5)[2:5]))+
  scale_fill_manual(values = c("tomato","royalblue",gg_color_hue(5)[2:5]),labels = type.str)+
  labs(fill = "Network type",)+
  guides(colour = "none")+
  xlab("Edge weight")+ylab("Frequency")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results$edge.dt |>
  subset(samp.eff == 50 & group.rep == 1) |>
  ggplot(aes(weight, colour = type))+
  facet_grid(i ~ j,scales = "free_y")+
  stat_ecdf()+
  stat_ecdf(data = results$edge.dt[type %in% c("SimuNet","real") & samp.eff == 50 & group.rep == 1],size = 1.2)+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(5)[2:5]))+
  labs(colour = "Network type",)+
  xlab("Edge weight")+ylab("Cumulative probability")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results$edge.comp.dt[order(mean.diff)] |>
  dplyr::mutate(type = paste0(y)) |>
  select(c("type","i","j","mean.diff","KS.stat","KS.p","KL","JS","EMD.e","EMD.m")) |>
  GGally::ggpairs(mapping = aes(colour = type,alpha = 0.2),columns = c("mean.diff","KS.stat","KS.p","KL","JS","EMD.e","EMD.m"))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal()

# Managing batches ----
library(ggridges)

results.multi$edge.comp.dt[,.N,by = .(i,j,samp.eff,y)]

results.multi$edge.comp.dt |>
  # subset(i == 4 & j == 7 ) |>
  ggplot(aes(x = KL,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.03)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.8,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"50"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

results.multi$edge.comp.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = JS,y = y,colour = as.factor(samp.eff),fill = as.factor(samp.eff)))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,samp.eff,y)),alpha = 0.1,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,samp.eff,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(8),"50"))+
  scale_fill_manual(values = paste0(gg_color_hue(8),"50"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,1))+
  ylab("")+xlab("Jensen-Shannon distance")+
  guides(fill = "none")+
  theme_minimal(12)

## Plotting results ----
x.labs <- paste0("... ",c("the real network","another network (similar generation)","an Erdős–Rényi graph (p = 0.5)","a Random network (fixed prob.)","a Random network (variable prob.)"))

KS.stat <-
  results.multi$edge.comp.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = KS.stat,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kolmogorov-Smirnoff statistic")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

KS.p <-
  results.multi$edge.comp.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = KS.p,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kolmogorov-Smirnoff p-value")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

KL <-
  results.multi$edge.comp.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = KL,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

JS <-
  results.multi$edge.comp.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = JS,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Jensen-Shannon distance")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

EMD.e <-
  results.multi$edge.comp.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = EMD.e,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Earth Mover Distance (Euclidean dist.)")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

EMD.m <-
  results.multi$edge.comp.dt |>
  subset(i <= 3 & j <= 3) |>
  ggplot(aes(x = EMD.m,y = y,colour = y,fill = y))+
  facet_grid(i ~ j)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),jittered_points = TRUE,alpha = 0.01,scale = .9)+
  geom_density_ridges(aes(group = interaction(i,j,group.rep,y)),fill = NA,alpha = 0.3,scale = .9,
                      jittered_points = TRUE,point_shape = "|",point_size = 2.5,
                      position = position_points_jitter(height = 0))+
  # geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_colour_manual(values = paste0(gg_color_hue(5),"13"))+
  scale_y_discrete(labels = x.labs)+
  scale_x_continuous(limits = c(0,NA))+
  ylab("")+xlab("Earth Mover Distance (Manhattan dist.)")+
  guides(colour = "none",fill = "none")+
  theme_minimal(12)

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))



results.multi$edge.comp.dt %>%
  .[,by = .(y,n,samp.eff,n.rep,group.rep),
     .(
       mean.diff = median(mean.diff),
       mean.diff.min = quantile(mean.diff,probs = 0.025),
       mean.diff.max = quantile(mean.diff,probs = 0.975),
       KS.stat   = median(KS.stat),
       KS.stat.min = quantile(KS.stat,probs = 0.025),
       KS.stat.max = quantile(KS.stat,probs = 0.975),
       KS.p      = median(KS.p),
       KS.p.min = quantile(KS.p,probs = 0.025),
       KS.p.max = quantile(KS.p,probs = 0.975),
       KL        = median(KL),
       KL.min = quantile(KL,probs = 0.025),
       KL.max = quantile(KL,probs = 0.975),
       JS        = median(JS),
       JS.min = quantile(JS,probs = 0.025),
       JS.max = quantile(JS,probs = 0.975),
       EMD.e     = median(EMD.e),
       EMD.e.min = quantile(EMD.e,probs = 0.025),
       EMD.e.max = quantile(EMD.e,probs = 0.975),
       EMD.m     = median(EMD.m),
       EMD.m.min = quantile(EMD.m,probs = 0.025),
       EMD.m.max = quantile(EMD.m,probs = 0.975)
     )
  ] |>
  ggplot(aes(y,KL,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge(jitter.width = .5,dodge.width = 2,seed = 42))+
  geom_errorbar(aes(ymin = KL.min, ymax = KL.max),alpha = 0.2,position = position_jitterdodge(jitter.width = .5,dodge.width = 2,seed = 42),width = 0.1)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  # coord_flip()+
  theme_minimal(12)



