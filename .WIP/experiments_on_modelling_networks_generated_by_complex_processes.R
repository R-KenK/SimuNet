# Load packages and custom functions ----
source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")

# Case study ----
## Generating main variables ----
set.seed(42)

## Simulations ----
n.rep <- 105L
n.group <- 3L
n.each <- 49L

### Running the simulations ----
param.n <-
  #seq(6,30,by = 3)
  c(5,8,10,15,20)[c(1,3)]
param.samp.eff <-
  # seq(20,500,by = 20)
  c(10,25,50,75,100,150,250,500)[c(4,5)]

param.netgen <-
  list(
    # data.table(netgen_name = "HRCMV",netgen_fun = list(highRankChooseMostValuable_generate_asso)),
    data.table(netgen_name = "GWC",netgen_fun = list(gregWithinClique_generate_asso))
  )

param.list <-
  generate_paramList(param.n = param.n,
                     param.samp.eff = param.samp.eff,
                     param.netgen = param.netgen,
                     n.rep = n.rep,
                     n.group = n.group,
                     n.each = n.each
  )[group.number <= 2 & group.rep <= 28]
param.list[]

### Exporting objects to parallel workers ----
cl <- prepare_clusters(7,param.list)
results <- copy(param.list)
start.time <- Sys.time()
start.time
results <- run_simulations(results,n.steps = 3)
results <- aggregate_edgeDT(results,cl = cl,n.steps = 2)
results <- prepare_for_distances(results)

snow::clusterCall(cl,\() {rm(results, pos = globalenv()); gc()})
snow::clusterExport(cl,list("results"),envir = environment())

results <- measure_all_distances(results,cl = cl)

end.time <- Sys.time()
end.time
end.time - start.time

results


saveRDS(results,".WIP/simulation.data/results.GWC.n_8.10.15.seff_75.100.150.rds")

results <- readRDS(".WIP/simulation.data/")
results

# test bench measure distances ----


test.results <- copy(results)
results <- prepare_for_distances(results,"real","SimuNet")
snow::clusterExport(cl,list("results"))
results <- measure_distances_vec(results,cl = cl)
results$distances[[220]]
testSN <- measure_distances2(results,"real","SimuNet",cl = NULL)
testrealB <- measure_distances2(results,"real","real.bis",cl = NULL)
testER <- measure_distances2(results,"real","ER",cl = NULL)

list(testSN,testrealB,testER) |>
  split_seq() |>
  lapply(do.call,what = rbind)

{
  do.call(mapply,list(FUN = rbind,... = .,SIMPLIFY = FALSE))
}

  mapply(FUN = rbind,SIMPLIFY = FALSE)

test.results[,edge.weight := mapply(cbind,
                                  n = n,samp.eff = samp.eff,n.rep = n.rep,group.number = group.number,group.rep = group.rep,edge.weight,SIMPLIFY = FALSE)]
test.results[]
test.results$edge.weight[[220]]

## Plotting distances ----
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

x.fill <- c("#FF6347","#4169E1",gg_color_hue(5)[2:5]) # "tomato" and "royalblue"
results[,edge.distance.dt :=lapply(edge.distance.dt,\(ddt) ddt$type <- ddt$type |>
                   factor(levels = c("real","real.bis","SimuNet","other","ER","fixed.rand","total.rand"))
          )]
### Single case as boxplots ----
KS.stat <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  ggplot(aes_string(x = "type",y = "KS.stat",colour = "type",fill = "type"))+
  geom_point(aes_string(group = "interaction(i,j,type)"),alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(aes_string(group = "interaction(i,j,type)"),alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_vline(xintercept = 1.5,size = 1.2)+
  # scale_x_discrete(labels = x.labs)+
  scale_fill_manual(values = x.fill)+
  scale_colour_manual(values = x.fill)+
  # xlab("")+ylab(label = ylab)+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)+theme(axis.text.y = element_text(face = x.labs.face))
  plot_distance(dist = "KS.stat",ylab = "Kolmogorov-Smirnoff statistic")

KS.p <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "KS.p",ylab = "Kolmogorov-Smirnoff p-value")

KL <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "KL",ylab = "Kullback-Leibler divergence")

JS <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "JS",ylab = "Jensen-Shannon distance")

EMD.e <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "EMD.e",ylab = "Earth Mover Distance (Euclidean dist.)")

EMD.m <-
  results |>
  subset(samp.eff == 75) |>
  subset(group.number == 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance(dist = "EMD.m",ylab = "Earth Mover Distance (Manhattan dist.)")

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))


### Multiple cases as densities ----
KS.stat <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "KS.stat",xlab = "Kolmogorov-Smirnoff statistic",x.lims = c(0,1))
  # facet_grid(. ~ samp.eff)

KS.p <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "KS.p",xlab = "Kolmogorov-Smirnoff p-value",x.lims = c(0,1))#+
# facet_grid(. ~ samp.eff)

KL <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "KL",xlab = "Kullback-Leibler divergence",x.lims = c(0,NA))#+
# facet_grid(. ~ samp.eff)

JS <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "JS",xlab = "Jensen-Shannon distance",x.lims = c(0,NA))#+
  # facet_grid(. ~ samp.eff)

EMD.e <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "EMD.e",xlab = "Earth Mover Distance (Euclidean dist.)",x.lims = c(0,NA))#+
# facet_grid(. ~ samp.eff)

EMD.m <-
  results |>
  subset(n == 10) |>
  subset(samp.eff == 75) |>
  subset(group.number <= 1) |>
  extract_column("edge.distance.dt") |>
  plot_distance_distrib(dist = "EMD.m",xlab = "Earth Mover Distance (Manhattan dist.)",x.lims = c(0,NA))#+
# facet_grid(. ~ samp.eff)

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...",
                                             gp=grid::gpar(fontsize = 16,fontface = "bold")))

# Graphic exploration ---------------------------------------------------------------------------------------
## single sim ----
type.str <- c("SimuNet network",
              "the real network",
              "the real network bis",
              "another network (similar generation)",
              "an Erdős–Rényi graph (p = 0.5)",
              "a Random network (fixed prob.)",
              "a Random network (variable prob.)")
n.single <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(n)

Adj.obs.dt <-
  expand.grid(i = 1:n.single,j = 1:n.single) |>
  subset(i < j) |>
  as.matrix() %>%
  {
    cbind(
      .,
      weight = {Adj <- results |>
        subset(samp.eff == 50) |>
        subset(group.number == 1 & group.rep == 1) %>%
        {pull(.,dists)[[1]]} %>%
        {.$real.dist[[1]][,,1]};Adj[.]}
    )
  } |>
  data.table()

real.dt <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(edge.dt) |>
  rbindlist() |>
  subset(type == "real") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

SimuNet.dt <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(edge.dt) |>
  rbindlist() |>
  subset(type == "SimuNet") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

## plots ----
results.single <-
  results |>
  subset(samp.eff == 50) |>
  subset(group.number == 1 & group.rep == 1) |>
  pull(edge.dt) |>
  rbindlist()
results.single |>
  ggplot(aes(weight, colour = type, fill = type))+
  facet_wrap(. ~ paste0("edge ",i,"-",j),nrow = n.single,scales = "free_y")+
  geom_histogram(binwidth = 1,alpha = 0.15,position = "identity",colour = NA)+
  geom_histogram(data = results.single[type %in% c("SimuNet","real")],
                 binwidth = 1,alpha = 0.15,position = "identity")+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]))+
  scale_fill_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]),labels = type.str)+
  labs(fill = "Network type",)+
  guides(colour = "none")+
  xlab("Edge weight")+ylab("Count")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results.single |>
  ggplot(aes(weight, colour = type, fill = type))+
  facet_grid(i ~ j,scales = "free_y")+
  geom_density(bw = 1,alpha = 0.25,position = "identity",colour = NA)+
  geom_density(data = results.single[type %in% c("SimuNet","real")],
               bw = 1,alpha = 0.15,position = "identity")+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]))+
  scale_fill_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]),labels = type.str)+
  labs(fill = "Network type",)+
  guides(colour = "none")+
  xlab("Edge weight")+ylab("Frequency")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results.single |>
  ggplot(aes(weight, colour = type))+
  facet_grid(i ~ j,scales = "free_y")+
  stat_ecdf()+
  stat_ecdf(data = results.single[type %in% c("SimuNet","real")],size = 1.2)+
  geom_vline(data = SimuNet.dt,aes(xintercept = weight),
             colour = "tomato",lty = "dashed",size = 1.1)+
  geom_vline(data = real.dt,aes(xintercept = weight),
             colour = "royalblue",lty = "dashed",size = 1.1)+
  geom_point(data = Adj.obs.dt,aes(x = weight,y = 0),fill = "white",
             colour = "tomato",shape = 21,stroke = 1.5,size = 2)+
  scale_colour_manual(values = c("tomato","royalblue",gg_color_hue(6)[2:6]))+
  labs(colour = "Network type",)+
  xlab("Edge weight")+ylab("Cumulative probability")+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal(16)+theme(legend.position = "top")

results$edge.distance.dt[order(mean.diff)] |>
  dplyr::mutate(type = paste0(y)) |>
  select(c("type","i","j","mean.diff","KS.stat","KS.p","KL","JS","EMD.e","EMD.m")) |>
  GGally::ggpairs(mapping = aes(colour = type,alpha = 0.2),columns = c("mean.diff","KS.stat","KS.p","KL","JS","EMD.e","EMD.m"))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  theme_minimal()

# Managing batches ----
results |>
  pull(edge.distance.dt) |>
  rbindlist() %>%
  .[,.N,by = .(i,j,samp.eff,y)]

results.multi$edge.distance.dt |>
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

results.multi$edge.distance.dt |>
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
  results.multi$edge.distance.dt |>
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
  results.multi$edge.distance.dt |>
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
  results.multi$edge.distance.dt |>
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
  results.multi$edge.distance.dt |>
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
  results.multi$edge.distance.dt |>
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
  results.multi$edge.distance.dt |>
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



results.multi$edge.distance.dt %>%
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



