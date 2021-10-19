# Load packages and custom functions ----
source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")

# Case study ----
## Generating main variables ----
set.seed(42)
n <- 10L

samp.eff <- 50L
n.rep <- 140L
n.each <- 10L
# ### Generating the group ----
# ind <- data.table(
#   name = as.character(1:n),
#   value = runif(n),
#   rank = sample(1:n,n)
# )
#
# ### Generating another group ----
# ind.other <- data.table(
#   name = as.character(1:n),
#   value = runif(n),
#   rank = sample(1:n,n)
# )

## Simulations ----
### Exporting objects to parallel workers ----
cl <- snow::makeCluster(7)
snow::clusterCall(
  cl,function() {
    source(".WIP/experiments_on_modelling_networks_generated_by_complex_processes_tools.R")
  }
)

### Running the simulations ----
param.samp.eff <- c(10,25,50,75,100,150,250,500)
param.list <-
  data.table(
    n = 10L,
    samp.eff = rep(param.samp.eff,each = n.each),
    n.rep = n.rep,
    ind = rep(each = n.each,
              lapply(
                1:length(param.samp.eff),
                \(r) data.table(name = as.character(1:n),value = runif(n),rank = sample(1:n,n))
              )
    ),
    ind.rep = rep(1:n.each,times = length(param.samp.eff))
  )

results.list <-
  lapply(
    1:nrow(param.list),
    function(r) {
      n <- param.list$n[r]; samp.eff <- param.list$samp.eff[r]; n.rep <- param.list$n.rep[r]
      ind <- param.list$ind[[r]]; ind.rep <- param.list$ind.rep[r]
      cat(paste0("\nn = ",n," - samp.eff = ",samp.eff," - ind.rep = ",ind.rep,
                 " (params ",r,"/",nrow(param.list),")\n"))
      cat("________________________________\n\n")
      results <-
        run_experiment(
          netgen_fun = highRankChooseMostValuable_generate_asso,
          netgen_args = list(samp.eff = samp.eff,ind = ind),
          netgen_other = list(samp.eff = samp.eff,ind = ind.other),
          n = n,samp.eff = samp.eff,names = ind$name,n.rep = n.rep,cl = cl
        )
      results$edge.dt <- cbind(results$edge.dt,ind.rep = ind.rep)
      results$edge.comp.dt <- cbind(results$edge.comp.dt,ind.rep = ind.rep)
      results
    }
  )
snow::stopCluster(cl)

results <- do.call(rbind_results,results.list)

# results <-
#   run_experiment(
#     netgen_fun = highRankChooseMostValuable_generate_asso,
#     netgen_args = list(samp.eff = samp.eff,ind = ind),
#     netgen_other = list(samp.eff = samp.eff,ind = ind.other),
#     n = n,samp.eff = samp.eff,names = ind$name,n.rep = 21,cl = cl
#   )

## Plotting results ----
x.labs <- paste0("... ",c("the real network","another network (similar generation)","an Erdős–Rényi graph (p = 0.5)","a Random network (fixed prob.)","a Random network (variable prob.)"))

KS.stat <-
  results$edge.comp.dt |>
  ggplot(aes(y,KS.stat,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Kolmogorov-Smirnoff statistic")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)

KS.p <-
  results$edge.comp.dt |>
  ggplot(aes(y,KS.p,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Kolmogorov-Smirnoff p-value")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)

KL <-
  results$edge.comp.dt |>
  ggplot(aes(y,KL,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Kullback-Leibler divergence")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)

JS <-
  results$edge.comp.dt |>
  ggplot(aes(y,JS,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Jensen-Shannon distance")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)

EMD.e <-
  results$edge.comp.dt |>
  ggplot(aes(y,EMD.e,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Earth Mover Distance (Euclidean dist.)")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)

EMD.m <-
  results$edge.comp.dt |>
  ggplot(aes(y,EMD.m,colour = y,fill = y))+
  geom_point(alpha = 0.2,position = position_jitterdodge())+
  geom_boxplot(alpha = 0.4)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  scale_x_discrete(labels = x.labs)+
  xlab("")+ylab("Earth Mover Distance (Manhattan dist.)")+
  guides(colour = "none",fill = "none")+
  coord_flip()+
  theme_minimal(12)

gridExtra::grid.arrange(KS.stat,KS.p,KL,JS,EMD.e,EMD.m,ncol = 2,
                        top = grid::textGrob("SimuNet compared to...", gp=grid::gpar(fontsize = 16,fontface = "bold")))

# Graphic exploration ---------------------------------------------------------------------------------------
## Data handling and wrangling ----
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
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
      weight = results$real.dist[,,1][.]
    )
  } |>
  data.table()

real.dt <-
  results$edge.dt |>
  subset(type == "real") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

SimuNet.dt <-
  results$edge.dt |>
  subset(type == "SimuNet") %>%
  .[,by = .(i,j,n,samp.eff,n.rep,type),
    .(weight = mean(weight))
  ]

## plots ----
results$edge.dt |>
  ggplot(aes(weight, colour = type, fill = type))+
  facet_wrap(. ~ paste0("edge ",i,"-",j),nrow = n,scales = "free_y")+
  geom_histogram(binwidth = 1,alpha = 0.15,position = "identity",colour = NA)+
  geom_histogram(data = results$edge.dt[type %in% c("SimuNet","real")],
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
  ggplot(aes(weight, colour = type, fill = type))+
  facet_grid(i ~ j,scales = "free_y")+
  geom_density(bw = 1,alpha = 0.15,position = "identity",colour = NA)+
  geom_density(data = results$edge.dt[type %in% c("SimuNet","real")],
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
  ggplot(aes(weight, colour = type))+
  facet_grid(i ~ j,scales = "free_y")+
  stat_ecdf()+
  stat_ecdf(data = results$edge.dt[type %in% c("SimuNet","real")],size = 1.2)+
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
test <- rbind(
  cbind(results[[7]],test.type = "original"),
  cbind(results.2[[7]],test.type = "second"),
  cbind(results.3[[7]],test.type = "provided"),
  cbind(results.4[[7]],test.type = "provided bis")
)


test$test.type <- test$test.type |> factor(levels = c("original","second","provided","provided bis"))
test |>
subset(type %in% c("real","SimuNet")) |>
ggplot(aes(weight, colour = test.type, fill = test.type))+
facet_grid(i ~ j,scales = "free_y")+
geom_density(aes(lty = type),bw = 1,alpha = 0.3,position = "identity")+
labs(fill = "Network type",)+
guides(colour = "none")+
xlab("Edge weight")+ylab("Frequency")+
geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
theme_minimal(16)+theme(legend.position = "top")


