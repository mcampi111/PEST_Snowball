library(igraph)
library(sna)
library(DirectedClustering)
library(snowboot)
library(network)
library(intergraph)
library(ggplot2)
library(tidyr)
library(reshape2)
library(dplyr)
library(purrr)
library(data.table)

####################
#  SET DIRECTORIES #
####################

my_dir_figs = "C:\\Users\\Marta\\Desktop\\New_Project_Dimitris\\PEST\\paper\\figs\\"


#############
#  SET SEED #
#############

set.seed(1)


###############
#  FUNCTIONS  #
###############

#This function computes the snowball functions with the same package snowboot but 
#each seed, at each wave, recruits a percentage of its own edges, i.e. delta parameter 


#\delta is a function of each node that recruits, i.e. it is a function of 

#NOTE THAT  I USED CEILING RATHER THAN ROUND SO IT ROUNDS UP


my_snowball_fun<- function(G, n_seed, lsmi_union, n_waves, delta){
  
  
  #G_rp = round(n_seed*delta) # this recruits this multiplication number of nodes
  #G_rp = delta # this recruits exactly this number of nodes
  
  first_wave<- lapply(1:n_seed,function(l) lsmi_union(G, n.seed = l, n.wave = 1) )
  first_wave_lsmi <- lapply(1:n_seed,function(o) first_wave[[o]]$lsmi_big)
  
  G_rp = lapply(1:n_seed, function(t)  ceiling(G$degree[match(mapply(`[[`, first_wave_lsmi[[t]], 1), 1:G$n)]*delta) ) 
  
  
  
  first_wave_refs<- c(
    list(sample(lapply(1:n_seed, 
                       function(t) mapply(`[[`, first_wave_lsmi[[t]], 2) )[[1]], G_rp[[1]] , FALSE)),
    
    lapply(2:n_seed, function(x)
      lapply(1:length(G_rp[[x]]), function(s)
        sample( lapply(1:n_seed,function(y) mapply(`[[`, first_wave_lsmi[[y]], 2) )[[x]][[s]], 
                G_rp[[x]][[s]],  FALSE) )) )
  
  
  
  # first_wave_refs<- c(list(sample(lapply(1:n_seed, function(t) mapply(`[[`, first_wave_lsmi[[t]], 2) )[[1]], G_rp[[1]] , FALSE)),
  #                     lapply(2:n_seed, function(x)
  #                       lapply(lapply(1:n_seed, 
  #                                     function(y) mapply(`[[`, first_wave_lsmi[[y]], 2) )[[x]], sample, G_rp[[x]][[y]], FALSE) ) )
  # 
  
  seeds = map(lapply(1:n_seed, function(seed) first_wave[[seed]]$sequence_seeds),unlist)
  
  wave_refs_fin<- map(first_wave_refs, unlist)
  
  
  refs_other_waves<- list()
  refs_other_waves[[1]]<- seeds
  refs_other_waves[[2]]<- wave_refs_fin
  
  for (i in 1:(n_waves-1)) {
    
    seed_wave<- wave_refs_fin
    
    wave_s = lapply(1:(n_seed), function(j) lapply(1:length(seed_wave[[j]]), function(h) 
      lsmi_union(G, 
                 n.wave = 1, 
                 seeds= seed_wave[[j]][[h]],
                 n.seeds = length(seed_wave[[j]][[h]]) ) ))
    
    
    wave_s_lsmi<- lapply(1:(n_seed), function(j)
      lapply(1:length(seed_wave[[j]]), function(h) wave_s[[j]][[h]]$lsmi_big ))
    
    
    G_rp = lapply(1:(n_seed), function(j)
      lapply(1:length(seed_wave[[j]]), function(h) 
        ceiling(G$degree[match(lapply(wave_s_lsmi[[j]][[h]], `[[`, 1), 1:G$n)]*delta)  ))
    
    wave_s_refs<- lapply(1:(n_seed), function(j)
      lapply(1:length(seed_wave[[j]]), function(h)  
        lapply( lapply(wave_s_lsmi[[j]][[h]], `[[`, 2) , sample, G_rp[[j]][[h]], FALSE ) ) )  #you could this with mapply but change the order of the arguments
    
    wave_refs_fin<- map(wave_s_refs, unlist)
    
    refs_other_waves[[i+2]]<- wave_refs_fin
    
  }
  
  
  
  all_waves_refs<- do.call(Map, c(f = list, refs_other_waves))
  
  return(all_waves_refs)
  
}

#function running my_snowball_fun returning list divided by seeds and in each seeds
#the identified nodes per wave - NB: the firs entry, i.e. prova[[1]][[1]], prova[[2]][[1]], etc. provides the
#index of the seeds while the rest are the  nodes identified at each wave
prova<- my_snowball_fun(G_d_to_net[[4]][[3]],
                        n_seed = n_seed,
                        lsmi_union, 
                        n_waves = n_waves , 
                        delta = delta)






#################################################################################
# Generate (?) Weighted and Undirected graph with Erdos-Renyi distribution      #
#################################################################################

# Note about the simulated network:
# studied network is the recruited network which is different from the nomination network
# although we recognize that there is a difference between recruitment and nomination
# we simulate the network of reputation and only study that


#F. and Snjder algorithm (find that on the other packages)


#DIMITRIS: PLAY WIT THIS BELOW to identify the correct configuration

gne = barabasi.game(n_nodes[1], 
              power = 1, 
              m = 10, 
              out.dist = NULL, 
              out.seq = NULL, 
              out.pref = FALSE, 
              zero.appeal = 1, 
              directed = FALSE)

plot(gne)


#PARAMETERS OF THE SIMULATED NETWORK

n_nodes = c(50,100,150,200)

n_edges = list(c(10, 20, 40), 
               c(20, 40, 80),
               c(30, 60, 120),
               c(60, 120, 240) )

#MARTA PLAYS WITH THIS
n_edges = list(c(5, 10, 20), 
               c(5, 10, 20),
               c(5, 10, 20),
               c(5, 10, 20) )

n_nodes[1]/n_edges[[1]]
n_nodes[2]/n_edges[[2]]
n_nodes[3]/n_edges[[3]]
n_nodes[4]/n_edges[[4]]

n_ratio = n_nodes[1]/n_edges[[1]]


G_d<- lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j)
    barabasi.game(n_nodes[i], 
                  power = 1, 
                  m = n_edges[[i]][[j]], 
                  out.dist = NULL, 
                  out.seq = NULL, 
                  out.pref = FALSE, 
                  zero.appeal = 1, 
                  directed = FALSE) ))


G_d[[1]][[1]] #first index select nodes, second index select p of drawing an edge


#play with the below function to check how to recruit

min_degree = do.call(rbind,map( sapply(1:length(n_nodes), function(i)
  sapply(1:length(n_edges[[i]]),function(j) min(igraph::degree(G_d[[i]][[j]])) )), unlist ))

max_degree = do.call(rbind,map( sapply(1:length(n_nodes), function(i)
  sapply(1:length(n_edges[[i]]),function(j) max(igraph::degree(G_d[[i]][[j]])) )), unlist ))

mean_degree = do.call(rbind,map( sapply(1:length(n_nodes), function(i)
  sapply(1:length(n_edges[[i]]),function(j) mean(igraph::degree(G_d[[i]][[j]])) )), unlist ))

df_degree = data.frame(min_degree,max_degree,mean_degree)

df_degree_r_names = do.call(rbind, map(sapply(1:length(n_nodes), function(i)
  sapply(1:length(n_edges[[i]]),function(j) 
    paste("n. nodes ", n_nodes[i], " n. edges ", n_edges[[i]][[j]], sep = "") )), unlist))

rownames(df_degree)<- df_degree_r_names

df_degree

df_degree$approx_mean_degree = rep( round(df_degree$mean_degree[1:3],0), length(n_nodes))

mean_degree_label<- round(df_degree$mean_degree[1:3],0)


df_degree$n_nodes = rep(n_nodes, each =3)

df_degree_melt = reshape2::melt(df_degree[,-3], id = c("n_nodes", "approx_mean_degree"))

gg = ggplot(df_degree_melt, aes(x = factor(approx_mean_degree) , y = factor(n_nodes), fill = value))
gg = gg +   facet_grid(~variable )
gg = gg +  geom_tile(aes(fill = value))
gg = gg +  geom_text(aes(label = round(value, 3)), size = 5)
gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
gg = gg + theme_bw()
#gg = gg + scale_x_discrete(labels= approx_mean_degree)
gg = gg + scale_y_discrete(labels= n_nodes)
gg = gg + xlab("Approx. Mean Degree") + ylab("N. Nodes")
gg = gg + theme(legend.position = "none")
gg = gg + scale_fill_gradient(low = "white", high = "steelblue")
gg


ggsave(filename =  paste(my_dir_figs,"s_f_min_max_degree.pdf", sep = ""),
       plot = gg)



#########
#WEIGHTS# (we decide to drop this)
#########

# w<- lapply(1:length(n_nodes), function(i)
#   sapply(1:length(n_edges[[i]]),function(j)   runif(length(E(G_d[[i]][[j]])), 0, 1) ))
# 
# 
# #plug the weights into the configuration of the graphs
# for (i in 1:length(n_nodes)){
#   for (j in 1:length(n_edges[[i]])){
#     
#     E(G_d[[i]][[j]])$weight <- w[[i]][[j]]
#     
#   }
# }


#######################
#COMPUTE THE ADJACENCY#
#######################

# A<- lapply(1:length(n_nodes), function(i)
#   lapply(1:length(n_edges[[i]]),function(j)  get.adjacency(G_d[[i]][[j]], sparse=FALSE, attr="weight") ))

##########################################
#COMPUTE Barrat et al. (2004) coefficient# --> you can use this only if you weight the graphs
##########################################
###################################################################
#reference: https://www.pnas.org/content/pnas/101/11/3747.full.pdf#
###################################################################

# BarratClust<- lapply(1:length(n_nodes), function(i)
#   lapply(1:length(n_edges[[i]]),function(j)    ClustBCG(A[[i]][[j]], "undirected")  ))
# 
# 
# clust_coef_global = lapply(1:length(n_nodes), function(i)
#   sapply(1:length(n_edges[[i]]),function(j)  BarratClust[[i]][[j]]$GlobalCC ))
# 
# clust_coef_global_df = do.call(cbind, clust_coef_global)
# colnames(clust_coef_global_df)<- n_nodes
# rownames(clust_coef_global_df)<- n_ratio
# 
# df_cc = reshape2 ::melt(clust_coef_global_df)
# colnames(df_cc)<- c("n_ratio", "n_nodes", "value")
# 
# 
# gg = ggplot(df_cc, aes(x = factor(n_ratio) , y = factor(n_nodes), fill = value)) 
# gg = gg +  geom_tile(aes(fill = value))  
# gg = gg +  geom_text(aes(label = round(value, 3)), size = 5) 
# gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
# gg = gg + theme_bw()
# gg = gg + scale_x_discrete(labels= n_ratio)
# gg = gg + scale_y_discrete(labels= n_nodes)
# gg = gg + xlab("ratio nodes/edges") + ylab("n. nodes")
# gg = gg + labs(fill = "Cluster Coefficient")
# gg = gg + scale_fill_gradient(low = "white", high = "steelblue") 
# gg

##########################
# DENSITY of the network #
##########################

g_dens<- lapply(1:length(n_nodes), function(i)
  sapply(1:length(n_edges[[i]]),function(j)  edge_density(G_d[[i]][[j]], loops=TRUE) ))


g_dens_df = do.call(cbind, g_dens)
colnames(g_dens_df)<- n_nodes
rownames(g_dens_df)<- mean_degree_label

df_dens = reshape2::melt(g_dens_df)
colnames(df_dens)<- c("mean_degree", "n_nodes", "value")


gg = ggplot(df_dens, aes(x = factor(mean_degree) , y = factor(n_nodes), fill = value)) 
gg = gg +  geom_tile(aes(fill = value))  
gg = gg +  geom_text(aes(label = round(value, 3)), size = 5) 
gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
gg = gg + theme_bw()
gg = gg + scale_x_discrete(labels= mean_degree)
gg = gg + scale_y_discrete(labels= n_nodes)
gg = gg + xlab("Approx. Mean Degree") + ylab("N. Nodes")
gg = gg + labs(fill = "Density")
gg = gg + scale_fill_gradient(low = "white", high = "steelblue") 
gg



ggsave(filename =  paste(my_dir_figs,"s_f_dens_nets.pdf", sep = ""),
       plot = gg)


#####################
# SNOWBALL SAMPLING #
#####################

#Convert the class of functions to the network package which incorporates the snowball sampling
G_d_to_net <-  lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j) igraph_to_network(G_d[[i]][[j]]) ))

G_d_to_net[[1]][[3]]


#PARAMETRS snow
n_seed = 5
n_waves = 3

delta = 0.3


#snowball sampling - function snowboot 
#here, all the neighbors are selected - ASSUMPTION: 100% of the neighbors of a node are recruited
snow_samp = lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j)
    lapply(1:n_seed, function(h) 
      lapply(1:n_waves, function(k) 
        lsmi_union(G_d_to_net[[i]][[j]], 
                   n.seed = h, 
                   n.wave = k, 
                   seeds = NULL) ))))



# ASSUMPTION: delta of the neighbors of a node are recruited
snow_samp_recr = lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j)
    my_snowball_fun(G_d_to_net[[i]][[j]],
                    n_seed = n_seed,
                    lsmi_union, 
                    n_waves = n_waves , 
                    delta = delta) )) #use delta[j] if a vector is employed otherwise
#if a scalar is used --> delta



n_seed_per_wave_unique = lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j)
    lapply(1:n_seed, function(h) 
      lapply(1:(n_waves), function(k)  
        length(unique(do.call(c, 
                              as.list(sapply(snow_samp[[i]][[j]][[h]][[k]]$lsmi_big, "[[",k+1))) )) ))))


n_seed_per_wave_recr =  lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j)
    lapply(1:n_seed, function(h)
      lapply(1:(n_waves+1), function(k)
        length(snow_samp_recr[[i]][[j]][[h]][[k]])   ))))


#FOR MARTA TO CHECK
# snow_samp[[1]][[1]][[1]][[1]] # first - n_nodes, second: n_edges, third seed and fourth for the wave
# 
# #1 seed
# length(do.call(c, as.list(sapply(snow_samp[[1]][[1]][[1]][[1]]$lsmi_big, "[[",2))) )#wave 1
# 
# length(do.call(c, as.list(sapply(snow_samp[[1]][[1]][[1]][[2]]$lsmi_big, "[[",2))) ) #wave 2
# length(do.call(c, as.list(sapply(snow_samp[[1]][[1]][[1]][[2]]$lsmi_big, "[[",3))) )
# 
# 
# length(do.call(c, as.list(sapply(snow_samp[[1]][[1]][[1]][[3]]$lsmi_big, "[[",2))) )# wave 3
# length(do.call(c, as.list(sapply(snow_samp[[1]][[1]][[1]][[3]]$lsmi_big, "[[",3))) )
# length(do.call(c, as.list(sapply(snow_samp[[1]][[1]][[1]][[3]]$lsmi_big, "[[",4))) )
# 
# 
# #2 seed
# length(do.call(c, sapply(snow_samp[[1]][[1]][[2]][[1]]$lsmi_big, "[[",2)) )#wave 1
# length(do.call(c, as.list(sapply(snow_samp[[1]][[1]][[2]][[1]]$lsmi_big, "[[",2))) )#wave 1
# 
# 
# length(do.call(c, sapply(snow_samp[[1]][[1]][[2]][[2]]$lsmi_big, "[[",2)) )# wave 2
# length(do.call(c, sapply(snow_samp[[1]][[1]][[2]][[2]]$lsmi_big, "[[",3)) )
# 
# length(do.call(c, sapply(snow_samp[[1]][[1]][[2]][[3]]$lsmi_big, "[[",2)) )# wave 3
# length(do.call(c, sapply(snow_samp[[1]][[1]][[2]][[3]]$lsmi_big, "[[",3)) )
# length(do.call(c, sapply(snow_samp[[1]][[1]][[2]][[3]]$lsmi_big, "[[",4)) )
# 
# #3 seed
# length(do.call(c, sapply(snow_samp[[1]][[1]][[3]][[1]]$lsmi_big, "[[",2)) )#wave 1
# 
# length(do.call(c, sapply(snow_samp[[1]][[1]][[3]][[2]]$lsmi_big, "[[",2)) )# wave 2
# length(do.call(c, sapply(snow_samp[[1]][[1]][[3]][[2]]$lsmi_big, "[[",3)) )
# 
# length(do.call(c, sapply(snow_samp[[1]][[1]][[3]][[3]]$lsmi_big, "[[",2)) )# wave 3
# length(do.call(c, sapply(snow_samp[[1]][[1]][[3]][[3]]$lsmi_big, "[[",3)) )
# length(do.call(c, sapply(snow_samp[[1]][[1]][[3]][[3]]$lsmi_big, "[[",4)) )



n_seed_wave_unique = lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j)  do.call(rbind, n_seed_per_wave_unique[[i]][[j]]) ))

n_seed_wave_recr = lapply(1:length(n_nodes), function(i)
  lapply(1:length(n_edges[[i]]),function(j)   do.call(rbind, n_seed_per_wave_recr[[i]][[j]])[,-1] ))




n_seed_wave_new_unique =  lapply(1:length(n_nodes), function(i) do.call(rbind, n_seed_wave_unique[[i]]))

n_seed_wave_new_recr =  lapply(1:length(n_nodes), function(i) do.call(rbind, n_seed_wave_recr[[i]]))



n_seed_df_unique =  lapply(1:length(n_nodes), function(i) data.frame(n_seed_wave_new_unique[[i]]))

n_seed_df_recr =  lapply(1:length(n_nodes), function(i) data.frame(n_seed_wave_new_recr[[i]]))


r_name_df_seed = lapply(1:length(n_edges), function(h) lapply(1:length(n_edges[[h]]), function(j) 
  sapply(1:5, function(i) paste("seed ", i, ", n. edges ", n_edges[[h]][[j]], sep = "" ))))



row_names_df_seed = lapply(1:length(r_name_df_seed), function(j) do.call(rbind, 
                                                                         lapply(1:length(r_name_df_seed[[j]]), function(i) cbind(r_name_df_seed[[j]][[i]]) )))



for (i in 1:length(n_nodes)) {
  
  colnames(n_seed_df_unique[[i]]) <- 1:(n_waves)
  rownames(n_seed_df_unique[[i]]) <-  row_names_df_seed[[i]]
  colnames(n_seed_df_recr[[i]]) <- 1:(n_waves)
  rownames(n_seed_df_recr[[i]]) <-  row_names_df_seed[[i]]
  
}




for (i in 1:length(n_nodes)){ n_seed_df_unique[[i]]$n_nodes = n_nodes[i] }
for (i in 1:length(n_nodes)){ n_seed_df_unique[[i]]$n_seed = rep(1:n_seed, length(n_edges[[1]]) ) }
for (i in 1:length(n_nodes)){ n_seed_df_unique[[i]]$n_edges = rep(n_edges[[i]], each = n_seed) }
for (i in 1:length(n_nodes)){ n_seed_df_unique[[i]]$mean_degree =  rep(mean_degree_label, each =n_seed) }

for (i in 1:length(n_nodes)){ n_seed_df_recr[[i]]$n_nodes = n_nodes[i] }
for (i in 1:length(n_nodes)){ n_seed_df_recr[[i]]$n_seed = rep(1:n_seed, length(n_edges[[1]]) ) }
for (i in 1:length(n_nodes)){ n_seed_df_recr[[i]]$n_edges = rep(n_edges[[i]], each = n_seed) }
for (i in 1:length(n_nodes)){ n_seed_df_recr[[i]]$mean_degree =  rep(mean_degree_label, each =n_seed)}

n_seed_df_new_unique = lapply(1:length(n_nodes), function(i)
  mutate_all(n_seed_df_unique[[i]], function(x) as.numeric(as.character(x)))  )

n_seed_df_new_recr = lapply(1:length(n_nodes), function(i)
  mutate_all(n_seed_df_recr[[i]], function(x) as.numeric(as.character(x)))  )



df_snow_unique = reshape2::melt(n_seed_df_new_unique, id = c("n_nodes","mean_degree" ,"n_seed", "n_edges"))
df_snow_unique = df_snow_unique %>%  mutate(value2 = value/n_nodes) 


df_snow_recr = reshape2::melt(n_seed_df_new_recr, id = c("n_nodes","mean_degree" ,"n_seed", "n_edges"))
df_snow_recr = df_snow_recr %>%  mutate(value2 = value/n_nodes) 
df_snow_recr = df_snow_recr %>%  mutate(value3 = ifelse(value2 > 1, 1, value2)) 


############
# Heatmaps #
############

#########
#numbers# --> DO NOT LOOK AT THIS
#########


#100% recruited (unique)
# gg = ggplot(df_snow_unique, aes(x = factor(variable) , y = factor(n_seed), fill = value)) 
# gg = gg +  geom_tile(aes(fill = value))  
# gg = gg +   facet_grid(n_nodes ~ mean_degree )
# gg = gg +  geom_text(aes(label = round(value, 3)), size = 3) 
# gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
# gg = gg + theme_bw()
# gg = gg + scale_x_discrete(labels= c(1:(n_waves)))
# gg = gg + scale_y_discrete(labels= c(1:n_seed))
# gg = gg + xlab("n. waves") + ylab("n. seed")
# gg = gg + labs(fill = "N. Recruited")
# gg = gg + scale_fill_gradient(low = "white", high = "steelblue") 
# gg
# 
# 
# 
# #delta% recruited 
gg = ggplot(df_snow_recr, aes(x = factor(variable) , y = factor(n_seed), fill = value))
gg = gg +  geom_tile(aes(fill = value))
gg = gg +   facet_grid(n_nodes ~ mean_degree )
gg = gg +  geom_text(aes(label = round(value, 3)), size = 3)
gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
gg = gg + theme_bw()
gg = gg + scale_x_discrete(labels= c(1:(n_waves)))
gg = gg + scale_y_discrete(labels= c(1:n_seed))
gg = gg + xlab("n. waves") + ylab("n. seed")
gg = gg + labs(fill = "N. Recruited")
gg = gg + scale_fill_gradient(low = "white", high = "steelblue")
gg



############
#percentage#
############


#100% recruited (unique)
gg = ggplot(df_snow_unique, aes(x = factor(variable) , y = factor(n_seed), fill = value2)) 
gg = gg +  geom_tile(aes(fill = value2))  
gg = gg +   facet_grid(n_nodes ~ mean_degree )
gg = gg +  geom_text(aes(label = round(value2, 3)), size = 3) 
gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
gg = gg + theme_bw()
gg = gg + scale_x_discrete(labels= c(1:(n_waves)))
gg = gg + scale_y_discrete(labels= c(1:n_seed))
gg = gg + xlab("n. waves") + ylab("n. seed")
gg = gg + labs(fill = "N. Recruited")
gg = gg + scale_fill_gradient(low = "white", high = "steelblue") 
gg



#delta% recruited 
gg = ggplot(df_snow_recr, aes(x = factor(variable) , y = factor(n_seed), fill = value3)) 
gg = gg +  geom_tile(aes(fill = value3))  
gg = gg +   facet_grid(n_nodes ~ mean_degree )
gg = gg +  geom_text(aes(label = round(value3, 2)), size = 3) 
gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
gg = gg + theme_bw()
gg = gg + scale_x_discrete(labels= c(1:(n_waves)))
gg = gg + scale_y_discrete(labels= c(1:n_seed))
gg = gg + xlab("n. waves") + ylab("n. seed")
gg = gg + labs(fill = "Perc. Recruited")
gg = gg + scale_fill_gradient(low = "white", high = "steelblue") 
gg


ggsave(filename =  paste(my_dir_figs,"s_f_snow_0.3_recr.pdf", sep = ""),
       plot = gg)

#########################
#  PLOT OF THE NETWORK  #
#########################

#example for one graph
df = melt(snow_samp[[1]][[1]][[1]])
df = df %>%  mutate(col = ifelse(L1 == 1, "black", ifelse(L1 ==2, "cyan", ifelse(L1 == 3, "green", "blue")))   ) 
df = df %>%  mutate(shape = ifelse(duplicated(df$value) == TRUE, "square", "circle" ))

miss_nodes = as.vector( V(G_d[[1]][[1]]))[which(!as.vector( V(G_d[[1]][[1]])) %in% sort(unique(df$value)))]
new_rows = cbind(miss_nodes, rep("5",length(miss_nodes)),rep("white",length(miss_nodes)), rep("circle",length(miss_nodes)) )
colnames(new_rows)<- colnames(df)

df_col = rbind(df, new_rows)
df_col_2 = df_col[order(as.numeric(df_col$value)),]
#df_col_ord =  data.table(df_col_2, key = 'value')
#df_col_ord_fin = df_col_ord[unique(df_col_ord[,list(value)]), mult= 'last']
df_col_ord_fin = df_col_2 %>%
  group_by(value) %>%
  slice(c(1, n())) %>% mutate( first_col= first(col), last_shape = last(shape)) %>%
  slice(1)

plot(G_d[[1]][[1]], vertex.color = df_col_ord_fin$first_col , 
     vertex.frame.color="red", 
     vertex.shape=df_col_ord_fin$last_shape,
     vertex.lty=2, 
     edge.label.color = "darkgrey")


#functions for more graphs

df_plot_net<- function(snow_samp_type, G_d_correspondent){
  
  n_seed_h = length(snow_samp_type)
  
  df_final = list()
  
  for (i in 1:n_seed_h) {
    
    df = reshape2::melt(snow_samp_type[[i]])
    df = df %>%  mutate(col = ifelse(L1 == 1, "black", ifelse(L1 ==2, "cyan", ifelse(L1 == 3, "green", "blue")))   ) 
    df = df %>%  mutate(shape = ifelse(duplicated(df$value) == TRUE, "square", "circle" ))
    
    miss_nodes = as.vector( V(G_d_correspondent))[which(!as.vector( V(G_d_correspondent)) %in% sort(unique(df$value)))]
    new_rows = cbind(miss_nodes, rep("5",length(miss_nodes)),rep("white",length(miss_nodes)), rep("circle",length(miss_nodes)) )
    colnames(new_rows)<- colnames(df)
    
    df_col = rbind(df, new_rows)
    df_col_2 = df_col[order(as.numeric(df_col$value)),]
    df_col_ord_fin = df_col_2 %>%
      group_by(value) %>%
      slice(c(1, n())) %>% mutate( first_col= first(col), last_shape = last(shape)) %>%
      slice(1)
    
    df_final[[i]] = df_col_ord_fin
    
  }
  
  return(df_final)
  
}

#how to pick snow_samp_type and G_d
#REMARK that 
length(G_d) #number of nodes according to:
n_nodes

length(G_d[[1]]) #we are in thel ist of graph with 50 nodes
#there are 4 configuration which are based on 
p_edge[1]
#to study the different results

#the same indices apply to the selected snowball sampling scheme
#Note that snow_samp assumes 100% recruitment and snow_samp_recr assumes 50% recruited 

index_nodes = 1
index_edges = 1
snow_samp_type = snow_samp_recr

toplot_snow_seeds = df_plot_net(snow_samp_type[[index_nodes]][[index_edges]], G_d[[index_nodes]][[index_edges]])

attach(mtcars)
par(mfrow=c(3,2), mar = c(1.1, 1.1, 3.1, 1.1) )

for (i in 1:n_seed) {
  
  plot(G_d[[index_nodes]][[index_edges]], vertex.color = toplot_snow_seeds[[i]]$first_col , 
       vertex.frame.color="red", 
       vertex.shape=toplot_snow_seeds[[i]]$last_shape,
       vertex.lty=2, 
       edge.label.color = "darkgrey")
  
  
}
title('Snowball with 5 seeds', outer = TRUE, line = -3)


#TODO - produce the same heatmaps of the obtained per number of seeds and number of waves
#############################################
# COMPUTE QUANTITIES OF THE INDUCED NETWORK #
#############################################





###################################################################################################################
###################################################################################################################
#                                FURTHER RESEARCH & OTHER RESEARCH QUESTIONS                                      #
###################################################################################################################
###################################################################################################################



##################################################################
# COMPUTE IN HOW MANY WAVES WE REACHED CENTRAL NODES OF THE NETS # - status reputation status 
##################################################################

#RESEARCH QUESTION HERE:
#what is the prob of you recruiting from the list of other actors and the process of recruiting
#ends up the underlying network


prova = igraph::betweenness(G_d[[1]][[1]])
prova

igraph::estimate_betweenness(G_d[[1]][[1]])


gne = rev(order(prova)) #from the most central to the last one
gne[1]

snow_samp_recr[[1]][[1]][[1]]

sapply(1:n_seed, function(j)
  unlist(sapply(1:(n_waves+1), function(i) which(gne[1] == snow_samp_recr[[1]][[1]][[j]][[i]]))) )



#we can think about this more
####################
# SEEDS CENTRALITY # centrality of seeds & number of seeds - the second one will already be answered
####################

#RESEARCH QUESTION HERE:
#PLAY WITH THE SEEDs, I.E. WHAT HAPPENS TO THE SNOWBALL SAMPLING IF THE SEEDS IS MORE OR LESS CENTRAL 


#####################
#   NEW TECHNIQUE   #
#####################

#nominate randomly whoever without a network structure and compare the two approaches
#do this and compare the two ideas


###############################
###############################
# SOME THOUGHTS & DISCUSSION ##
###############################
###############################

#M: check notes -->mynotebook

#at what level have we saturated and if we can claim --> at what level do we reach 50% --> to make claim 
#about the population


# what is the probability of responding (reached/sampled) given the centrality/degree of the node - this is 
#related to the fact if you are central - to provide a profile  your degree given this type of clustering 
#implies this probability of being selected (side line to our story)


#try to compute the centrality of the nodes
#try to study the relationship parent-children
#try to re-fit an Erdos-Renyi model and then study the properties of the identified network

#TODO:

#1) plot the network per wave  color
#2) compare the underlying network and the estimated network (quantities like density of the net and the
# cluster coeff.)
#3) rank the nodes according to their centrality
#4) check after how many waves we capture the above nodes



#find the central nodes of the nodes (most connected)
#after how many waves do I get them



#superforecasters (check this in the literature)










