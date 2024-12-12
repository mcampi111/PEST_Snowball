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

my_dir_figs = "path/"

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


my_snowball_fun<- function(G, n_seed, lsmi_union, n_waves, delta){
  
  
  #G_rp = round(n_seed*delta) # this recruits this multiplication number of nodes
  #G_rp = delta # this recruits exactly this number of nodes
  
  first_wave<- lapply(1:n_seed,function(l) lsmi_union(G, n.seed = l, n.wave = 1) )
  first_wave_lsmi <- lapply(1:n_seed,function(o) first_wave[[o]]$lsmi_big)
  
  G_rp = lapply(1:n_seed, function(t)  round(G$degree[match(mapply(`[[`, first_wave_lsmi[[t]], 1), 1:G$n)]*delta,0) )
  
  first_wave_refs<- c(
                      list(sample(lapply(1:n_seed, 
                                    function(t) mapply(`[[`, first_wave_lsmi[[t]], 2) )[[1]], G_rp[[1]] , FALSE)),
                     
                       lapply(2:n_seed, function(x)
                         lapply(1:length(G_rp[[x]]), function(s)
                           sample( lapply(first_wave_lsmi[[x]], `[[`, 2)[[s]] , G_rp[[x]][[s]], FALSE ) )) )

   seeds = map(lapply(1:n_seed, function(seed) first_wave[[seed]]$sequence_seeds),unlist)
  
  wave_refs_fin<- map(first_wave_refs, unlist)
  
  
  refs_other_waves<- list()
  refs_other_waves[[1]]<- seeds
  refs_other_waves[[2]]<- wave_refs_fin
  
  for (i in 1:(n_waves - 1)) {
    seed_wave <- wave_refs_fin
    
    wave_s <- lapply(1:(n_seed), function(j) {
      if (length(seed_wave[[j]]) == 0) {
        return(list())
      } else {
        lapply(1:length(seed_wave[[j]]), function(h) 
          lsmi_union(G, 
                     n.wave = 1, 
                     seeds = seed_wave[[j]][[h]],
                     n.seeds = length(seed_wave[[j]][[h]])))
      }
    })
    
    wave_s_lsmi <- lapply(1:(n_seed), function(j) {
      if (length(seed_wave[[j]]) == 0) {
        return(list())
      } else {
        lapply(1:length(seed_wave[[j]]), function(h) wave_s[[j]][[h]]$lsmi_big)
      }
    })
    
    G_rp <- lapply(1:(n_seed), function(j) {
      if (length(seed_wave[[j]]) == 0) {
        return(list())
      } else {
        lapply(1:length(seed_wave[[j]]), function(h) 
          round(G$degree[match(lapply(wave_s_lsmi[[j]][[h]], `[[`, 1), 1:G$n)] * delta, 0))
      }
    })
    
    wave_s_refs <- lapply(1:(n_seed), function(j) {
      if (length(seed_wave[[j]]) == 0) {
        return(list())
      } else {
        lapply(1:length(seed_wave[[j]]), function(h)  
          lapply(lapply(wave_s_lsmi[[j]][[h]], `[[`, 2), sample, G_rp[[j]][[h]], FALSE))
      }
    })
    
    wave_refs_fin <- map(wave_s_refs, unlist)
    #adds a placeholder empty list if null
    wave_refs_fin<-lapply(wave_refs_fin, function(x) if (is.null(x)) list() else x)
    
    refs_other_waves[[i + 2]] <- wave_refs_fin
  }
  
  
  
  all_waves_refs<- do.call(Map, c(f = list, refs_other_waves))
  
  return(all_waves_refs)
  
}

#function running my_snowball_fun returning list divided by seeds and in each seeds




# Note about the simulated network:
# studied network is the recruited network which is different from the nomination network
# although we recognize that there is a difference between recruitment and nomination
# we simulate the network of reputation and only study that





#PARAMETERS OF THE SIMULATED NETWORK 
n_nodes = c(50,100,200,400)
n_edges = list(c(1000, 750,  500), 
               c(2000, 1500, 1000),
               c(4000, 3000, 2000),
               c(8000, 6000, 4000) )

n_nodes[1]/n_edges[[1]]
n_nodes[2]/n_edges[[2]]
n_nodes[3]/n_edges[[3]]
n_nodes[4]/n_edges[[4]]

n_ratio = n_nodes[1]/n_edges[[1]]

# Create a function that encapsulates the entire code to be repeated
run_snowball_sampling <- function(G_d_to_net,mean_degree_label) {
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

# ASSUMPTION: recr_perc of the neighbors of a node are recruited
snow_samp_recr = lapply(1:length(n_nodes), function(i)
                    lapply(1:length(n_edges[[i]]),function(j)
                       my_snowball_fun(G_d_to_net[[i]][[j]],
                                       n_seed = n_seed,
                                       lsmi_union, 
                                       n_waves = n_waves , 
                                       delta = delta) )) 




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
                                 length(unique(snow_samp_recr[[i]][[j]][[h]][[k]]) )  ))))





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



gc()
return(list(df_snow_unique = df_snow_unique, df_snow_recr = df_snow_recr))
}





#PARAMETRS snow
n_seed = 5
n_waves = 3
delta = 0.3


# Define the number of networks
num_net<-100
# Define the number of iteration
num_iterations <- 100


network_sampling <- function() {
  G_d<- lapply(1:length(n_nodes), function(i)
    lapply(1:length(n_edges[[i]]),function(j)
      erdos.renyi.game(n_nodes[i], n_edges[[i]][[j]], type="gnm", directed = FALSE, loops = FALSE) ))
  

  
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
  

  
  
  mean_degree_label<- c(unique(mean_degree))
  df_degree$n_nodes = rep(n_nodes, each =3)
  df_degree_melt = reshape2::melt(df_degree, id = c("n_nodes", "mean_degree"))
  
  
  
  
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
  
  
  #####################
  # SNOWBALL SAMPLING #
  #####################
  
  
  #Convert the class of functions to the network package which incorporates the snowball sampling
  G_d_to_net <-  lapply(1:length(n_nodes), function(i)
    lapply(1:length(n_edges[[i]]),function(j) igraph_to_network(G_d[[i]][[j]]) ))
  
  #G_d_to_net[[1]][[3]]
  
  
  start_time <- Sys.time()
  # Iterations for snowball sampling
  for (ni in 1:num_iterations) {
    print(paste("SB Sampling number:", ni))
    result <- run_snowball_sampling(G_d_to_net,mean_degree_label)
    
    if (ni == 1) {
      df_snow_unique <- result$df_snow_unique
      df_snow_recr <- result$df_snow_recr
      df_snow_unique <-df_snow_unique %>% rename(!!paste0("value2_iter_", ni) := value2)
      df_snow_recr <-df_snow_recr %>% rename(!!paste0("value2_iter_", ni) := value2)
    } else {
      df_snow_unique[[paste0("value2_iter_", ni)]] <- result$df_snow_unique$value2
      df_snow_recr[[paste0("value2_iter_", ni)]] <- result$df_snow_recr$value2
    }
  }
  
  end_time <- Sys.time()
  runtime <- end_time - start_time
  print(runtime)
  
  
  
  
  df_snow_unique<-df_snow_unique %>%
    mutate(value2_iter_avg = rowMeans(select(., starts_with("value2_iter_"))))%>%
    mutate(value2=value2_iter_avg) %>% select(., c(-starts_with("value2_iter_"),-value))
  df_snow_recr<-df_snow_recr %>%
    mutate(value2_iter_avg = rowMeans(select(., starts_with("value2_iter_"))))%>%
    mutate(value2=value2_iter_avg) %>% select(., c(-starts_with("value2_iter_"),-value))
  gc()
  return(list(df_snow_unique = df_snow_unique, df_snow_recr = df_snow_recr))
  
} 


#Itearions for differnt networks
for (net_num in 1:num_net) {
  print(paste("Network number:", net_num))
  result_comb <- network_sampling()
    if (net_num == 1) {
    df_snow_unique_comb <- result_comb$df_snow_unique
    df_snow_recr_comb  <- result_comb$df_snow_recr
    df_snow_unique_comb  <-df_snow_unique_comb  %>% rename(!!paste0("value2_net_", net_num) := value2)
    df_snow_recr_comb  <-df_snow_recr_comb  %>% rename(!!paste0("value2_net_", net_num) := value2)
  } else {
    df_snow_unique_comb [[paste0("value2_net_", net_num)]] <- result_comb$df_snow_unique$value2
    df_snow_recr_comb [[paste0("value2_net_", net_num)]] <- result_comb$df_snow_recr$value2
  }
  gc()
}


df_snow_unique_comb<-df_snow_unique_comb %>%
  mutate(value2_net_avg = rowMeans(select(., starts_with("value2_net_"))))%>%
  mutate(value2=value2_net_avg) %>% select(., c(-starts_with("value2_net_")))
df_snow_recr_comb<-df_snow_recr_comb %>%
  mutate(value2_net_avg = rowMeans(select(., starts_with("value2_net_"))))%>%
  mutate(value2=value2_net_avg) %>% select(., c(-starts_with("value2_net_")))



############
# Heatmaps #
############


############
#percentage#
############

#100% recruited (unique)
#gg = ggplot(df_snow_unique, aes(x = factor(variable) , y = factor(n_seed), fill = value2)) 
#gg = gg +  geom_tile(aes(fill = value2))  
#gg = gg +   facet_grid(n_nodes ~ mean_degree )
#gg = gg +  geom_text(aes(label = round(value2, 3)), size = 3) 
#gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
#gg = gg + theme_bw()
#gg = gg + scale_x_discrete(labels= c(1:(n_waves)))
#gg = gg + scale_y_discrete(labels= c(1:n_seed))
#g = gg + xlab("n. waves") + ylab("n. seed")
#gg = gg + labs(fill = "N. Recruited")
#gg = gg + scale_fill_gradient(low = "white", high = "steelblue") 
#gg


#recr_perc% recruited 
gg = ggplot(df_snow_recr_comb, aes(x = factor(variable) , y = factor(n_seed), fill = value2)) 
gg = gg +  geom_tile(aes(fill = value2))  
gg = gg +   facet_grid(n_nodes ~ mean_degree )
gg = gg +  geom_text(aes(label = round(value2, 2)), size = 3) 
gg = gg +  scale_color_manual(guide = "none", values = c("gray", "black"))
gg = gg + theme_bw()
gg = gg + scale_x_discrete(labels= c(1:(n_waves)))
gg = gg + scale_y_discrete(labels= c(1:n_seed))
gg = gg + xlab("n. waves") + ylab("n. seed")
gg = gg + labs(fill = "Perc. Recruited")
gg = gg + scale_fill_gradient(low = "white", high = "steelblue") 
gg



ggsave(filename =  paste(my_dir_figs,"e_r_snow_0.3_recr.pdf", sep = ""),
       plot = gg)





# Define constants
Z <- 1.96  # Z-score for 95% confidence
E <- 0.03   # Margin of error (3%/10%)

# Calculate sample size and add it as a new column
df_snow_recr_comb<- df_snow_recr_comb %>%
  mutate(
    sample_size = ceiling((Z^2 * value2 * (1 - value2)) / E^2)
  )
df_snow_recr_comb<-df_snow_recr_comb %>%
  mutate(highlight = ifelse(sample_size <= n_nodes & value2>=0.50, "Highlight", "Normal")) 

#Heatmaps with cases satisfying the criteria highlighted
ggplot(df_snow_recr_comb, aes(x = factor(variable), y = factor(n_seed))) +
  geom_tile(aes(fill = value2), color = "white", size = 0.5) +  # Fill tiles with proportion value2, with a white border
  facet_grid(n_nodes ~ mean_degree) +  # Facet by n_nodes and mean_degree
  geom_text(aes(label = round(value2, 2), 
                fontface = ifelse(highlight == "Highlight", "bold", "plain")),  # Make text bold if highlighted
            size = 3) +  # Add text labels for value2
  scale_fill_gradient(low = "white", high = "steelblue") +  # Gradient for value2
  theme_bw() +  # Use black-and-white theme
  scale_x_discrete(labels = c(1:n_waves)) +  # X-axis labels
  scale_y_discrete(labels = c(1:n_seed)) +  # Y-axis labels
  xlab("Number of Waves") + ylab("Seed Number") +  # Axis labels
  labs(fill = "Percentage Recruited")  # Fill legend title




#Line graphs
df_with_stats <- df_snow_recr_comb %>%
  group_by(n_nodes,mean_degree,n_seed, n_edges, variable,L1)
# Plot the results
ggplot(df_with_stats, aes(x = variable, y = value2, color = as.factor(n_seed), group = n_seed)) +
  geom_line() + # Connect the points with lines for each seed
  geom_point() + # Add points for each seed
  facet_grid(n_nodes ~ n_edges) + # Facet by both n_nodes and n_edges
  labs(
    title = "Recruited percentage by Wave Number and Network Configuration",
    x = "Wave Number",
    y = "Recruited %",
    color = "Seed Number",
    fill = "Seed Number"
  ) +
  theme_minimal() # Use a clean theme


ggplot(df_with_stats, aes(x = as.factor(n_seed), y = value2, color = as.factor(variable), group = variable)) +
  geom_line() + # Connect the points with lines for each variable
  geom_point() + # Add points for each seed
  facet_grid(n_nodes ~ n_edges) + # Facet by both n_nodes and n_edges
  labs(
    title = "Recruited percentage by Seed Number and Network Configuration",
    x = "Seed Number",
    y = "Recruited %",
    color = "Wave Number",
    fill = "Wave Number"
  ) +
  theme_minimal()

############

























