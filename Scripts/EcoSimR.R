library(EcoSimR)


set.seed(1)

# ---- Load in the data


years <- c(2010, 2011, 2012, 2013, 2015, 2016, 2017, 2019, 2020, 2021) # List of years where we have wet season surveys

da.wet.list <- lapply(years, function(y) {
  read.csv(paste0("Processed Data/BaleMountainsNP_MeanAbundance_", y, "_WetSeason_res100m.csv"))
}) # Create list of all dataframes for each years



# Not all dataframes have records for mixed livestock, add an empty column if its absent so we can sum all domestic livestock together

livestock_names <- c("Cattle", "Sheep...goats", "Horse", "Donkey", "Mixed.Livestock") 
da.wet.list <- lapply(da.wet.list, function(df) {
  for (var in livestock_names) {
    if (!var %in% names(df)) {
      df[[var]] <- 0
    }
  }
  return(df)
})


da.wet.list <-lapply(da.wet.list, function(df) {
  df %>% mutate(Livestock=Cattle+Sheep...goats+Horse+Donkey+Mixed.Livestock)
}
)

# Create dataframe compatible with EcoSimR

Y.wet.list <- lapply(da.wet.list, function(df) {
  df_1<-data.frame(df %>% dplyr::select(id,Mountain.Nyala,Bohor.Reedbuck, Menelik.s.Bushbuck, Warthog, Livestock))
  df_final <- data.frame(df_1%>% dplyr::select(-1))
  rownames(df_final) <- df_1$id
  return(df_final)
}
) # Select dependent variables, i.e. abundance data for 4 most important wild ungulate species plus livestock
# Rownames are site names




# We create a dataframe that is purely presence-absence  
Y.wet.presab.list <- Y.wet.list 
Y.wet.presab.list <- lapply(Y.wet.presab.list, function(df) { 
  df[] <- lapply(df, function(x) ifelse(x > 0, 1, x))  # Replace values greater than 0 with 1
  return(df)
})

Y.wet.presab.list[[1]]

# Creat EcoSimR by transposing sites and species
EcoSIM.data <- lapply(Y.wet.presab.list, function(df) {df
  transposed.data<-as.data.frame(t(df)) })

cooc.model.list<-lapply(EcoSIM.data, cooc_null_model, algo = "sim2", metric = "c_score")

summary(cooc.model.list[[1]])

plot(cooc.model.list[[1]],type="hist") # All years show clear evidence of
plot(cooc.model.list[[10]],type="cooc") # All years show clear evidence of

sim.c_score<-as.data.frame(cooc.model.list[[1]]$Sim)
sim.c_score <- sim.c_score %>% rename(value = `cooc.model.list[[1]]$Sim`)

observed.c_score <-  cooc.model.list[[1]]$Obs
upper.quantile<- quantile(sim.c_score$value, probs = 0.975)
lower.quantile<- quantile(sim.c_score$value, probs = 0.025)


ggplot(sim.c_score, aes(x = value)) +
  geom_histogram(binwidth = 25, fill = "skyblue", color = "black") +
  geom_vline(xintercept = upper.quantile, color = "black", linetype= "dashed", size = 1) +
  geom_vline(xintercept = lower.quantile, color = "black", linetype= "dashed", size = 1) +
  geom_vline(xintercept = observed.c_score, color = "red", size = 1) +
  theme_bw() +
  annotate("text", x = (upper.quantile+lower.quantile)/2, y = 62, label = "Simulated", color = "skyblue", angle = 0, size = 5) +
  annotate("text", x = upper.quantile + 210, y = 25, label = "Upper 97.5%", color = "black", angle = 0, size = 4) +
  annotate("text", x = lower.quantile - 210, y = 25, label = "Lower 2.5%", color = "black", angle = 0, size = 4) +
  annotate("text", x = observed.c_score + 180, y = 25, label = "Observed", color = "red", angle = 0, size = 4) +
  xlim(min(sim.c_score$value-100), max(observed.c_score+500)) +
  labs(title = "Observed vs Simulated C-Scores - 2010",
       x = "C-Score", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
  
