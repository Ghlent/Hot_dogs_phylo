# Set working directory
setwd("Phylogenetics/Hot_dogs")

# Import libraries
library(ape)
library(igraph)

# Prepare data
data <- read.csv2("Data/Character matrix.csv", row.names = 1)
row.names(data) <- gsub(" ", "_", row.names(data))
str(data)

# Separate data into components
names <- row.names(data)
places <- data$Place
years <- data$Year
traits <- data[-c(1,2)]
traits <- apply(traits, 2, function(x) gsub(" \\| ", " ", x))
colnames(traits)
unique(traits[,8])

# Save years for dating
years_df <- data.frame(
    taxon = paste0(noquote(names), "_fossil"),
    min_age = 2025 - as.numeric(sub(".*-", "", years)),
    max_age = 2025 - as.numeric(sub("-.*", "", years))
)
years_df <- years_df[!is.na(years_df$max_age),]
years_df[years_df$min_age == years_df$max_age,"max_age"] <- years_df[years_df$min_age == years_df$max_age,"max_age"]+0.9
write.table(years_df, file="Data/fossils.csv", sep=",", row.names = FALSE, quote=FALSE)

# Save traits as nexus file, including a copy of the fossil of the hot dog species
traits <- apply(traits, c(1,2), function(x) {
    parts <- unlist(strsplit(as.character(x), " "))
    if (length(parts) == 1) {
      return(parts)  # Single state, keep as-is
    } else {
      return(paste0("(", paste0(parts, collapse=""), ")"))  # Polymorphic: wrap in ()
    }
  })
write_polymorphic_nexus <- function(data, filename, symbols = "01234567") {
  taxa <- rownames(data)
  nchars <- ncol(data)
  ntax <- nrow(data)

  f <- file(filename, "w")

  writeLines("#NEXUS", f)
  writeLines("BEGIN DATA;", f)
  writeLines(paste("  DIMENSIONS NTAX=", ntax+nrow(years_df), " NCHAR=", nchars, ";", sep=""), f)
  writeLines(paste0('  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS="', symbols, '" POLYMORPH=();'), f)
  writeLines("  MATRIX", f)

  for (i in 1:ntax) {
    line <- paste0(taxa[i], " ", paste0(data[i, ], collapse=""), sep="")
    writeLines(line, f)
    if(paste0(taxa[i], "_fossil") %in% years_df$taxon) {
        line <- paste(paste0(taxa[i], "_fossil "), paste0(data[i, ], collapse=""), sep="")
        writeLines(line, f)
    }
  }

  writeLines("  ;", f)
  writeLines("END;", f)
  close(f)
}
write_polymorphic_nexus(traits, file = "Data/morph_traits.nex")

# Pre-calculate distances between any two points
places_coordinates <- list(
  c(40.7128, -74.0060),
  c(40.7128, -74.0060),
  c(40.7357, -74.1724),
  c(38.9072, -77.0369),
  c(39.9526, -75.1652),
  c(41.5801, -71.4774),
  c(42.3601, -71.0589),
  c(39.9526, -75.1652),
  c(39.2904, -76.6122),
  c(40.8584, -74.1638),
  c(40.9168, -74.1718),
  c(39.9526, -75.1652),
  c(39.2904, -76.6122),
  c(45.2538, -69.4455),
  c(43.0000, -75.0000),
  c(43.0000, -75.0000),
  c(32.75, -86.67),
  c(32.460976, -84.987709),
  c(31.0000, -100.0000),
  c(35.1495, -90.0490),
  c(41.8781, -87.6298),
  c(40.7128, -74.0060),
  c(40.56, -82.05),
  c(41.4993, -81.6944),
  c(41.8781, -87.6298),
  c(42.3314, -83.0458),
  c(43.0125, -83.6875),
  c(32.2988, -90.1848),
  c(39.1031, -84.5120),
  c(46.7296, -94.6859),
  c(34.0489, -111.0937),
  c(36.1699, -115.1398),
  c(47.6062, -122.3321),
  c(39.7392, -104.9903),
  c(36.7783, -119.4179),
  c(44.0682, -114.7420),
  c(34.0522, -118.2437),
  c(21.3069, -157.8583),
  c(64.2008, -149.4937),
  c(39.5501, -105.7821)
)
places_coordinates <- setNames(places_coordinates, places)
haversine <- function(lat1, lon1, lat2, lon2) {
  # Convert degrees to radians
  to_rad <- function(deg) deg * pi / 180
  
  lat1 <- to_rad(lat1)
  lon1 <- to_rad(lon1)
  lat2 <- to_rad(lat2)
  lon2 <- to_rad(lon2)
  
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  
  R <- 6371  # Earth's radius in km
  d <- R * c
  return(d)
}
physical_dists <- c()
for (i in 1:(length(places)-1)) {
    for (j in (i+1):length(places)) {
        physical_dists <- append(physical_dists, haversine(places_coordinates[[i]][1], places_coordinates[[i]][2], places_coordinates[[j]][1], places_coordinates[[j]][2]))
    }
}

# Import selected tree
tree <- read.nexus("rev_Output/mk_gamma.mcc.tre")
no_fossils_tree <- drop.tip(tree, paste0(names, "_fossil"))

# Create function for calculating the geographical signal score
gs_score <- function(tree) {    
    g <- as.igraph(tree)
    topo_dist_matrix <- distances(g)[names,names]
    topo_dists <- c()
    for (i in 1:(length(places)-1)) {
        for (j in (i+1):length(places)) {
            topo_dists <- append(topo_dists, topo_dist_matrix[names[i], names[j]])
        }
    }
    sum(log10(physical_dists+10) * topo_dists)
}

# Permutation test
seed(1)
obs_score <- distance_score(no_fossils_tree)
permutations <- 100000
results <- c()
for (i in 1:permutations) {
    r_tree <- no_fossils_tree
    r_tree$tip.label <- sample(r_tree$tip.label)
    results[i] <- distance_score(r_tree)
}
results <- append(results, obs_score)
sum(results >= obs_score)/length(results)

# Save distribution of the permutation test as figure
svg("Figures/gs_scores.svg")
hist(results, main="Distribution of simulated GS Scores", xlab = "GS Scores", xlim=c(22650,23150))
abline(v = quantile(results, 0.95), col = "black", lwd=3)
abline(v = obs_score, col = "red", lwd=3)
legend(x="topleft", legend=c("95th percentile", "Observed score"), fill=c("black", "red"))
dev.off()
