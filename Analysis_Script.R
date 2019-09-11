# _______________ Reorganization of surviving mammal communities ________________#
#________________ after end-Pleistocene megafaunal extinction ________________ #

#________________ Script to replicate analyses   ________________ #
#________________     by Anikó B. Tóth    ________________ #

## This script runs all of the analyses presented in the manuscript and supplement.
## The script allows a user to choose certain parameters, 
  # such as minimum site richness (co) and subsample size (ss1), among others.

######## Environment ################
if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

p <- c("combinat", "reshape2", "sp", "hypervolume", 
       "cowplot", "grid", "gtable", "EcoSimR")

sapply(p, require, character.only = TRUE) %>% `[`(!.) %>% names %>% map(install.packages)
sapply(p, require, character.only = TRUE)

source('Helper_functions.R')

#### Load Data ####
load("Data/FML_dated.RData") # Faunmap locality data with calibrated calendar dates added
load("Data/FM_final.RData") # Faunmap occurrence data (Pleistocene and Holocene) with marine and indeterminate species removed
load("Data/MOM4.1_NewWorld.RData") # MOM database with body masses for mammal species

#### Prepare data ####
# merge occurrence data with relevant site data
x <- merge(fm, fml[,c('id', "DepositionalSystem", "LATDD","LONGDD", 'Epoch', "MinAge", "MaxAge", "MeanAge")], by.x = 'id', by.y = "id", all = T)
x <- merge(x, mom.nw[,c("taxon", "LogMass")], by.x = "name", by.y = "taxon", all.x = T, sort = F)

## Filters 
md <- 21000 # maximum date included (note that climate is not available for sites older than 21ka)
mbd <- 1000 # minimum body mass in grams
mlat <- 60  # maximum latitude
hr <- 2000  # Holocene-Recent cutoff
co <- 3     # minimum number of species per site
occ <- 1     # minimum number of occurrences per species

x <- filter(x, MeanAge <= md &           # filter out sites older than md
              LogMass >= log10(mbd) &    # filter out species smaller than mbd
              LATDD <= mlat &            # filter sites north of 60N
              !Epoch == "HOPL" &         # remove sites without definite epoch classification
              !Epoch == "")              # removes sites with no epoch

x[which(x$Epoch == "LPLE"),]$Epoch <- "PLEI"  #combine all Pleistocene sites into one category
x[which(x$MaxAge <= hr),]$Epoch <- "MOD"      # Implement cutoff between Holocene and Recent

PA <- split(x, x$Epoch) %>% 
  purrr::map(dcast, formula = id~name, value.var = 'observed', fun.aggregate = length) %>%
  purrr::map(namerows) %>% tobinary() %>% purrr::map(as.data.frame) %>% purrr::map(base::t) #structure data

PA <- purrr::map(PA, ~return(.[,which(colSums(.) >= co)]))    # remove sites with less than minimum species
PA <- purrr::map(PA, ~return(.[which(rowSums(.) >= occ),]))   # remove species with no occurrences

PA <- PA[c("MOD", "HOLO", "PLEI")]  # put time intervals in order.

### Niche analysis ####
  #### Prepare data ####
  load("Data/sitedat_clim1-4k_mm1-5.Rdata") # this file contains climate data for each site based on various averages.
  
  # sd.clim1k: 1ky average based on site mean age
  # sd.clim2k to 4k: multi thousand year average based on mean age and the next youngest one to three layers
  # sd.mm1 - 5: 1ky average based on the minimum or maximum age of sites, chosen randomly. Five random iterations are included. 
  
  # Results presented in the main text use sd.clim1k. 
  # Results presented in the supplement use sd.clim4k and sd.mm1
  
  # choose which climate data version to work with
  fml.sitedat <- sd.clim1k 
    #fml.sitedat <- sd.clim4k  # dating sensitivity results presented in supplement
    #fml.sitedat <- sd.mm1 # dating sensitivity results presented in supplement
  
  # categorize species and apirs
    survivors <- c(unique(c(rownames(PA[[1]]), rownames(PA[[2]]))), "Ovibos_moschatus") # under default parameters, the musk ox is the only extant species that was not sampled in the Holocene and Recent intervals
    extinct <- rownames(PA[[3]])[which(!rownames(PA[[3]]) %in% survivors)]
    
    # survivor-survivor pairs
    SS <- t(combn(survivors, 2))
    SS <- data.frame(pair = c(map2_chr(SS[,1], SS[,2], paste, sep = "-"), map2_chr(SS[,2], SS[,1], paste, sep = "-"))) 
     
    # extinct-extinct pairs
    EE <- t(combn(extinct, 2))
    EE <- data.frame(pair = c(map2_chr(EE[,1], EE[,2], paste, sep = "-"), map2_chr(EE[,2], EE[,1], paste, sep = "-")) )
     
    # survivor-extinct pairs
    SE <- t(combn(c(survivors, extinct), 2))
    SE <-data.frame(pair = c(map2_chr(SE[,1], SE[,2], paste, sep = "-"), map2_chr(SE[,2], SE[,1], paste, sep = "-")))
    SE <- data.frame(pair = SE[!SE$pair %in% SS$pair & !SE$pair %in% EE$pair,])
    
    pairs <- bind_rows(list(SS=SS, EE=EE, SE=SE), .id = "type")
    
  #### Return PA data to long format ####
    sitebyspecies <-  map(PA, as.matrix) %>% map(melt) %>% bind_rows(.id = "tbn") %>% filter(value == 1) %>% setNames(c("tbn", "name", "id", "value"))
    # Add climate data 
    sitebyspecies <- merge(sitebyspecies, fml.sitedat[,c("id","LATDD", "LONGDD", "DepositionalSystem", "MinAge", "MaxAge", "MeanAge", "MAP", "MAT", "bio4", "bio15")], by = "id", all = T) %>% na.omit()
    # Add species status data
    sitebyspecies$status[sitebyspecies$name %in% survivors] <- "survivor"
    sitebyspecies$status[sitebyspecies$name %in% extinct] <- "victim"

    #Change coordinates to equal area
    equal.area <- SpatialPoints(sitebyspecies[,c("LONGDD", "LATDD")], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>%
      spTransform(CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")) %>%
      data.frame() %>% setNames(c("laea.long", "laea.lat"))
    
    sitebyspecies <- cbind(sitebyspecies, equal.area)
    
    # Square root transform precipitation variables
    sitebyspecies$MAP <- sqrt(sitebyspecies$MAP)
    sitebyspecies$bio15 <- sqrt(sitebyspecies$bio15)
    
    #z-transform all abiotic variables
    sitebyspecies.z <- sitebyspecies
    abvar <- c("MAP", "MAT", "bio4", "bio15", "laea.long", "laea.lat")
    sitebyspecies.z[,abvar] <- apply(sitebyspecies.z[,abvar], 2, scale)
    
  #### site table
    sites <- sitebyspecies.z %>% dplyr::select(id, tbn, laea.lat, laea.long, MAP, MAT, bio4, bio15, MeanAge, MinAge, MaxAge) %>% unique()
  
  #### Calculate separate geographic and climatic Niche and background volumes/areas  #####
    # OPTION 1: with hypervolumes ####
        # **THIS TAKES A LONG TIME TO RUN**  
        # source("Hypervolume_script.R")
        # it can also be partially run by a careful user -- intermediate named objects are recycled throughout.
        # There is a block of code at the bottom that produces the desired product v even if you only have some of the hypervolumes calculated.
    # OPTION 2: with convex hulls ####
          # total area by timebin
        ch_geog <- sitebyspecies.z %>% split(.$tbn) %>% map(find_hull.geog) %>% map(select, laea.long, laea.lat) %>% map(Polygon, hole=F) %>% map_dbl(function(x) x@area)
        ch_clim <- sitebyspecies.z %>% split(.$tbn) %>%  map(find_hull.clim) %>% map(select, MAP, MAT) %>% map(Polygon, hole=F) %>% map_dbl(function(x) x@area)
        ch_seas <- sitebyspecies.z %>% split(.$tbn) %>% map(find_hull.seas) %>% map(select, bio4, bio15) %>% map(Polygon, hole=F) %>% map_dbl(function(x) x@area)
        
          # niche area and proportion of total area by species
        ssv_geog <- sitebyspecies.z %>% split(interaction(.$tbn, .$status)) %>% 
          purrr::map(niche_areas, type = "geog") %>% bind_rows(.id = 'group')
        ssv_geog <- strsplit(ssv_geog$group, fixed = TRUE, split = ".") %>% 
          reduce(rbind) %>% data.frame(ssv_geog) %>% setNames(c("tbn", "status", "group", "name", "chull.area"))
        ssv_geog$fill = ssv_geog$chull.area / rep(ch_geog, times = table(ssv_geog$tbn))
        
        ssv_clim <- sitebyspecies.z %>% split(interaction(.$tbn, .$status)) %>% 
          purrr::map(niche_areas, type = "climate") %>% bind_rows(.id = 'group')
        ssv_clim <- strsplit(ssv_clim$group, fixed = TRUE, split = ".") %>% 
          reduce(rbind) %>% data.frame(ssv_clim) %>% setNames(c("tbn", "status", "group", "name", "chull.area"))
        ssv_clim$fill = ssv_clim$chull.area / rep(ch_clim, times = table(ssv_clim$tbn))
        
        ssv_seas <- sitebyspecies.z %>% split(interaction(.$tbn, .$status)) %>% 
          purrr::map(niche_areas, type = "seas") %>% bind_rows(.id = 'group')
        ssv_seas <- strsplit(ssv_seas$group, fixed = TRUE, split = ".") %>% 
          reduce(rbind) %>% data.frame(ssv_seas) %>% setNames(c("tbn", "status", "group", "name", "chull.area"))
        ssv_seas$fill = ssv_seas$chull.area / rep(ch_seas, times = table(ssv_seas$tbn))
    
    
  # Calculate overlap of hypervolumes with sites (use only ONE of the two options below) ####
    # OPTION 1: with hypervolumes ####
      # by species (across timebins)
          niche.geog <- get_niche_table_hypervolume(sitebyspecies.z, dims = c("laea.lat", "laea.long"), min.occ = 4) 
          niche.clim <- get_niche_table_hypervolume(sitebyspecies.z, dims = c("MAP", "MAT", "bio4", "bio15"), min.occ = 4)
          niche <- niche.geog[sort(rownames(niche.geog)),] + niche.clim[sort(rownames(niche.clim)),]
          niche[niche==1] <- 0
          niche[niche==2] <- 1

      # by species by timebin
          niche.geog.sp.tbn <- sitebyspecies.z %>% split(.$tbn) %>% map(get_niche_table_hypervolume, dims = c("laea.lat", "laea.long"), min.occ = 4)
          niche.geog.sp.tbn <- niche.geog.sp.tbn[c("MOD", "HOLO", "PLEI")]
    # OPTION 2: with convex hulls ####
        niche <- get_niche_table_6D(sitebyspecies)

### Co-occurrence analysis with subsampling #### 
  # not presented in paper except in Fig. S9 to contrast randomization tests
  rps = 100 #number of subsamples
  ss1 = 60 # subsample size, must be less than minimum number of sites: min(map_int(PA, ncol))
  
  nm <- map(PA, resamp, reps = rps, sites = ss1, replace = F) # subsamples (nested list)
  nz <- map(nm, map, simpairs) # FETmP on subsamples
  nel <- map(1:length(nz), function(x) map2(nz[[x]], nm[[x]], dist2edgelist)) %>%   # change to data frame from dist objects
    map(bind_rows, .id = "subsample") %>%              # consolidate subsample results
    setNames(c("Recent", "Holocene", "Pleistocene")) %>%             # name output list
    bind_rows(.id = "timebin")                         # consolidate timebins to one data frame.
  
### Biotic / Abiotic co-occurrence analysis #####
  # Set parameters
  ss1 <- 60      # size of subsample for each repetition, must be less than min number of sites: min(map_int(PA, ncol))
  ss2 <- 10      # minimum number of mutual niche sites to be included (Mij)
  reps <- 100   # number of iterations
  
  # Calculate biotic/abiotic components
  out <- list()
  for(i in 1:reps){
    PAs <- lapply(PA, function(x) return(x[,sample(1:ncol(x), ss1, replace = F)]))
    PAs <- lapply(PAs, function(x) return(x[which(rowSums(x) > 0),]))  # remove any empty rows
    y <- lapply(PAs, t) # transpose
    out[[i]] <- abio_single_niche(y, niche, ss2)
    print(i) # to keep track of progress
  }
    
  names(out) <- c(1:reps)  
  out <- bind_rows(out, .id = "subsample")
  out$type <- pairs$type[match(out$id, pairs$pair)]
    
#### Randomization analyses #####
  library(EcoSimR)
  Sim2 <- map(nm, map, sim2)              # Fixed row (species occupancy), equiprobable column (site richness) sums
  zSim2 <- map(Sim2, map, simpairs)       # FETmP scores for fixed-equiprobable (FE) randomization
  zSim2el <- map2(zSim2, nm, map2, dist2edgelist) %>%     # change dist to dataframe
    setNames(c("Recent", "Holocene", "Pleistocene")) %>%  # name output list
    flat2(layer1 = "timebin", layer2 = "subsample")       # flatten, specify what the nested layers are. (flat2 is in helper functions)
  
  Sim3 <- lapply(nm, lapply, sim3)        # Equiprobable row (species occupancy), fixed column (site richness) sums
  zSim3 <- lapply(Sim3, lapply, simpairs) # FETmP scores for equiprobable-fixed (EF) randomization
  zSim3el <- map2(zSim3, nm, map2, dist2edgelist) %>%     # change dist 2 data frame
    setNames(c("Recent", "Holocene", "Pleistocene")) %>%  # name output list
    flat2(layer1 = "timebin", "subsample")                # flatten, specify what the nested layers are. (flat2 is in helper functions)
  
  d <- reduce(list(nel, zSim2el, zSim3el), merge, all = T, by = c("timebin", "subsample", "id")) %>%
    dplyr::select(timebin, subsample, id, Z.Score.x, Z.Score.y, Z.Score) %>% 
    setNames(c("timebin", "subsample", "id", "Z.Score", "FE", "EF"))
  
  
  ##########

    
    
    