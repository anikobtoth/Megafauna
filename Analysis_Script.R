# _______________ End-Pleistocene extinction caused a fundamental ________________#
#________________ shift in mammal survivor community structure ________________ #

#________________ Script to replicate analyses   ________________ #
#________________     by Anikó B. Tóth    ________________ #
#________________    November 8, 2018    ________________ #

## This script runs all of the analyses presented in the manuscript and supplement.
## The script allows a user to choose certain parameters, 
  # such as minimum site richness (co) and subsample size (ss1), among others.

######## Environment ################
library(tidyverse)
library(reshape2)

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
co <- 5     # minimum number of species per site
occ <- 1     # minimum number of occurrences per species

x <- filter(x, MeanAge <= md &           # filter out sites older than md
              LogMass >= log10(mbd) &    # filter out species smaller than mbd
              LATDD <= mlat &            # filter sites north of 60N
              !Epoch == "HOPL" &         # remove sites without definite epoch classification
              !Epoch == "")              # removes sites with no epoch

x[which(x$Epoch == "LPLE"),]$Epoch <- "PLEI"  #combine all Pleistocene sites into one category
x[which(x$MaxAge <= hr),]$Epoch <- "MOD"      # Implement cutoff between Holocene and Recent

PA <- split(x, x$Epoch) %>% 
  map(dcast, formula = id~name, value.var = 'observed', fun.aggregate = length) %>%
  map(namerows) %>% tobinary %>% map(as.data.frame) %>% map(t) #structure data

PA <- map(PA, ~return(.[,which(colSums(.) >= co)]))    # remove sites with less than minimum species
PA <- map(PA, ~return(.[which(rowSums(.) >= occ),]))   # remove species with no occurrences

PA <- PA[c("MOD", "HOLO", "PLEI")]  # put time intervals in order.

### Niche analysis ####
  #### Prepare data ####
  load("Data/sitedat_clim1-4k_mm1-5.Rdata") # this file contains climate data for each site based on various averages.
  
  # sd.clim1k: 1ky average based on site mean age
  # sd.clim2k to 4k: multi thousand year average based on mean age and the next youngest one to three 1k layer
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
    SS <- data.frame(pair = map2_chr(t(combn(survivors, 2))[,1], t(combn(survivors, 2))[,2], paste, sep = "-")) 
    
    # extinct-extinct pairs
    EE <- data.frame(pair = map2_chr(t(combn(extinct, 2))[,1], t(combn(extinct, 2))[,2], paste, sep = "-")) 
    
    # survivor-extinct pairs
    SE <-data.frame(pair = map2_chr(t(combn(c(survivors, extinct), 2))[,1], t(combn(c(survivors, extinct), 2))[,2], paste, sep = "-")) 
    SE <- SE[!SE %in% SS & !SE %in% EE]
    
    pairs <- bind_rows(list(SS=SS, EE=EE, SE=SE), .id = "type")
  
  #### Return PA data to long format ####
    sitebyspecies <-  map(PA, as.matrix) %>% map(melt) %>% bind_rows(.id = "tbn") %>% filter(value == 1)
    names(sitebyspecies) <- c("tbn", "name", "id", "value")
    # Add climate data 
    sitebyspecies <- merge(sitebyspecies, fml.sitedat[,c("id","LATDD", "LONGDD", "DepositionalSystem", "MinAge", "MaxAge", "MeanAge", "MAP", "MAT", "bio15")], by = "id", all = T) %>% na.omit()
    # Add species status data
    sitebyspecies$status[sitebyspecies$name %in% survivors] <- "survivor"
    sitebyspecies$status[sitebyspecies$name %in% extinct] <- "victim"
    sitebyspecies$name <- as.character(sitebyspecies$name)
    # change timebin to factor
    sitebyspecies$tbn <- factor(sitebyspecies$tbn, levels = c("MOD", "HOLO", "PLEI"))


  #### Calculate Niche areas #####
    sv_geog <- sitebyspecies %>% filter(tbn == "PLEI") %>% split(.$status) %>% 
      purrr::map(niche_areas, typ = "geog") %>% bind_rows(.id = 'status')
      sv_geog$tbn <- "PLEI"
    sv_clim <- sitebyspecies %>% filter(tbn == "PLEI") %>% split(.$status) %>% 
      purrr::map(niche_areas, typ = "climate") %>% bind_rows(.id = 'status')
      sv_clim$tbn <- "PLEI"
    ss_geog <- sitebyspecies %>% filter(status == "survivor") %>% split(.$tbn) %>% 
      purrr::map(niche_areas, type = "geog") %>% bind_rows(.id = 'tbn')
      ss_geog$tbn <- factor(ss_geog$tbn, levels = c("MOD", "HOLO", "PLEI"))
      ss_geog$status <- "survivor"
    ss_clim <- sitebyspecies %>% filter(status == "survivor") %>% split(.$tbn) %>% 
      purrr::map(niche_areas, type = "climate") %>% bind_rows(.id = 'tbn')
      ss_clim$tbn <- factor(ss_clim$tbn, levels = c("MOD", "HOLO", "PLEI"))
      ss_clim$status <- "survivor"

 
### Simple co-occurrence analysis with subsampling #### 
  # not presented in paper except in Fig. S6 to contrast randomization tests
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
  reps <- 1000   # number of iterations
  
  # Calculate combined niche sites (Ni) 
  niche <- get_niche_table(sitebyspecies)
  
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

    
    
    