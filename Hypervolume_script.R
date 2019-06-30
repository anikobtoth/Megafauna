# _______________ Reorganization of surviving mammal communities ________________#
#________________ after end-Pleistocene megafaunal extinction ________________ #

#________________ Script to create hypervolumes   ________________ #
#________________     by Anikó B. Tóth    ________________ #


# This script can be used to generate hypervolume calculations for species in PA (presence-absence table) in various different ways. 
# The hypervolumes are used to calculate geographic and climate envelopes for species occurrences and backgrounds. 

library(hypervolume)

# Terminology: 
# tbn timebins: MOD (Recent), HOLO (Holocene), PLEI (end-Pleistocene)
# sp species  : list of species names across the entire dataset
# hvm         : hypervolume


# GEOGRAPHIC HYPERVOLUMES ######
  # a: Geo by tbn     : hypervolumes to estimate the total amount of geographic space available in each timebin
      hvm.geog.tbn <- sites %>% split(.$tbn) %>% map(select,laea.lat, laea.long) %>% map(hypervolume_svm)
      
  # b: Geo by sp      : hypervolumes to estimate the total amount of geographic space a species could access in the entire time
      hvm.geog.sp <- sitebyspecies.z %>% get_hypervolume(dims = c("laea.lat", "laea.long"))
      
  # c: Geog by sp, tbn: hypervolumes to estimate the amount of geographic space a species covered in each timebin. 
      hvm.geog.sp.tbn <- sitebyspecies.z %>% split(.$tbn) %>% map(get_hypervolume, dims = c("laea.lat", "laea.long"))
      

# CLIMATIC HYPERVOLUMES #######
  # d: Clim by tbn    : hypervolumes estimate the total amount of climatic space available in each timebin
      hvm.clim.tbn <- sites %>% split(.$tbn) %>% map(select, MAP, MAT, bio4, bio15) %>% map(hypervolume_svm) 
      
  # e: Clim by sp     : hypervolumes estimate the total amount of climatic space a species could survive in during the entire time
      hvm.clim.sp <- sitebyspecies.z %>% get_hypervolume(dims = c("MAP", "MAT", "bio4", "bio15"))
      
  # f: Clim by sp, tbn: hypervolumes estimate the amount of climatic space a species survived in during each timebin. 
      hvm.clim.sp.tbn <- sitebyspecies.z %>% split(.$tbn) %>% map(get_hypervolume, dims = c("MAP", "MAT", "bio4", "bio15"))

  
# CLIMATIC HYPERVOLUMES BY GEOG HYPERVOLUME INCLUSION ####      
  # g: Clim by geo by tbn    : Hypervolumes estimate the climatic space available in geographic space in each timebin 
          # actually same as d, as all sites in a timebin should be included in the geographic hypervolume. 
      hvm.clim.geog.tbn <- sites %>% split(.$tbn) %>% map(select, MAP, MAT, bio4, bio15) %>% map(hypervolume_svm) 
      
  # h: Clim by geo by sp     : Hypervolumes estimate the climatic space available in the geographic space a species could access in the entire time.
          # This is uniqe from e because it automatically excludes geographically inaccessible sites, even if they are climatically suitable. 
      # inclusion test    
      niche.geog.sp <- get_niche_table_hypervolume(sitebyspecies.z, sites = sites, dims = c("laea.lat", "laea.long"), min.occ = 4) 
      s <- apply(niche.geog.sp, 1, function(x) names(x[x==1]))
      
      # hypervolume
      hvm.clim.geog.sp <- purrr::map(s, function(x) return(sites %>% filter(id %in% x))) %>% bind_rows(.id = "name") %>% 
        get_hypervolume(dims = c("MAP", "MAT", "bio4", "bio15"))
      
  # i: Clim by geo by sp, tbn: Hypervolumes estimate the climatic space available in the geographic space a species could access during each timebin. 
          # This is unique from f because it automatically excludes geographically inaccessible sites, even if they are climatically suitable.
      # inclusion test
      niche.geog.sp.tbn <- sitebyspecies.z %>% split(.$tbn) %>% map(get_niche_table_hypervolume, dims = c("laea.lat", "laea.long"), min.occ = 4)
      niche.geog.sp.tbn <- niche.geog.sp.tbn[c("MOD", "HOLO", "PLEI")]
      s <- map(niche.geog.sp.tbn, apply, 1, function(x) names(x[x==1]))
      
      # hypervolume
      hvm.clim.geog.sp.tbn <- purrr::map(s, map, function(x) return(sites %>% filter(id %in% x))) %>% map(bind_rows, .id = "name") %>% 
        map(get_hypervolume, dims = c("MAP", "MAT", "bio4", "bio15"))
      
# CLIMATIC HYPERVOLUMES ON ABSENCES BY GEOG HYPERVOLUME INCLUSION ####          
  # j: Absent clim by sp by tbn     : hypervolume estimates the climatic space in each timebin from which focal species is absent.
      # get absences
      s.abs <- map(PA, apply, 1, function(x) names(x[x==0]))
      # hypervolume
      hvm.clim.sp.tbn.abs <- map(s.abs, map, function(x) return(sites %>% filter(id %in% x))) %>% map(bind_rows, .id = "name") %>% 
        map(get_hypervolume, dims = c("MAP", "MAT", "bio4", "bio15"))
      
  # k: Absent clim by geo by sp      : same as h but run on sites where focal species absent only. 
      # Isolate included absences
      PA.ALL <- PA %>% reduce(multimerge, by = 0, all = T)
      PA.ALL <- PA.ALL[rownames(niche.geog.sp),colnames(niche.geog.sp)]
      PA.ALL[is.na(PA.ALL)] <- 0
      background <-  niche.geog.sp - PA.ALL # sites in the geo range where the species is absent.
      s.abs <- apply(background, 1, function(x) names(x[x==1]))
      # Hypervolume
      hvm.clim.geog.sp.abs <- purrr::map(s.abs, function(x) return(sites %>% filter(id %in% x))) %>% bind_rows(.id = "name") %>% 
        get_hypervolume(dims = c("MAP", "MAT", "bio4", "bio15"))
      
  # l: Absent clim by geo by sp, tbn : same as i but run on sites where focal species absent only. 
      # Isolate included absences
      
      PA2 <- map2(PA, niche.geog.sp.tbn, function(x,y) x[rownames(y),colnames(y)])
      background.geog.sp.tbn <- map2(niche.geog.sp.tbn, PA2, `-`)
      s.abs <- map(background.geog.sp.tbn, apply, 1, function(x) names(x[x==1]))
      # Hypervolume
      hvm.clim.geog.sp.tbn.abs <- purrr::map(s.abs, map, function(x) return(sites %>% filter(id %in% x))) %>% map(bind_rows, .id = "name") %>% 
        map(get_hypervolume, dims = c("MAP", "MAT", "bio4", "bio15"))
      
# CLIMATIC HYPERVOLUMES BY SP BY TIMEBIN GEOG INCLUSION, POOLED ####    
      niche.geog.sp.tbn.pooled <- reduce(niche.geog.sp.tbn, multimerge, by = 0, all = T) 
      niche.geog.sp.tbn.pooled[is.na(niche.geog.sp.tbn.pooled)] <- 0 
  # m: Uses Presences in by-timebin by-species geographic hvms to select sites
      s <- apply(niche.geog.sp.tbn.pooled, 1, function(x) names(x[x==1]))
      hvm.clim.geog.sp.tbn.pooled <- map(s, function(x) return(sites %>% filter(id %in% x))) %>% bind_rows(.id = "name") %>% 
        get_hypervolume(dims = c("MAP", "MAT", "bio4", "bio15"))
      
  # n: Absence
      s <- apply(niche.geog.sp.tbn.pooled, 1, function(x) names(x[x==0]))
      hvm.clim.geog.sp.tbn.abs.pooled <- map(s, function(x) return(sites %>% filter(id %in% x))) %>% bind_rows(.id = "name") %>% 
        get_hypervolume(dims = c("MAP", "MAT", "bio4", "bio15"))
      
      
# CALCULATE VOLUMES AND COLLATE DATA ####
hyper.list <- grep(pattern = "hvm.", ls(), value = TRUE) %>% map(get) 
layer <- hyper.list %>% map(map_chr, class) %>% sapply(table) %>% 
  names %>% `==`("Hypervolume")

hyper.list <- setNames(hyper.list, grep(pattern = "hvm.", ls(), value = TRUE))
layer2 <- which(map(hyper.list[layer] %>% map(map_dbl, ~.@Volume), length) > 3) %>% names
layer3 <- which(map(hyper.list[layer] %>% map(map_dbl, ~.@Volume), length) == 3) %>% names

v0 <- hyper.list[layer3] %>% map(map_dbl, ~.@Volume) %>% 
  map(cbind) %>% map(data.frame) %>% 
  reduce(multimerge, by = 0, all = T) %>% 
  setNames(layer3) %>% rownames_to_column("tbn") 

v1 <- hyper.list[layer2] %>% map(map_dbl, ~.@Volume) %>% 
  map(cbind) %>% map(data.frame) %>% 
  reduce(multimerge, by = 0, all = T) %>% 
  setNames(layer2) %>% rownames_to_column("name") 
  
v2 <- hyper.list[!layer] %>% map(map, map_dbl, ~.@Volume) %>% 
  map(map, cbind) %>% map(map, data.frame) %>%
   map(reduce, multimerge, by = 0, all = T) %>%  
   map2(hyper.list[!layer] %>% map(names), function(x,y) setNames(x, y)) %>% 
   map(rownames_to_column, "name") %>% map(melt) %>% 
   bind_rows(.id = "hvm") %>% 
   dcast(name+variable~hvm, value.var = "value")

v <- merge(v2, v1, by = "name", all = T) %>% merge(v0, by.x = "variable", by.y = "tbn")

# extinct only
ve <- v %>% filter(name %in% extinct) %>% filter(variable == "PLEI") %>% map_df(replace_na, 0) %>%  mutate(status = "victim")
# survivors only
vs <- v %>% filter(name %in% survivors) %>% map_df(replace_na, 0) %>% mutate(status = "survivor")

v <- rbind(vs, ve)

### NOTES ON FOREGROUND/BACKGROUND SELECTION ####   
    # combination of fore and background hvms depends on the biological question. 
# For instance: 
# c:a How much of the geographic space available in each time interval has the species accessed in each time interval?
# c:b How much of the total geographic space accessible to a species has the species accessed in each time interval?
# b:a (not too logical because b has nsp and a has 3 (n tbns), but could factorialize) 
    #How much of the geographic space in each time interval could the species potenially access? (could be more than 100%)
# f:d How much of the climatic space available in each time interval has the species accessed in each time interval?
# f:e How much of the total climatic space accessible to a species has the species accessed in each time interal?
# e:d (not too logical because e has nsp and d has 3 (n tbns), but could factorialize) 
    #How much of the climatic space in each time interval could the species potenially access? (could be more than 100%)


# i:f How much of the climatic space accessed by a species in each timebin was also part of the species' geographic range in that interval?
# h:e How much of the total climatic space accessible to a species was also in the species' total geographic range?
# i:e How much of the total climatic space accessed by a species was also in its geographic range in each timebin?
# i:h How much of the total climatically AND geographically accessible space of a species was accessed by the species in each timebin?
# i:k How does the total background (absent) climatic and geographic space available to each species compare with the climatic and geographic space it accessed in each timebin?

# etc., etc. so on and so forth. 
    
    
      