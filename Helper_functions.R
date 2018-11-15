#__________ Helper functions ______________________
#__________ by Anikó B. Tóth ______________________


require(tidyverse)
require(sp)

############### DATA MANIPULATION #########
# replaces the rownames of the data frame with the first column of that df.
namerows <- function(table){
  rownames(table) <- table[,1]
  table <- table[,2:ncol(table)]
  return(table)
}

# changes abundance data to presence-absence
tobinary <- function(PA.LIST){
  binary <- lapply(PA.LIST, function(x) {
    x <- x/x
    x[is.na(x)] <- 0
    return(x)
  })
  return(binary)
}

# facilitates merging multiple data frames using reduce() and by = 0 by restoring rownames after merge.
multimerge <- function(x, y, ...){
  return(merge(x, y, ...) %>% namerows())
}

# produces multiple subsamples of one presence-absence table, 
# specify number of samples and number of reps, removes empty rows

resamp <- function(PA, sites = 50, reps = 100, replace = F){
  samp <- list()
  for(j in 1:reps){
    samp[[j]] <- PA[, sample(1:ncol(PA), sites, replace = replace)]
    samp[[j]] <- samp[[j]][which(rowSums(samp[[j]]) > 0),]
  }
  return(samp)
}

#dist object to edge list
dist2edgelist <- function(z, sppDat, cname = ""){  
  k3 = as.matrix(z)
  dimnames(k3) <- list(rownames(sppDat), rownames(sppDat)) 
  
  xy <- t(combn(colnames(k3), 2))
  k3 <- data.frame(xy, dist=k3[xy], stringsAsFactors = F)

  k3 <- data.frame(k3, id = paste(k3$X1, k3$X2, sep = "-"), stringsAsFactors = F)
  colnames(k3) <- c("Sp1", "Sp2", "Z.Score", "id")
  
  return(k3)
}

# nested list with 2 layers to data frame
flat2 <- function(list2, layer1, layer2){
  map(list2, bind_rows, .id = layer2) %>% # consolidate subsample results
  bind_rows(.id = layer1) %>% return()
}


########### ANALYSIS FUNCTIONS ###########
# Convex hull calculation based on climate or biogeography 
find_hull.clim <- function(df) df[chull(df$MAP, df$MAT), ]
find_hull.geog <- function(df) df[chull(df$LATDD, df$LONGDD), ]

# Calculates the convex hull areas of each species occupancies in geographic ("geog") or climatic ("climate") space
niche_areas <- function(sitebyspecies, type = "geog"){
  if(!type %in% c("geog","climate")) stop("invalid type")
  if(type == "geog"){
    hulls <- lapply(unique(sitebyspecies$name), function(x) return(find_hull.geog(sitebyspecies[sitebyspecies$name==x,])))
    hulls <- lapply(hulls, function(x) return(cbind(x$LONGDD, x$LATDD)))
    hulls <- lapply(hulls, SpatialPoints, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
    hulls <- lapply(hulls, spTransform, CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
    chull.poly <- lapply(hulls, Polygon, hole=F) #will throw warnings for polygons with 3 or fewer vertices
    chull.area <- sapply(chull.poly, function(x) x@area/1e6)} #convert to km2
  if(type == "climate"){
    hulls <- lapply(unique(sitebyspecies$name), function(x) return(find_hull.clim(sitebyspecies[sitebyspecies$name==x,])))
    hulls <- lapply(hulls, function(x) return(cbind(x$MAP, x$MAT)))
    chull.poly <- lapply(hulls, Polygon, hole=F) #will throw warnings for polygons with 3 or fewer vertices
    chull.area <- sapply(chull.poly, function(x) x@area)
  }
  
  chull.area <- data.frame(unique(sitebyspecies$name), chull.area)
  
  chull.area
}

# Produces estimated combined niche Ni for each species, returned as a presence-absence table
# 1's indicate that species is present or that site parameters match species dispersal ability and environmental preferences.
get_niche_table <- function(sitebyspecies){
  hulls.clim <- lapply(unique(sitebyspecies$name), function(x) return(find_hull.clim(sitebyspecies[sitebyspecies$name==x,])))
  hulls.geog <- lapply(unique(sitebyspecies$name), function(x) return(find_hull.geog(sitebyspecies[sitebyspecies$name==x,])))
  
  sites <- unique(sitebyspecies[,c('id', "MAP", "MAT", "LATDD", "LONGDD")])
  require(sp)
  niche.clim <- do.call(rbind, lapply(hulls.clim, function(x) point.in.polygon(sites$MAT, sites$MAP, x$MAT, x$MAP, mode.checked=FALSE)))
  niche.geog <- do.call(rbind, lapply(hulls.geog, function(x) point.in.polygon(sites$LATDD, sites$LONGDD, x$LATDD, x$LONGDD, mode.checked=FALSE)))
  niche.clim <- tobinary(list(niche.clim))[[1]]
  niche.geog <- tobinary(list(niche.geog))[[1]]
  niche <- niche.clim + niche.geog
  colnames(niche) <- sites$id
  rownames(niche) <- unique(sitebyspecies$name)
  
  niche[niche==1] <- 0
  niche[niche==2] <- 1
  niche
}


# biotic/abiotic analysis, one iteration.
abio_single_niche <- function(y, niche, ss){
  out <- list()
  for(l in 1:length(y)){
    x <- y[[l]]
    id <- character()
    SFS <- numeric()        # full association calculated with all sites in ss1
    SFS.ss2 <- numeric()    # full association calculated with ss2 sites
    SFSa <- numeric()       # biotic association calculated with all Mij sites
    SFSa.ss2 <- numeric()   # biotic association inside Mij calculated with ss2

    c <- t(combn(colnames(x), 2)) # all pairwise combinations of species
    a <- numeric()
    b <- numeric()
    num <- numeric()
    for(i in 1:nrow(c)){  # for all pairwise combinations of species
      
      id[i] <- paste(c[i,1], c[i,2], sep = "-") # pair id
      n <- rownames(x)[which(rownames(x) %in% and.range.niche(c[i,1], c[i,2], niche))] # derive overlapping sites in species niches (Mij)
      a[i] <- sum(x[n,c[i,1]]) # number of occurrences of species a in Mij
      b[i] <- sum(x[n,c[i,2]]) # number of occurrences of species b in Mij
      num[i] <- length(n) # number of sites in Mij
      
      SFS[i] <- simpairs(t(x[,c(c[i,1],c[i,2])])) # full association with ss1 sites
      SFS.ss2[i] <- simpairs(t(x[sample(1:nrow(x), size = ss, replace = F),c(c[i,1],c[i,2])])) # full association with ss2 sites sampled from ss1
      
      if(length(n) >= ss) {
        
        SFSa.ss2[i] <- simpairs(t(x[n[sample(1:length(n), size = ss, replace = F)],c(c[i,1],c[i,2])])) # biotic association with ss2 sites sampled from ss1
        SFSa[i] <- simpairs(t(x[n,c(c[i,1],c[i,2])])) # biotic association with all Mij
      }
      
      if(length(n) < ss) { 
        SFSa.ss2[i] <- NA
        if(length(n) < 2) {SFSa[i] <- 0
          } else { SFSa[i] <- simpairs(t(x[n,c(c[i,1],c[i,2])]))} #SFS overlapping niche
        
      }
    }
    
    tb <- data.frame(timebin = names(y)[l], id = id, SFS, SFS.ss2, SFSa, SFSa.ss2, a, b, num, 
                     abtc.ss2 = SFS.ss2 - SFSa.ss2)
    head(tb)
    out[[l]] <- tb
  }
  out <- do.call(rbind, out)
  return(out)
}

# finds mutual range Mij of two species based on niche table
and.range.niche <- function(sp1, sp2, niche){
  
  if(sp1 %in% rownames(niche)){
    Sp1 <- names(which(niche[sp1,] == 1))
  }else(Sp1 <- NULL)
  
  
  if(sp2 %in% rownames(niche)){
    Sp2 <- names(which(niche[sp2,] == 1))
  }else(Sp2 <- NULL)
  
  return(Sp1[which(Sp1 %in% Sp2)])
}

## Co-occurence analysis functions ##
simpairs <- function(x){ # FETmP function
  samples = ncol(x)  #S
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  
  occs = array()
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #FETmP Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # B
      
      for (k in 0:a)
        z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      if(z[i,j]>=1) {z[i,j] <- 0.99999999999999994}
      z[i,j] = qnorm(z[i,j])
      z[j,i] = z[i,j]
    }
  }
  return(as.dist(z, diag = F, upper = F))
} 

simpairspd <- function(A, B, nsites){
  cscore <- rep(0, min(A, B)+1)
  for(a in max(0, A+B-nsites):min(A,B)){
    for (k in max(0, A+B-nsites):a){
      cscore[a+1] = 
      choose(B , k) * choose(nsites - B , A - k) / choose(nsites , A)}
      }
    return(cscore)
}  

##### Simple operations ####3
####
percpos <- function(x)
  length(which(x>0))/length(x)

percneg <- function(x)
  length(which(x<0))/length(x)

posnegzero <- function(x){
  out <- x > 0
  out[which(x==0)] <- "ZERO"
  out
}




######### PLOTTING ####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
