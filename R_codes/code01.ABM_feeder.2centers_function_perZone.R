# << Incorporating feeder attraction to ABM >> (crs:32618) # changed 11:38pm
# modification -- 'get ready for the ABM' section
#                             1. feederweightworld
#                             2.draw.ta.sl() function!

my_packages <- c("NetLogoR", "tidyr", "dplyr", "lubridate", "sf", "raster", "circular", 
                 "tidyverse", "sp", "CircStats", "REdaS", "terra")
lapply(my_packages, library, character.only = TRUE)

# choose one
setwd("Z:/sl5145/ABM")
setwd("D:/github/ABM_4-poster") # home computer


# import variables 
landras <- raster("Connetquot_nlcd_01.tif")
feederloc <- st_read("4posters_Connetquot.shp") %>% as_Spatial()
camloc <- st_read("Cams_Connetquot.shp") %>% as_Spatial()
cambuffer <- st_read("Cams_Connetquot.shp") %>% st_buffer(dist = 100) %>% as_Spatial()
inbuffer <- st_read("Connetquot_boundary-200.shp", crs = 32618)

# create 1) a world layer that contains the location of feeder
#        2) distance from the nearest feeders
#        3) landscape for deer initialization 
{
  emptyr <- landras
  values(emptyr) <- NA
  
  # 1) feeder location 01 raster 
  feedras <- emptyr
  feedras[which(landras[] == 1 | landras[] == 0)] <- 0
  feedras[cellFromXY(feedras, feederloc)] <- 1
  plot(feedras)
  
  # 1-2) cam location 01 raster
  camras <- emptyr
  camras[which(landras[] == 1 | landras[] == 0)] <- 0
  camras[cellFromXY(camras, camloc)] <- 1
  plot(camras)
  
  # 1-3) cam buffer w/ unique id raster
  cambuffras <- emptyr
  cambuffras[which(landras[] == 1 | landras[] == 0)] <- 0
  original_na <- is.na(cambuffras[]) # Remember which cells were originally outside your valid landscape.
  
  cam_ids <- rasterize(
    x = vect(cambuffer),    # buffer polygons
    y = rast(cambuffras),    # existing raster template
    field = "camera",              # ID to assign: 1–13
    background = 0,               # cells outside buffers
    touches = F                # include any cell touched by a buffer
  )
  
  # Convert back to RasterLayer for compatibility with your existing code
  cambuffras <- raster::raster(cam_ids)
  cambuffras[original_na] <- NA # Restore the original NA mask
  plot(cambuffras)
  
  # 2) distance from the nearest feeders
  feederloc_vect <- terra::vect(feederloc) # Convert the SpatialPointsDataFrame to an SpatVector (for terra distance function)
  distras <- terra::distance(terra::rast(emptyr), feederloc_vect) %>% raster()
  distras[which(is.na(landras[]))] <- NA
  plot(distras)
  
  # 2-1) Zone map
  zoneras <- emptyr
  distvec <- seq(250, by = 250, length.out = 7)
  zoneras[] <- ifelse(distras[] <= distvec[1], 1, ifelse(distras[] <= distvec[2], 2, ifelse(distras[] <= distvec[3], 3, ifelse(distras[] <= distvec[4], 4, ifelse(distras[] <= distvec[5], 5, ifelse(distras[] <= distvec[6], 6, 7))))))
  plot(zoneras)
  points(feederloc); points(camloc, pch = 18)
  
  # 3) landscape for deer initialization (smaller)
  initras <- terra::rasterize(terra::vect(inbuffer), terra::rast(emptyr), 1) %>% raster()
  initras[which(landras[] == 0)] <- 0
  plot(initras)
  
  # 4) Multiple feeder overlap plot
  # Feeder location -- feeder ID 1 from the top & right
  feederloc.sf <- feederloc %>% st_as_sf()
  feeder.buffer <- st_buffer(feederloc.sf, dist = 300) # 300m
  Nfeederras <- terra::rasterize(terra::vect(feeder.buffer), terra::rast(emptyr), fun = sum) %>% raster()
}

# List of worlds!
{
  landras[which(landras[] == 0)] <- NA
  myworld <- raster2world(landras)
  feedworld <- raster2world(feedras)
  camworld <- raster2world(camras)
  cambuffworld <- raster2world(cambuffras)
  distworld <- raster2world(distras)
  zoneworld <- raster2world(zoneras)
  initworld <- raster2world(initras)
  }


# ---------------------------------------------------------------------------------------------- #
# [1] Create turtle groups randomly (The location is the HR center!) -------
# ---------------------------------------------------------------------------------------------- #
{
  N <- 176 # the number of deer groups !
  greenpatch <- patches(initworld) %>% as_tibble() %>% 
    mutate(lu = of(initworld, patches(initworld))) %>% # adds each cell’s value 
    filter(lu == 1) %>% 
    dplyr::select(-lu)
  
  forsprout <- greenpatch %>% sample_n(N) %>% as.matrix()
  t1 <- sprout(patches = forsprout, n = N) # randomly
  t1 <- turtlesOwn(t1, tVar = "sex", tVal = sample(c(rep("F", 56), rep("M", 120)), 
                                                   size  = N)) # female:male = 7:3 (210 F/ 90 M)
  
  plot(initworld); points(t1, pch = 23, col = t1$sex %>% recode("F" = "red", "M" = "blue"), bg = "white", lw = 2)
}

# ---------------------------------------------------------------------------------------------- #
# [2] Each group is at least 200m away from each other -------------
# ---------------------------------------------------------------------------------------------- #
{
  # Eligible starting locations: one row per patch, with x/y coordinates.
  coords <- as.matrix(greenpatch)
  # eligible patches
  available <- seq_len(nrow(coords))
  
  n_groups <- nrow(t1)
  selected <- integer(n_groups)
  
  # Required separation in world units, assuming 30 m per unit.
  min_gap <- 200 / 30
  
  for (k in seq_len(n_groups)) {
    # Stop if this random arrangement runs out of eligible locations.
    if (length(available) == 0L) {
      stop("Placement ran out of candidates; retry or reduce group density.")
    }
    
    # Randomly select one available patch for group k.
    chosen <- available[sample.int(length(available), 1L)]
    selected[k] <- chosen
    
    # Calculate x/y differences from the chosen patch.
    offsets <- sweep(
      coords[available, , drop = FALSE],
      2,
      coords[chosen, ],
      "-")
    
    # Keep only patches at least 200 m from this group center.
    # Previously excluded patches stay excluded.
    available <- available[
      rowSums(offsets^2) >= min_gap^2
    ]
  }
  
  # Move the existing turtles; do not create new ones.
  t1 <- moveTo(t1, coords[selected, , drop = FALSE])
  t1 <- turtlesOwn(t1, tVar = "group", tVal = seq_len(n_groups))
}

# ---------------------------------------------------------------------------------------------- #
# [3] female group size = 5 -------
# ---------------------------------------------------------------------------------------------- #
{
  # Read the existing turtles' coordinates and stored attributes.
  group_data <- as.data.frame(t1@.Data)
  
  # Use the readable sex labels to identify female groups.
  # Keep the stored sex values for copying into t3 later.
  group_data$is_female <- turtles2sf(t1)$sex == "F"
  
  # Repeat female rows five times and male rows once.
  # Sorting keeps all members of each group together.
  deer_init <- group_data %>%
    mutate(group_size = if_else(is_female, 5L, 1L)) %>%
    tidyr::uncount(weights = group_size) %>%
    arrange(group)
  
  # Identify the expanded female rows.
  female_rows <- which(deer_init$is_female)
  n_females <- length(female_rows)
  
  # Give every female an independent offset around her group center.
  # As before, offsets are ±1 world unit along each axis.
  # Male coordinates remain unchanged.
  deer_init$xcor[female_rows] <- deer_init$xcor[female_rows] +
    runif(n_females, -1, 1)
  
  deer_init$ycor[female_rows] <- deer_init$ycor[female_rows] +
    runif(n_females, -1, 1)
  
  # Create the final population with fresh turtle IDs.
  # Select coordinates by name rather than by column position.
  t3 <- createTurtles(
    n = nrow(deer_init),
    coords = as.matrix(deer_init[, c("xcor", "ycor")])
  )
  
  # Restore the attributes needed by leader assignment and the ABM.
  t3 <- turtlesOwn(t3, tVar = "sex", tVal = deer_init$sex)
  t3 <- turtlesOwn(t3, tVar = "group", tVal = deer_init$group)
}

# therefore, 

# ---------------------------------------------------------------------------------------------- #
# [4] Assign leader for each group (random)  ---------------
# ---------------------------------------------------------------------------------------------- #
{
  # add 'leaderwho': who is leader of that group?
  leaderwhovec <- t3@.Data %>% as.data.frame() %>% group_by(group) %>% mutate(leaderwho = first(who)) %>% pull(leaderwho)
  t3 <- turtlesOwn(t3, tVar = "leaderwho", tVal = leaderwhovec)
  # add 'leaderYN': is that turtle leader of that group?
  leaderYNvec <- t3@.Data %>% as.data.frame() %>% mutate(leaderYN = ifelse(who == leaderwho, 1, 0)) %>% pull(leaderYN)
  t3 <- turtlesOwn(t3, tVar = "leaderYN", tVal = leaderYNvec)
}

# ---------------------------------------------------------------------------------------------- #
# Get ready for the ABM --------
# ---------------------------------------------------------------------------------------------- #
{
  # FEEDER ATTRACTION -- two centers approach ###########
  {
    # Move based on the distributions of TA and SL; turtles do not move to landuse == 0 with infection procedure
    # number of next step candidates
    deer.t <- t3
    deer.t.sf <- turtles2sf(deer.t)
    deer.t.leader.sf <- deer.t.sf %>% filter(leaderYN == 1)
    deer.t.leader <- sf2turtles(deer.t.leader.sf)
    
    # Feeder location -- feeder ID 1 from the top & right
    feederloc.sf <- rasterToPoints(world2raster(feedworld), fun = function(x){x == 1}, spatial = T) %>% st_as_sf() %>% mutate(N = 1:n())
    nearestfeedervec <- st_nearest_feature(deer.t.leader.sf, feederloc.sf) # Identify the ids of feeder that is the nearest to each turtle's HR center (per group and expand)
    groupfeederdf <- data.frame(group = 1:N, nearfeederID = nearestfeedervec)
    
    # cam location
    camloc.sf <- rasterToPoints(world2raster(camworld), fun = function(x){x == 1}, spatial = T) %>% st_as_sf()
    
    # Specify each deer's zonal information
    zone.pergroup <- of(world = zoneworld, agents = patchHere(zoneworld, deer.t.leader))
    groupzonedf <- data.frame(group = 1:N, zone = zone.pergroup)
    
    # join nearest feeder ID and zone information to turtles (generate t4)
    deer.t.sf <- deer.t.sf %>% left_join(groupfeederdf, by = "group") %>% left_join(groupzonedf, by = "group") 
    
    # new turtle with nearest feeder & current zone information 
    t4 <- sf2turtles(deer.t.sf)
    deer.t.leader.sf <- deer.t.sf %>% filter(leaderYN == 1)
    deer.t.leader <- sf2turtles(deer.t.leader.sf)
    
    # for movement, save the info of feeder centers of each turtle
    feedercoordsdf <- data.frame(nearfeederID = 1:6, X = st_coordinates(feederloc.sf)[,1], 
                                 Y = st_coordinates(feederloc.sf)[,2])
    
    feeder.center <- deer.t.leader@.Data %>% as.data.frame() %>% left_join(feedercoordsdf, by = "nearfeederID") %>% 
      dplyr::select(X, Y, who, group, leaderYN, nearfeederID) # this goes into the movement process
    feeder.center.sf <- st_as_sf(feeder.center, coords = c("X", "Y"))
  }
  
  # Home range center save!
  hr.center <- t4@.Data %>% as.data.frame() %>% filter(leaderYN == 1) %>% as.data.frame()
  hr.center.sf <- st_as_sf(hr.center, coords = c("xcor", "ycor"))
}

t4 # all turtles with nearest feeder & current zone info
deer.t.leader # leader turtles with nearest feeder & current zone info
feeder.center.sf # only leaders -- feeder location sf
hr.center.sf # only leaders -- HR center location sf

plot(myworld, main = "timestep = 0 (initial)"); points(t4, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = "white", lw = 2)


# ---------------------------------------------------------------------------------------------- #
# [5] runABM function  ------
# ---------------------------------------------------------------------------------------------- #
# number of simulation days (2-hr each step), need burn-in period
t.end <- 12 * 30 # 30 days (1 month)
# 


totalNsimul <- 5
select.feeder.prob.vec <- c(0, 0.1, 0.2, 0.3, 0.2, 0.1, 0)  # control -- value is 0

runABM <- function(t.end, select.feeder.prob.vec, totalNsimul = 5){
  # MOVEMENT RANDOM SAMPLING
  { # [1] LEADERS 
    # TA -- wrapped cauchy distribution (focused: the most likely TA (peak prob) is angle.to.hr & 
    #                                             the degree of spread (scale; rho here)) increases with the increasing distance from hr!)
    # SL -- weibull distribution (Scaled: scale changes with increasing distance from their home range center)
    sample.ta <- function(k, distances){
      rho.0.x <- rbeta(1, 4.24, 5.96)
      rho.0 <- (rho.0.x + 1)/2
      rho.inf.x <- rbeta(1, 28.48, 32.63)
      rho.inf <- sapply(1, function(i) {
        sample(c(1, rho.inf.x), size = 1, prob = c(0.3, 0.7))
      })
      y <- (-2.79 * rho.inf) - 4.32
      MSE <- 1.21
      gamma.rho <- rnorm(1, y, MSE) %>% exp()
      # intend: as dist increases, rho.t increases                       #################### weird?
      rho.t <- rho.inf + (rho.0 - rho.inf)*exp(-gamma.rho*distances[k])
      
      val <- rwrpcauchy(1, rad.to.hr[k], rho.t) 
      prob <- dwrpcauchy(val, rad.to.hr[k], rho.t)
      
      #curve(dwrpcauchy(x, rad.to.hr[k], 0.5), from = 0, to = 6) 
      
      return(c(val = val, prob = prob))
    }
    
    sample.sl <- function(k, distances){
      alpha.t <- rlnorm(1, -0.077, 0.058)
      
      beta.0 <- rlnorm(1, 4.70, 0.23)
      gamma.beta <- rnorm(1, 0.17, 0.085)
      beta.t <- beta.0 + gamma.beta*distances[k]
      while(beta.t < 0){
        beta.0 <- rlnorm(1, 4.70, 0.23)
        gamma.beta <- rnorm(1, 0.17, 0.085)
        beta.t <- beta.0 + gamma.beta*distances[k]
      }
      
      val <- rweibull(1, alpha.t, beta.t)
      prob <- dweibull(val, alpha.t, beta.t)
      
      #curve(dweibull(x, alpha.t, beta.t), from = 0, to = 1000)
      while(val > 5421 | is.na(val)){
        val <- rweibull(1, alpha.t, beta.t)
        prob <- dweibull(val, alpha.t, beta.t)
      }
      return(c(val = val, prob = prob))
    }
    
    # Ultimate draw TA and SA combined for leaders
    # k: number of candidates
    # distances: distance from hr center
    draw.ta.sl <- function(k, distances){
      ta.val <- sample.ta(k, distances)
      ta.val[1] <- ta.val[1] %>% rad2deg()
      sl.val <- sample.sl(k, distances)
      sl.val[1] <- sl.val[1]/30
      test.t <- zone.deer.t.leader.hr[k] %>% right(ta.val[1]) %>% fd(sl.val[1]) # zonal
      
      while (is.na(of(world = myworld, agents = patchHere(myworld, test.t)))) {
        ta.val <- sample.ta(k, distances)
        ta.val[1] <- ta.val[1] %>% rad2deg()
        sl.val <- sample.sl(k, distances)
        sl.val[1] <- sl.val[1]/30
        test.t <- zone.deer.t.leader.hr[k] %>% right(ta.val[1]) %>% fd(sl.val[1]) # zonal
      }
      
      return(c(ta = ta.val[1], 
               sl = sl.val[1], 
               joint = ta.val[2]*sl.val[2]) # joint probability 
             )
    }
    
    # [2] FOLLOWERS
    sample.sl.follow <- function(k){
      sl.val <- rexp(1, rate = 0.00406)
      test.t <- deer.t.follower[k] %>% fd(sl.val/30)
      while (is.na(of(world = myworld, agents = patchHere(myworld, test.t))) | sl.val > 5421) {
        sl.val <- rexp(1, rate = 0.00406)
        test.t <- deer.t.follower[k] %>% fd(sl.val/30)
      }
      
      return(sl.val)
    }
  }
  
  # save results for N simulations
  reslist <- vector(mode = "list", length = totalNsimul)
  
  simN <- 1
  # Run the model for totalNsimul times!
  while(simN < totalNsimul + 1){
    # Initialization (turtles and passage world)
    deer.t <- t4
    deer.t.sf <- turtles2sf(deer.t)
    
    deer.t.leader.sf <- deer.t.sf %>% filter(leaderYN == 1)
    deer.t.follower.sf <- deer.t.sf %>% filter(leaderYN == 0)
    deer.t.leader <- sf2turtles(deer.t.leader.sf)
    deer.t.follower <- sf2turtles(deer.t.follower.sf)
    
    passageworld <- createWorld(minPxcor(myworld), maxPxcor(myworld), minPycor(myworld), maxPycor(myworld), data = 0)
    
    i <- 1 # timestep
    fp <- select.feeder.prob.vec # feeder selection probability 
    N.zone <- 7
    
    for(i in 1:t.end){ # t.end (for each time step)
      # MOVE GROUP LEADERS =======================================================================================
      
      for (ii in 1:N.zone) { # zonal loop.. for each zone,
        zone.deer.t.leader.sf <- deer.t.leader.sf %>% filter(zone == ii)
        
        turtleN.zone <- nrow(zone.deer.t.leader.sf)
        v <- runif(turtleN.zone, 0, 1)
        
        N.hr <- length(which(v >= fp[ii])) # leaders heading HR center
        which.hr <- which(v >= fp[ii])
        N.fd <- length(which(v < fp[ii])) # leaders heading to feeders
        which.fd <- which(v < fp[ii])
        
        if(N.hr > 0){ # HR center use
          zone.deer.t.leader.sf.hr <- zone.deer.t.leader.sf[which.hr,]
          zone.deer.t.leader.hr <- sf2turtles(zone.deer.t.leader.sf.hr)
          
          # calculate distance from HR center and the current location (leaders)
          hr.center.sf.zone <- hr.center.sf %>% filter(zone == ii) %>% slice(which.hr)
          distances <- st_distance(zone.deer.t.leader.sf.hr, hr.center.sf.zone, by_element = T)*30
          
          # calculate the turning angle that make it face HR center (leaders)
          angle.to.hr <- towards(zone.deer.t.leader.hr, cbind(st_coordinates(hr.center.sf.zone)[,1], st_coordinates(hr.center.sf.zone)[,2]))
          turning.to.hr <- sapply(1:N.hr, function(x){
            currentH <- zone.deer.t.leader.hr@.Data[x,"heading"]
            
            if(currentH > angle.to.hr[x]){
              angle.to.hr[x] + (365 - currentH)
            }else if(currentH < angle.to.hr[x]){
              angle.to.hr[x] - currentH
            }else if(currentH == angle.to.hr[x]){
              0
            }
          })
          
          rad.to.hr <- deg2rad(turning.to.hr) # required for TA sampling
          
          # TA and SL selection -- 1:N
          cand.n <- 10 # Number of candidate steps!
          movemat.hr <- sapply(1:N.hr, function(x){
            candsmat <- replicate(cand.n, draw.ta.sl(x, distances))
            candsmat[, which.max(candsmat[3,])]
          })
          
          # turn and go straight
          zone.deer.t.leader.hr.upd <- zone.deer.t.leader.hr %>% right(movemat.hr[1,]) %>% fd(movemat.hr[2,])
         
        }
        
        if(N.fd > 0){ # Feeder center use
          zone.deer.t.leader.sf.fd <- zone.deer.t.leader.sf[which.fd,]
          zone.deer.t.leader.fd <- sf2turtles(zone.deer.t.leader.sf.fd)
          
          # TA -- It directly faces the feeder
          feeder.center.zone <- feeder.center %>% filter(group %in% zone.deer.t.leader.sf.fd$group)
          angle.to.fd <- towards(zone.deer.t.leader.fd, as.matrix(feeder.center.zone[, 1:2]))
          turning.to.fd <- sapply(1:N.fd, function(x){
            currentH <- zone.deer.t.leader.fd@.Data[x,"heading"]
            
            if(currentH > angle.to.fd[x]){
              angle.to.fd[x] + (365 - currentH)
            }else if(currentH < angle.to.fd[x]){
              angle.to.fd[x] - currentH
            }else if(currentH == angle.to.fd[x]){
              0
            }
          })
          
          # SL selection
          # calculate distance from HR center and the current location (leaders)
          distances <- st_distance(zone.deer.t.leader.sf.fd, feeder.center.sf[which.fd, ], by_element = T)*30
          SLs <- sapply(1:N.fd, function(x){
            sl.val <- sample.sl(x, distances)/30
            test.t <- zone.deer.t.leader.fd[x] %>% right(turning.to.fd[x]) %>% fd(sl.val[1])
            while(is.na(of(world = myworld, agents = patchHere(myworld, test.t)))){
              sl.val <- sample.sl(x, distances)/30
              test.t <- zone.deer.t.leader.fd[x] %>% right(turning.to.fd[x]) %>% fd(sl.val[1])
            }
            
            return(sl.val[1])
          })
          
          # Create movemat
          movemat.fd <- as.matrix(data.frame(TA = turning.to.fd, SL = SLs) %>% t())
          
          # turn and go straight
          zone.deer.t.leader.fd.upd <- zone.deer.t.leader.fd %>% right(movemat.fd[1,]) %>% fd(movemat.fd[2,])
        }
        
        # Merge deer that moved towards HR center and feeders + movemat
        if(N.fd > 0){
          zone.deer.t.leader.sf <- rbind(turtles2sf(zone.deer.t.leader.hr.upd), turtles2sf(zone.deer.t.leader.fd.upd)) %>% arrange(who)
          zone.deer.t.leader <- sf2turtles(zone.deer.t.leader.sf)
          
          movemat <- data.frame(who = zone.deer.t.leader.sf$who, dist = 0)
          movemat[which.hr, "dist"] <- movemat.hr[2,]
          movemat[which.fd, "dist"] <- movemat.fd[2,]
        }else if(N.fd == 0){
          zone.deer.t.leader <- zone.deer.t.leader.hr.upd
          
          movemat <- data.frame(who = zone.deer.t.leader.sf$who, dist = 0)
          movemat$dist <- movemat.hr[2,]
        }
        
        assign(paste0("zone", ii, ".deer.t.leader"), zone.deer.t.leader)
        assign(paste0("zone", ii, ".movemat"), movemat)
        
        # plot(myworld, main = paste0("timestep = ", i)); points(deer.t.leader, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = "white", lw = 2)
        
      } # the end of zonal loop
      
      # merge zonal loop results (turtles and movemat)
      namelist <- lapply(1:N.zone, function(x) get(paste0("zone", x, ".deer.t.leader")))
      deer.t.leader.sf <- do.call(rbind, namelist) %>% turtles2sf() %>% arrange(who)
      deer.t.leader <- sf2turtles(deer.t.leader.sf)
      
      namelist2 <- lapply(1:N.zone, function(x) get(paste0("zone", x, ".movemat")))
      movemat <- do.call(rbind, namelist2) %>% arrange(who)
      
      
      # MOVE GROUP FOLLOWERS =======================================================================================
      # TA (towards the leader)
      leaderloc <- patchHere(myworld, deer.t.leader) %>% as_tibble() %>% mutate(group = deer.t.leader@.Data[,"group"])
      leaderloc <- data.frame(group = deer.t.follower@.Data[,"group"]) %>% left_join(leaderloc, by = "group") %>% dplyr::select(pxcor, pycor) %>% as.matrix()
      
      deer.t.follower <- face(deer.t.follower, leaderloc)
      
      # SL (random sampling from the exponential distribution with rate)
      sl.follow <- sapply(1:nrow(deer.t.follower), sample.sl.follow)/30
      
      deer.t.follower <- deer.t.follower %>% fd(sl.follow)
      
      # MERGE LEADERS AND FOLLOWERS ================================================================================
      deer.t.sf <- rbind(turtles2sf(deer.t.leader), turtles2sf(deer.t.follower)) %>% arrange(who)
      deer.t <- deer.t.sf %>% sf2turtles()
      
      # update
      deer.t.leader.sf <- deer.t.sf %>% filter(leaderYN == 1)
      
      # plot
      #plot(myworld, main = paste0("timestep = ", i)); points(deer.t, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = "white", lw = 2)
      #plot(myworld, main = paste0("timestep = ", i)); points(deer.t, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = t3$leaderYN %>% recode("0" = "white", "1" = "black"), lw = 2)
  
      # UPDATE PASSAGE WORLD =======================================================================================
      sldf1 <- data.frame(who = deer.t.leader.sf$who, sl = movemat[,"dist"])
      sldf2 <- data.frame(who = deer.t.follower.sf$who, sl = sl.follow)
      sldf.both <- rbind(sldf1, sldf2) %>% arrange(who)
      
      # Read turtle data once for this timestep.
      deer_data <- deer.t@.Data
      
      # Allocate one list entry per deer.
      routes <- vector("list", nrow(deer_data))
      
      for (j in seq_along(routes)) {
        
        # Previous coordinates, preserved as a one-row matrix.
        previous_xy <- deer_data[
          j, c("prevX", "prevY"), drop = FALSE
        ]
        
        # Keep the same route sampling interval as your original code.
        passed <- patchDistDir(
          myworld,
          previous_xy,
          dist = seq(0, sldf.both[j, "sl"], by = 0.1),
          angle = deer_data[j, "heading"]
        )
        
        # Each deer contributes at most one count per patch per timestep.
        routes[[j]] <- as.data.frame(unique(passed))
      }
      
      # Combine once, then calculate coordinates and counts together.
      passage_counts <- dplyr::bind_rows(routes) %>%
        filter(!is.na(pxcor), !is.na(pycor)) %>%
        count(pxcor, pycor, name = "n")
      
      if (nrow(passage_counts) > 0L) {
        passageworld[
          passage_counts$pxcor, passage_counts$pycor
        ] <- passageworld[
          passage_counts$pxcor, passage_counts$pycor
        ] + passage_counts$n
      }
      
      # routedf = data.frame()
      # for (j in 1:nrow(deer.t)) {
      #   dist.taken <- sldf.both[j, "sl"] # a vector of SLs taken at this time step
      #   
      #   tt <- deer.t@.Data[j,] 
      #   prevpatches <- tt[c("prevX", "prevY")] %>% as.matrix() %>% t() # extract previous patch coords
      #   pass <- patchDistDir(myworld, prevpatches, dist = seq(0, dist.taken, 0.1), 
      #                        angle = deer.t@.Data[j,"heading"]) %>% unique()
      #   pass # identify all patches that animal has passed
      #   
      #   routedf <- rbind(routedf, pass)
      # }
      # 
      # routedf <- routedf %>% na.omit()
      # passage.loc <- routedf %>% group_by(pxcor, pycor) %>% unique() %>% as.data.frame()
      # passage.n <- routedf %>% group_by(pxcor, pycor) %>% arrange(pxcor) %>% 
      #   summarise(n=n(), .groups ="drop_last") %>% pull(n)
      # 
      # 
      # passageworld[passage.loc[,1], passage.loc[,2]] <-  passageworld[passage.loc[,1], passage.loc[,2]] + passage.n
    
      cat(paste0(simN, "-", i, " "))
    } # end of each time step
    
    # plotting ---- # 
    colfunc <- colorRampPalette(c("white", "darkgreen", "red"))
    
    plot(passageworld, col = colfunc(30), main = "Number of passages")
    points(st_coordinates(feederloc.sf)[,"X"], st_coordinates(feederloc.sf)[,"Y"], 
           pch = 8, lwd = 2, col = "yellow")
    points(hr.center$xcor, hr.center$ycor, pch = 8, lwd = 2, col = "blue")
    points(camloc.sf, pch = 1, lwd = 2)
    # initial plot
    plot(myworld, main = "timestep = 0 (initial)"); points(t4, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = "white", lw = 2)
    # ------------- # 
    
    reslist[[simN]] <- passageworld
    
    simN <- simN + 1
  }
  
  return(reslist)
}

# Same deer initial population distribution
res.control <- runABM(t.end, select.feeder.prob.vec = rep(0, 7), totalNsimul = 3) # control
res.default <- runABM(t.end, select.feeder.prob.vec = c(0, 0.1, 0.2, 0.3, 0.2, 0.1, 0) , totalNsimul = 3) # deafult vector testing
res.strong <- runABM(t.end, select.feeder.prob.vec = c(0, 0.3, 0.4, 0.5, 0.2, 0.1, 0) , totalNsimul = 3) # deafult vector testing


colfunc <- colorRampPalette(c("white", "darkgreen", "red"))
plot(res.control[[1]], col = colfunc(30))
points(camloc.sf, pch = 1, lwd = 2)

lapply(1:13, function(x) {
  mean(
    res.control[[1]][][which(cambuffworld[] == x)],
    na.rm = TRUE
  )
})



# <to do list>
# change the structure of fp to curve (less parameters)
# calibration of fp .. based on a buffer around the cam location?



# export the passageworld and convert it into a raster with original extent and crs
back2ras <- function(x){
  # Convert the NetLogoR world into a RasterLayer.
  x_ras <- NetLogoR::world2raster(x)
  
  # Confirm it has the same grid dimensions as the original landscape.
  stopifnot(
    nrow(x_ras) == nrow(x),
    ncol(x_ras) == ncol(x)
  )
  
  # Restore the original geographic extent and coordinate system.
  # Setting the extent also restores the original cell resolution.
  raster::extent(x_ras) <- raster::extent(landras)
  raster::crs(x_ras) <- raster::crs(landras)
  
  return(x_ras)
}


raster::writeRaster(
  back2ras(res.control[[2]]),
  filename = "Results.prelim/passageworld.tif",
  format = "GTiff",
  overwrite = TRUE
)
raster::writeRaster(
  back2ras(cambuffworld),
  filename = "Results.prelim/cambuffworld.tif",
  format = "GTiff",
  overwrite = TRUE
)




