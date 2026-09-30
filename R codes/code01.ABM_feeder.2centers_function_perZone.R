# << Incorporating feeder attraction to ABM >> (crs:32618) # changed 11:38pm
# modification -- 'get ready for the ABM' section
#                             1. feederweightworld
#                             2.draw.ta.sl() function!

my_packages <- c("NetLogoR", "tidyr", "dplyr", "lubridate", "sf", "raster", "circular", 
                 "tidyverse", "sp", "CircStats", "REdaS")
lapply(my_packages, library, character.only = TRUE)
setwd("Z:/sl5145/ABM")

landras <- raster("Connetquot_nlcd_01.tif")
feederloc <- st_read("4posters_Connetquot.shp") %>% as_Spatial()
camloc <- st_read("Cams_Connetquot.shp") %>% as_Spatial()
inbuffer <- st_read("Connetquot_boundary-200.shp", crs = 32618)

# create 1) a world layer that contains the location of feeder
#        2) distance from the nearest feeders
#        3) landscape for deer initializing 
{
  emptyr <- landras
  values(emptyr) <- NA
  
  # 1) a world layer that contains the location of feeder
  feedras <- emptyr
  feedras[which(landras[] == 1 | landras[] == 0)] <- 0
  feedras[cellFromXY(feedras, feederloc)] <- 1
  #plot(feedras)
  
  # 2) distance from the nearest feeders
  feederloc_vect <- terra::vect(feederloc) # Convert the SpatialPointsDataFrame to an SpatVector (for terra distance function)
  distras <- terra::distance(terra::rast(emptyr), feederloc_vect) %>% raster()
  distras[which(is.na(landras[]))] <- NA
  #plot(distras)
  
  # 2-1) Zone map
  zoneras <- emptyr
  distvec <- seq(250, by = 250, length.out = 7)
  zoneras[] <- ifelse(distras[] <= distvec[1], 1, ifelse(distras[] <= distvec[2], 2, ifelse(distras[] <= distvec[3], 3, ifelse(distras[] <= distvec[4], 4, ifelse(distras[] <= distvec[5], 5, ifelse(distras[] <= distvec[6], 6, 7))))))
  plot(zoneras)
  points(feederloc); points(camloc, pch = 18)
  
  # 3) landscape for deer initializing (smaller)
  initras <- terra::rasterize(terra::vect(inbuffer), terra::rast(emptyr), 1) %>% raster()
  initras[which(landras[] == 0)] <- 0
  #plot(initras)
  
  # 4) Multiple feeder?
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
  distworld <- raster2world(distras)
  zoneworld <- raster2world(zoneras)
  initworld <- raster2world(initras)}



# [1] Create turtle groups randomly (The location is the HR center!) ----------------------------------------------
{
  N <- 176 # the number of groups !
  greenpatch <- patches(initworld) %>% as_tibble() %>% mutate(lu = of(initworld, patches(initworld))) %>% filter(lu == 1) %>% dplyr::select(-lu)
  forsprout <- greenpatch %>% sample_n(N) %>% as.matrix()
  t1 <- sprout(patches = forsprout, n = N) # randomly
  t1 <- turtlesOwn(t1, tVar = "sex", tVal = sample(c(rep("F", 56), rep("M", 120)), size  = N)) # female:male = 7:3 (210 F/ 90 M)
  
  plot(initworld); points(t1, pch = 23, col = t1$sex %>% recode("F" = "red", "M" = "blue"), bg = "white", lw = 2)
}

# [2] Each group is at least 200m away from each other ---------------------------------------------------------
# inRadius result: 'who' -- 'who' numbers from the second agent
#                  'id' -- agents from the first agent
{
  tryi <- 1
  while(!nrow(inRadius(t1, 6, t1) %>% filter(who+1 != id)) == 0){ 
    move.who <- inRadius(t1, 6, t1) %>% filter(who+1 != id) %>% pull(who) %>% unique()
    print(paste0(tryi, " try; n = ", length(move.who), " turtles have to be moved"))
    
    empty <- inRadius(greenpatch %>% as.matrix(), 6, t1, world = initworld) %>% pull(id)
    emptyid <- setdiff(1:nrow(greenpatch), empty)
    for (i in 1:length(move.who)) {
      t1[(move.who[i] + 1), ] <- moveTo(turtle(t1, who = move.who[i]), as.matrix(greenpatch[sample(emptyid, 1),]))
    }
    tryi <- tryi + 1
  }
  #plot(myworld); points(t1, pch = 23, col = t1$sex %>% recode("F" = "red", "M" = "blue"), bg = "white", lw = 2)
  
  # group assign
  t1 <- turtlesOwn(t1, tVar = "group", tVal = 1:N)
}

# [3] female group size = 5 --------------------------------------------------------------------------
{
  t1.sf <- turtles2sf(t1)
  ft1.sf <- t1.sf %>% filter(sex == "F")
  mt1.sf <- t1.sf %>% filter(sex == "M")
  
  whovec <- seq(200, 199 + (4*56))
  f.add <- ft1.sf[rep(seq_len(nrow(ft1.sf)), each = 5), ]
  f.add[-c(1, seq(6, nrow(f.add), 5)),]$who <- whovec # change the 'who' of added females
  shiftx <- runif(nrow(f.add), min = -1, max = 1)
  shifty <- runif(nrow(f.add), min = -1, max = 1)
  
  newx <- st_coordinates(f.add)[,"X"] + shiftx
  newy <- st_coordinates(f.add)[,"Y"] + shifty
  
  f.add <- f.add %>% mutate(newx = newx, newy = newy) %>% st_drop_geometry() %>% st_as_sf(coords = c("newx", "newy"))
  
  t2 <- sf2turtles(rbind(f.add,mt1.sf) %>% arrange(group))
  t2@.Data
  
  # the final deer turtles; group members start at the same place, 200m away from every group (including males)
  t3 <- createTurtles(n = nrow(t2), coords = t2@.Data[,1:2])
  t3 <- turtlesOwn(t3, tVar = "sex", tVal = t2@.Data[,"sex"])
  t3 <- turtlesOwn(t3, tVar = "group", tVal = t2@.Data[,"group"])
}

# [4] Assign leader for each group (random)  ---------------------------------------------------------------------------
{
  leadervec <- t3@.Data %>% as.data.frame() %>% group_by(group) %>% mutate(leader = first(who)) %>% pull(leader)
  t3 <- turtlesOwn(t3, tVar = "leader", tVal = leadervec)
  leaderYNvec <- t3@.Data %>% as.data.frame() %>% mutate(leaderYN = ifelse(who == leader, 1, 0)) %>% pull(leaderYN)
  t3 <- turtlesOwn(t3, tVar = "leaderYN", tVal = leaderYNvec)
}


# Get ready for the ABM ----------------------------------------------------------------------------------------------
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
    
    # Specify each deer's zonal information
    zone.pergroup <- of(world = zoneworld, agents = patchHere(zoneworld, deer.t.leader))
    groupzonedf <- data.frame(group = 1:N, zone = zone.pergroup)
    
    # join nearest feeder ID and zone information to turtles (generate t4)
    deer.t.sf <- deer.t.sf %>% left_join(groupfeederdf, by = "group") %>% left_join(groupzonedf, by = "group") 
    
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

plot(myworld, main = "timestep = 0 (initial)"); points(deer.t, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = "white", lw = 2)

t.end <- 12 * 30 # # of simulation days (2-hr each step), need burn-in period
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
      
      return(c(ta = ta.val[1], sl = sl.val[1], joint = ta.val[2]*sl.val[2]))
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
    
    i <- 1 
    fp <- select.feeder.prob.vec
    N.zone <- 7
    
    for(i in 1:t.end){
      for (ii in 1:N.zone) { # zonal loop start ()
        # MOVE GROUP LEADERS =======================================================================================
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
          angle.to.fd <- towards(zone.deer.t.leader.fd, as.matrix(feeder.center[which.fd, 1:2]))
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
        
        plot(myworld, main = paste0("timestep = ", i)); points(deer.t.leader, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = "white", lw = 2)
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
      
      routedf = data.frame()
      for (j in 1:nrow(deer.t)) {
        dist.taken <- sldf.both[j, "sl"] # a vector of SLs taken at this time step
        
        tt <- deer.t@.Data[j,] 
        prevpatches <- tt[c("prevX", "prevY")] %>% as.matrix() %>% t() # extract previous patch coords
        pass <- patchDistDir(myworld, prevpatches, dist = seq(0, dist.taken, 0.1), angle = deer.t@.Data[j,"heading"]) %>% unique()
        pass # identify all patches that animal has passed
        
        routedf <- rbind(routedf, pass)
      }
      
      routedf <- routedf %>% na.omit()
      passage.loc <- routedf %>% group_by(pxcor, pycor) %>% unique() %>% as.data.frame()
      passage.n <- routedf %>% group_by(pxcor, pycor) %>% arrange(pxcor) %>% summarise(n= n()) %>% pull(n)
      
      
      passageworld[passage.loc[,1], passage.loc[,2]] <-  passageworld[passage.loc[,1], passage.loc[,2]] + passage.n
    }
    
    
    colfunc <- colorRampPalette(c("white", "darkgreen", "red"))
    plot(passageworld, col = colfunc(30), main = "Number of passages")
    points(st_coordinates(feederloc.sf)[,"X"], st_coordinates(feederloc.sf)[,"Y"], pch = 8, lwd = 2, col = "yellow")
    points(hr.center$xcor, hr.center$ycor, pch = 8, lwd = 2, col = "blue")
    # initial plot
    #plot(myworld, main = "timestep = 0 (initial)"); points(deer.t, pch = 23, col = t3$sex %>% recode("1" = "red", "2" = "blue"), bg = "white", lw = 2)
    
    reslist[[simN]] <- sapply(1:7, function(x) mean(passageworld[which(zoneworld[] == x)]))
    
    simN <- simN + 1
  }
  
  return(reslist)
}

# Same deer initial population distribution
res0 <- runABM(t.end, select.feeder.prob = 0, totalNsimul = 5) # control
res0.1 <- runABM(t.end, select.feeder.prob = 0.1, totalNsimul = 5) # feeder 1
res0.3 <- runABM(t.end, select.feeder.prob = 0.3, totalNsimul = 5) # feeder 2

save(list = c("res0", "res0.1", "res0.3"), file = "Result/Twocenters_result_0_0.1_0.3.Rdata")


# Per zone plot
load("Result/Twocenters_result_0_0.1_0.3.Rdata")
map_dbl(res0, 1)
res0.mean <- sapply(1:6, function(x) mean(map_dbl(res0, x)))
res0.1.mean <- sapply(1:6, function(x) mean(map_dbl(res0.1, x)))
res0.3.mean <- sapply(1:6, function(x) mean(map_dbl(res0.3, x)))

dfdf <- data.frame(zone = 1:6, res0 = res0.mean, res0.1 = res0.1.mean, res0.3 = res0.3.mean)
dfdfvis <- dfdf %>% pivot_longer(cols = 2:4, names_to = "Type")
ggplot(dfdfvis, mapping = aes(zone, value, color = Type)) + geom_point() + geom_line() 

# do a buffer around the cam location?
# if they move small dist, then it means less passage rate in total than when it took a longer step -- how to tackle this?









