my_packages <- c("tidyr", "dplyr", "ggplot2", "Distance", "mrds", "readxl", "lubridate")
lapply(my_packages, library, character.only = TRUE)

setwd("D:/0.Columbia/0. Data/2.CDC_Suffolk/_Passage rate 2023")

camnames <- paste0("C", 1:13)
DDsheet <- read.csv("DDsheet.csv") %>% mutate(camN = paste0(substr(Park, 1, 1), camera)) %>% filter(Park == "Connetquot")

# [1] Calculate EDD -----------------------------------------------------------------------------------------------
edd.df <- as.data.frame(matrix(nrow = length(camnames), ncol = 3)); colnames(edd.df) <- c("camN", "Pa", "EDD")
edd.df$camN <- camnames
for (i in 1:length(camnames)) {
  # CSV dataset import
  cam <- camnames[i]
  df <- read.csv(paste0(cam, ".csv")) %>% dplyr::select(Date, Time, n_images, Species, N, Passage) %>%
    filter(!(Species == ""))
  
  # select only deer and mark with/without feeders
  df2 <- df %>% 
    filter(Species == "Deer", Date >= "2023-10-01") %>% 
    filter(Passage != "12.5m~") %>% 
    mutate(Passage = factor(Passage, levels = c("0m~", "2.5m~", "5m~", "7.5m~", "10m~")))

  # Connetquot 4-poster removal date: 12/8
  df2 <- df2 %>% mutate(Date = as.POSIXct(Date)) %>% 
    mutate(feeder = ifelse(Date <= "2023/12/08", 1, 0),
           month = format(Date, "%m"),
           day = format(Date, "%d"))
  
  # TEMPORAL WINDOW FOR INDEPENDENT DETECTION ------------ # 
  df3 <- df2 %>% distinct(pick(Date, Time), .keep_all = T) %>%
    mutate(dhms = as_datetime(paste0(Date, Time)),
           dt = dhms - lag(dhms,1))
  df3[1, "dt"] <- 0
  
  # Temporal window value -----
  window.secs <- 60 * 10 # amount of seconds (Decision here!)
  
  # using the above time window, identify unique individuals or groups,
  # while keeping track of their passage zones!!! 
  df4 <- df3 %>% mutate(new_group = cumsum(dt > window.secs),
                        unique_passage = !duplicated(paste(new_group, Passage))) %>% 
    filter(unique_passage) %>%
    dplyr::select(-c("new_group", "unique_passage"))
  df4 <- rbind(df3[1,], df4) # add the first row from df3 bc it is removed due to time window thing
  # ------------------------------------------------------ # 
  
  # Duplicate rows if N > 1 and append
  Nlist <- setdiff(na.omit(as.numeric(unique(df4$N))), 1)
  if(any(Nlist)){
    add.df <- as.data.frame(matrix(ncol = 11))
    colnames(add.df) <- c("Date", "Time", "n_images", "Species", "N", "Passage", 
                          "feeder", "month", "day", "dhms", "dt")
    for (ii in 1:length(Nlist)) {
      sub.df <- df4 %>% filter(N == Nlist[ii])
      val <- Nlist[i]
      
      sub.df <- sub.df[rep(seq_len(nrow(sub.df)), each = val), ]
      add.df <-  rbind(add.df, sub.df)
      add.df <- add.df[-1,]
    }
    df5 <- df4 %>% filter(N == 1) %>% rbind(add.df) %>% arrange(Date, Time)
  }else{
    df5 <- df4
  }
  
  
  # Fit a probability distribution function using {mrds} package
  dfds <- df5 %>% mutate(distance = (as.numeric(Passage)-1)*2.5,
                                      observer = 1,
                                      detected = 1,
                                      object = 1:nrow(df5),
                                      distbegin = distance,
                                      distend = distance+2.5)
  
  barplot(summary(df4$Passage))
  
  # hn: half-normal function
  # hr: hazard-rate function
  hn.model <- ddf(data = dfds, 
                  dsmodel = ~cds(key = "hn"), 
                  meta.data = list(point = F, binned = T, breaks = 0:5*2.5))  ########
  hr.model <- ddf(data = dfds, 
                  dsmodel = ~cds(key = "hr"), 
                  meta.data = list(point = F, binned = T, breaks = 0:5*2.5))  ########
  # identify which model shows the lower AIC
  nn <- which.min(c(hn.model$criterion, hr.model$criterion))
  
  if(nn == 1){
    dpf <- hn.model
    Pa <- predict(dpf)$fitted[1]
  }else{
    dpf <- hr.model
    Pa <- predict(dpf)$fitted[1]
  }
  
  png(filename = paste0("Plots/fitted_", cam, ".png"),
    width = 5, height = 5, units = "in",
    res = 300)
  plot(dpf)
  title(main = if (nn == 1) {
    paste0(cam, "; Half normal function")
  } else {
    paste0(cam, "; Hazard rate key function")
  })
  dev.off()
  
  
  # Calculate passage rate by setting w (truncation distance)
  # w either from 1) distance frequency data or 2) field-measured max. detection distance
  DDval <- DDsheet %>% filter(camN == cam) %>% pull(detection_)
  
  w <- 12.5
  EDD <- w * sqrt(Pa)
  
  edd.df[i, "Pa"] <- Pa
  edd.df[i, "EDD"] <- EDD
}

edd.df 


# [2] Compute passage rate before and after the feeder ---------------------------------------------------------------
camnames <- paste0("C", 1:13)
camnames <- setdiff(camnames, "C10")

pvdf <- as.data.frame(matrix(nrow = length(camnames), ncol = 9)); colnames(pvdf) <- c("camN", "withfeeder", "nofeeder", "t.stats", "t.pval", "glmnb.coef", "glmnb.pval", "Count.nofeed", "Count.withfeed")
pvdf$camN <- camnames
for (i in 1:length(camnames)) {
  {
  cam <- camnames[i]
  
  df <- read.csv(paste0(cam, ".csv")) %>% dplyr::select(Date, Time, n_images, Species, N, Passage) %>% filter(Date >= "2023-10-01")
  initDate <- df[1,"Date"]
  df2 <- df %>% filter(Species == "Deer") %>%
    mutate(Passage = factor(Passage, levels = c("0m~", "2.5m~", "5m~", "7.5m~", "10m~", "12.5m~")))
  df2 <- df2 %>% mutate(Date = as.POSIXct(Date)) %>% mutate(feeder = ifelse(Date <= "2023/12/08", 1, 0),
                                                            month = format(Date, "%m"),
                                                            day = format(Date, "%d"))
  df3 <- df2
  
  # TIME WINDOW
  df3 <- df3 %>% distinct(pick(Date, Time), .keep_all = T) %>%
    mutate(dhms = as_datetime(paste0(Date, Time)),
           dt = dhms - lag(dhms,1))
  df3[1, "dt"] <- 0

  window.mins <- 60 * 30 # amount of seconds (or mins) (Decision here!)
  df4 <- df3 %>%  mutate(new_group = cumsum(dt > window.mins),
                         unique_passage = !duplicated(paste(new_group, Passage))) %>% filter(unique_passage) %>%
    dplyr::select(-c("new_group", "unique_passage"))
  df4 <- rbind(df3[1,], df4) # add the first row from df3 bc it is removed due to time window thing
  
  # COUNT ZEROS
  dates <- seq(as.Date(initDate), as.Date("2023-12-21"), by="day") ######### "2023-10-01" to "initDate"
  cam.dates <- df4$Date %>% unique() %>% as_date()
  nodeerdates <- setdiff(dates, cam.dates) %>% as_date()
  zeronumb.feeder1 <- sum(nodeerdates < as_date("2023-12-08"))
  zeronumb.feeder0 <- sum(nodeerdates > as_date("2023-12-08"))
  
  # Calculate # of passage/day before and after the feeder removal
  passage.feeder0 <- df4 %>% filter(feeder == 0)
  passage.feeder1 <- df4 %>% filter(feeder == 1)
  
  pv1 <- passage.feeder1 %>% group_by(month, day) %>% tally() %>% pull(n) %>% append(rep(0, zeronumb.feeder1))
  pv2 <- passage.feeder0 %>% group_by(month, day) %>% tally() %>% pull(n) %>% append(rep(0, zeronumb.feeder0))
  
  # EDD
  EDDval <- edd.df %>% filter(camN == cam) %>% pull(EDD)
  pv1.edd <- pv1/EDDval
  pv2.edd <- pv2/EDDval
  
  # t.test
  testres <- t.test(pv1.edd, pv2.edd)

  # Negative binomial
  testdf <- rbind(data.frame(PR = pv1, group = "yfeed"), data.frame(PR = pv2, group = "nfeed"))
  }
  model <- glm(PR ~ group, testdf, family = poisson(link = "log"))
  summary(model)
  performance::check_overdispersion(model)
  model.nb <- summary(MASS::glm.nb(PR ~ group, testdf))
  
  
  pvdf[i, "withfeeder"] <- mean(pv1.edd)
  pvdf[i, "nofeeder"] <- mean(pv2.edd)
  pvdf[i, "t.stats"] <- testres$statistic
  pvdf[i, "t.pval"] <- testres$p.value
  #pvdf[i, "glmnb.coef"] <- coef(summary(model))[2,1]
  #pvdf[i, "glmnb.pval"] <- coef(summary(model))[2,4]
  pvdf[i, "glmnb.coef"] <- model.nb$coefficients[2,1]
  pvdf[i, "glmnb.pval"] <- model.nb$coefficients[2,4]
  
  
  if(model.nb$coefficients[2,4] < 0.05){
    pvdf[i, "Count.nofeed"] <- exp(model.nb$coefficients[1,1])
    pvdf[i, "Count.withfeed"] <- exp(model.nb$coefficients[1,1]) + exp(model.nb$coefficients[2,1])
  }else{
    pvdf[i, "Count.nofeed"] <- exp(model.nb$coefficients[1,1])
    pvdf[i, "Count.withfeed"] <- exp(model.nb$coefficients[1,1])
  }

}

pvdf <- pvdf %>% mutate(t.sig = ifelse(t.pval < 0.05, "yes", ""),
                        glm.sig = ifelse(glmnb.pval < 0.05, "yes", ""))
pvdf <- pvdf %>% left_join(DDsheet %>% dplyr::select(camN, zone, NearDist), by = "camN") %>% arrange(zone)
pvdf

# [3] Analyze the results ---------------------------------------------------------------------------------------
#save(list = c("edd.df", "pvdf"), file = "Results.Rdata")
#load("Results.Rdata")

pvdf2 <- pvdf %>% left_join(edd.df, by = "camN")
pvdf2 <- pvdf2 %>% mutate(MultipleF = c(0, 1, 1, rep(0, 9)))

fit1 <- lm(withfeeder ~ NearDist, data = pvdf2) # no sig. relationship at all
fit0 <- lm(nofeeder ~ NearDist, data = pvdf2) # no sig. relationship at all
summary(fit1)
summary(fit0)

ggplot(pvdf2, mapping = aes(NearDist, withfeeder, label = camN)) +
  geom_point() + geom_smooth(method = lm) +
  labs(title = "Deer activity ~ Distance, with feeders") + ggrepel::geom_text_repel() + ylim(0, 0.7) +
  geom_vline(xintercept = seq(250, 1250, 250), linetype = "dotted") + 
  xlab("Distance from feeders") + ylab("Passage rate (# of passage/day/m)")
ggplot(pvdf2, mapping = aes(NearDist, nofeeder, label = camN)) + geom_point() + geom_smooth(method = lm) +
  labs(title = "Deer activity ~ Distance, without feeders") + ggrepel::geom_text_repel() + ylim(0, 0.7) + 
  geom_vline(xintercept = seq(250, 1250, 250), linetype = "dotted") +
  xlab("Distance from feeders") + ylab("Passage rate (# of passage/day/m)")

# 'feeder removal' as a covariate
pvdf3 <- pvdf2 %>% pivot_longer(cols = withfeeder:nofeeder, names_to = "PR")
summary(lm(value ~ NearDist + MultipleF + PR, data = pvdf3))  
  # this means with feeder, 50% increase in passage rate in entire area
  # without feeder, one feeder zone -- 0.11/ but multiple feeder zone, 0.25/ 



# Passage rate plot before/ after the feeder removal
pvdf.long <- pvdf2 %>% pivot_longer(cols = withfeeder:nofeeder, names_to = "Feeder?") %>% 
  mutate(glm.sig = if_else(glm.sig == "", "no", glm.sig)) %>%
  mutate(sig = factor(glm.sig, levels = c("yes", "no")),
         `Feeder?` = factor(`Feeder?`, levels = c("withfeeder", "nofeeder"))) %>%
  dplyr::select(camN, zone, `Feeder?`, value, sig)
ggplot(pvdf.long, mapping = aes(`Feeder?`, value, color = as.factor(zone), linetype = sig, group = camN, label = camN)) + 
  geom_point() + geom_line(linewidth = 0.6) + labs(y = 'Passage rate (# of passage/day/m)') +
  scale_color_discrete(name = "Zone") + scale_linetype_discrete(name = "significant?")


# Calculate percent reduction
pvdf2 %>% dplyr::select(camN, withfeeder, nofeeder, zone, EDD, t.sig, MultipleF) %>%
  mutate(percR = - ((withfeeder - nofeeder)/withfeeder) * 100)



