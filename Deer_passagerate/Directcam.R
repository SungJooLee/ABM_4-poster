my_packages <- c("tidyr", "dplyr", "ggplot2", "lubridate")
lapply(my_packages, library, character.only = TRUE)

setwd("D:/0.Columbia/0. Data/2.CDC_Suffolk/_Passage rate 2023")
cd1 <- read.csv("Direct cams/CD1.csv")

cd1 <- cd1 %>% dplyr::select(Date, Time, L_species, L_Sex, R_Species, R_Sex) %>% filter(L_species != "" | R_Species != "") %>% 
  mutate(DT = paste(Date, Time)) %>% mutate(DT = as.POSIXct(DT),
                                            Time = hms(Time))

# Deer density histogram analysis ---------------------------------------------------------------------------------------------
cd1.deer <- cd1 %>% filter(L_species == "Deer" | R_Species == "Deer")
cd1.deer %>% mutate(LS = recode(L_Sex, 'Female' = 1, 'Male' = 10),
                    RS = recode(R_Sex, 'Female' = 1, 'Male' = 10),
                    LS2 = replace_na(LS, 0),
                    RS2 = replace_na(RS, 0),
                    Sum = LS2 + RS2,
                    Type = ifelse(Sum == 2, "Female", ifelse(Sum == 1, "Female", 
                              ifelse(Sum == 10, "Male", ifelse(Sum == 11, "Mixed", 
                                ifelse(Sum == 20, "Male", "IDK")))))) -> cd1.deer

ggplot(cd1.deer %>% filter(Date == "2023-06-24"), mapping = aes(DT, after_stat(density), fill = Type)) +
  geom_density(alpha = 0.5) + 
  scale_x_datetime(breaks = "4 hour", date_labels =  "%H:%M")
  

ggplot(cd1 %>% filter(Date == "2023-06-22"), aes(x = DT, y = L_species)) + 
  ggridges::geom_density_ridges(fill = "#00AFBB", rel_min_height = 0.01, alpha = 0.5)

# Non-target hosts ---------------------------------------------------------------------------------------------4
df2 <- cd1 %>% mutate(Group = hour(Time)) %>% mutate(across(L_species:R_Sex, ~ na_if(.x, ""))) %>%
  mutate(L = ifelse(is.na(L_species), NA, ifelse(L_species == "Deer", paste(L_species, L_Sex), L_species)),
         R = ifelse(is.na(R_Species), NA, ifelse(R_Species == "Deer", paste(R_Species, R_Sex), R_Species)))
df2 <- df2 %>% group_by(Date, Group) %>%
  summarise(DL = ifelse(all(is.na(table(L))), NA, names(table(L))[which.max(table(L))]),
            DR = ifelse(all(is.na(table(R))), NA, names(table(R))[which.max(table(R))])) %>%
  mutate(Date = as.Date(Date))

# tile plot
ggplot(df2 %>% dplyr::select(-DR) %>% na.omit) + geom_tile(mapping = aes(x = Group, y = Date, fill = DL)) +
  scale_fill_manual(values = c("Chipmunk" = "yellow", "Deer Female" =  "brown", "Deer Male" = "darkblue", "Geese" = "pink","Racoon"= "darkgray", "Turkey" = "lightgreen"))+
  scale_x_continuous(0:23, name = "Hour of day", breaks = seq(0, 23, 2)) + scale_y_date() + theme_bw() + labs(fill = "Left")
ggplot(df2 %>% dplyr::select(-DL) %>% na.omit) + geom_tile(mapping = aes(x = Group, y = Date, fill = DR)) +
  scale_fill_manual(values = c("Chipmunk" = "yellow", "Deer Female" =  "brown", "Deer Male" = "darkblue", "Geese" = "pink","Racoon"= "darkgray", "Turkey" = "lightgreen"))+
  scale_x_continuous(0:23, name = "Hour of day", breaks = seq(0, 23, 2)) + scale_y_date() + theme_bw() + labs(fill = "Right")

# pie chart
ggplot(df2 %>% filter(!is.na(DL))) + geom_bar(mapping = aes(x = "", fill = DL)) + coord_polar(theta='y') + labs(fill = "Left") + 
  scale_fill_manual(values = c("Chipmunk" = "yellow", "Deer Female" =  "brown", "Deer Male" = "darkblue", "Geese" = "pink","Racoon"= "darkgray", "Turkey" = "lightgreen"))
ggplot(df2 %>% filter(!is.na(DR))) + geom_bar(mapping = aes(x = "", fill = DR)) + coord_polar(theta='y') + labs(fill = "Right")  + 
  scale_fill_manual(values = c("Chipmunk" = "yellow", "Deer Female" =  "brown", "Deer Male" = "darkblue", "Geese" = "pink","Racoon"= "darkgray", "Turkey" = "lightgreen"))





