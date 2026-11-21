output_dir <- "D:/github/ABM_4-poster/Results.climate"
weather_daily <- read.csv(file.path(output_dir, "halsey_weather_daily.csv"))
code <- readLines("R_codes/code03.Halsey_climate.R")
start <- grep("^seasonal_weather <-", code)
end <- which(seq_along(code) > start & code == "dev.off()")[1]
eval(parse(text = code[start:end]))