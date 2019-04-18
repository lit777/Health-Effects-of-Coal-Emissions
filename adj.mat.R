#### Create the matrix, W, for identifying adjacent zipcodes

library(tigris)
library(rgeos)
library(sp)

df <- zctas(cb = FALSE)
id <- gIntersects(df, df, byid = TRUE)

# read crosswalk between zip and zcta
crosswalk <- read.csv("Zip_to_ZCTA_Crosswalk_JSI2014.csv")
index <- NULL
zcta <- formatC(crosswalk$ZCTA_crosswalk, width = 5, format = "d", flag = "0")
zip <- formatC(crosswalk$ZIP, width = 5, format = "d", flag = "0")
temp_id <- Data2$zip

for(i in 1:12370){
  index[i] <- zcta[which(zip==temp_id[i])]
}

index1 <- NULL
temp_id <- Data2$zip
for(i in 1:12370){
  index1[i] <- which(df$GEOID10 == index[i])
}

adj.mat <- id[index1,index1]

save(adj.mat, file="adj.mat.RData")
