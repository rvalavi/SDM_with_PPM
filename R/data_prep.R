library(sf)
library(dplyr)

# load Fithian et al (2015) species data
load("data/moddat.RData")

lif <- ls(pattern = "moddat")
lif
for(i in 1:length(lif)){
  get(lif[i]) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
    st_transform(crs = 32756) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    mutate(po = 1) %>% 
    write.csv(paste0("data/species/", substr(lif[i], 1, nchar(lif[i]) - 7), ".csv"), row.names = FALSE)
  print(i)
}


st_read("data/towns/ecologist_towns.shp", crs = 4236) %>% 
  st_transform(crs = 32756) %>% 
  st_write("data/towns.shp")
