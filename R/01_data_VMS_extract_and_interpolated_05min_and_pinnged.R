library(sp)
library(gisland)
library(dplyrOracle)
library(mar)
library(padr)
library(lubridate)
library(stringr)
library(multidplyr)
library(tidyverse)
attach("data/logbooks.rda")
attach("data/shapes.rda")

print("### Allocate processing among 35 cores")
n.clusters <- 35
cluster <- create_cluster(n.clusters)

GID <- Gear$gid[Gear$gclass %in% c(6, 9, 14, 15)]
VID <- sort(unique(Base$vid[Base$gclass %in% c(6, 9, 14, 15)]))
#YEARS <- 2016:2009
YEARS <- 2012:2009

db <- src_oracle("mar")

for(y in 1:length(YEARS)) {
  
  print(YEARS[y])
  
  toga_db <-
    afli_toga(db) %>%
    select(visir, ibotni, togtimi)
  
  tows <-
    afli_stofn(db) %>%
    filter(ar %in% YEARS[y],
           veidarf %in% GID) %>%
    select(visir, skipnr, veidarf, vedags) %>%
    left_join(toga_db) %>%
    rename(id = visir,
           gid = veidarf,
           vid = skipnr) %>%
    collect(n = Inf) %>%
    mutate(ibotni = str_pad(ibotni, 4, "left", pad = "0"),
           ibotni = paste0(str_sub(ibotni, 1, 2),
                           ":",
                           str_sub(ibotni, 3, 4)),
           t1 = ymd_hm(paste(as.character(vedags), ibotni)),
           t2 = t1 + minutes(togtimi)) %>%
    select(t1, t2, vid, gid, id) %>%
    filter(!is.na(t2))
  
  VMS <-
    tbl_mar(db, "stk.stk_vms_v") %>%
    mutate(year = to_number(to_char(posdate, 'YYYY')),
           lon = poslon * 45 / atan(1),
           lat = poslat * 45 / atan(1),
           heading = heading * 45 / atan(1),
           speed = speed * 1.852) %>%
    filter(year == YEARS[y],
           between(lon, -179.9, 179.9),
           between(lat, -89.9, 89.9)) %>%
    select(year, mobileid, vid = skip_nr, date = posdate, lon, lat, speed, heading, in_out_of_harbor, harborid) %>%
    collect(n = Inf) %>%
    mutate(vid = as.integer(vid)) %>%
    filter(!is.na(vid),
           vid %in% VID) %>% 
           #is.na(in_out_of_harbor)) %>% 
    arrange(vid, date) %>% 
    distinct(date, vid, lon, lat, .keep_all = TRUE) %>% 
    group_by(vid) %>% 
    mutate(p.interval = (date - lag(date))/dminutes(1),
           p.interval = ifelse(is.na(p.interval), lead(p.interval), p.interval)) %>% 
    ungroup() #%>% 
    #select(-c(in_out_of_harbor, harborid))
  
  # exclude vessel with only one ping
  vessel.with.two.pings <-
    VMS %>% 
    group_by(vid) %>% 
    summarise(n = n()) %>% 
    arrange(n) %>% 
    filter(n <= 2)
  VMS <-
    VMS %>%
    filter(!vid %in% vessel.with.two.pings$vid) %>% 
    mutate(p = 1:n())
  
  
  
  VID <- unique(VMS$vid)
  
  res <- list()
  
  for (v in 1:length(VID)) {
    
    vms <- 
      VMS %>% filter(vid == VID[v])
    
    print(paste("VID:",c(VID[v]), "Pings:", nrow(vms)))
             # check only intervals within enpoints
    if(any(vms$p.interval[2:(nrow(vms)-1)] > 12)) {
      
      counter <- 0
      res.tmp <- list()
      for (i in 2:(nrow(vms) - 1)) {
        #print(i)
                                        
        if(vms$p.interval[i + 1] > 9 & 
           # pings not in harbour
           (is.na(vms$in_out_of_harbor[i]) | is.na(vms$in_out_of_harbor[i + 1]))) {
          
          # print(paste0("v: ", v, ", i: ", i)) 
          counter <- counter + 1
          x <- vms %>% slice(i:(i + 1))
          x2 <-
            x %>% 
            pad(by = "date", interval = "6 min") %>% 
            mutate(p.intpl = TRUE) %>% 
            bind_rows(x) %>% 
            mutate(p.intpl = ifelse(is.na(p.intpl), FALSE, p.intpl)) %>% 
            arrange(date, p.intpl) %>%
            fill(p) %>% 
            distinct(date, lon, lat, .keep_all = TRUE)
          n <- nrow(x2)
          x2$lon[2:(n - 1)] <- seq(x2$lon[1], x2$lon[n], length.out = n)[2:(n-1)]
          x2$lat[2:(n - 1)] <- seq(x2$lat[1], x2$lat[n], length.out = n)[2:(n-1)]
          x2$speed[2:(n - 1)] <- seq(x2$speed[1], x2$speed[n], length.out = n)[2:(n-1)]
          
          res.tmp[[counter]] <-
            x2 %>% 
            fill(year:vid)
        } # interpolation loop
        
      } # next point sets
      
      res.tmp <-
        bind_rows(res.tmp)
      
      res[[v]] <-
        vms %>% 
        filter(!p %in% unique(res.tmp$p)) %>% 
        bind_rows(res.tmp) %>% 
        arrange(vid, date) #%>% 
        #mutate(p.intpl = ifelse(is.na(p.intpl), FALSE, p.intpl))
      
    } else { # If all less than 15 minutes, pass the original
      
      
      res[[v]] <- vms #%>% mutate(p.intpl = FALSE)
      
    }
  } # vessel loop
  
  vms <-
    bind_rows(res) %>% 
    mutate(p.intpl = ifelse(is.na(p.intpl), FALSE, TRUE)) %>% 
    arrange(vid, date) %>% 
    mutate(id = lookup_interval_ids(date, vid, tows, cn = c("t1", "t2", "vid", "id")))
  
  print("done linking to tows")
  #save(vms, file = paste0("data/vms_interpolate_", YEARS[y], ".rda"))
  
  # ----------------------------------------------------------------------------
  # Now for the ping classification
  group <- rep(1:n.clusters, length.out = nrow(vms))
  vms <- bind_cols(tibble(group), vms)
  
  # the data partitioning
  vms <-
    vms %>%
    partition(group, cluster = cluster)
  
  # the objects (can this be done outside the TWO-looop????)
  cluster_library(vms, "gisland")
  cluster_library(vms, "sp")
  cluster_assign_value(vms, "fao_27",      fao_27)
  cluster_assign_value(vms, "ospar",       ospar)
  cluster_assign_value(vms, "harbours.sp", harbours.sp)
  cluster_assign_value(vms, "coast",       coast)
  cluster_assign_value(vms, "regl",        regl)
  cluster_assign_value(vms, "skip3",       skip3)
  print("starting")
  
  vms <-
    vms %>%
    mutate(p.ices = geo_inside(lon, lat, fao_27))
  print("done fao")
  
  vms <-
    vms %>%
    mutate(p.ospr = geo_inside(lon, lat, ospar))
  print("done ospar")
  
  vms <-
    vms %>%
    mutate(p.harb = geo_inside(lon, lat, harbours.sp))
  print("done harbour")
  
  vms <-
    vms %>%
    mutate(p.land = geo_inside(lon, lat, coast))
  print("done coast")
  
  vms <-
    vms %>%
    mutate(p.eez = geo_inside(lon, lat, eez))
  print("done eez")
  
  vms <-
    vms %>%
    mutate(p.skip3 = geo_inside(lon, lat, skip3))
  print("done skipaflokkur 3")
  
  vms <-
    vms %>%
    mutate(p.regl = geo_inside(lon, lat, regl))
  print("done reglugerdir")
  
  vms <-
    vms %>%
    collect() %>%
    ungroup()
  print("done collecting")
  
  save(vms, file = paste0("data/vms_05min_", YEARS[y], ".rda"))
  
}
