
<!-- README.md is generated from README.Rmd. Please edit that file -->

# trafficEstimatr

<!-- badges: start -->

<!-- badges: end -->

The goal of trafficEstimatr is to
…

``` r
u = "http://data.dft.gov.uk/road-traffic/dft_traffic_counts_raw_counts.zip"
f = tempfile()
download.file(url = u, destfile = f)
unzip(zipfile = f, exdir = tempdir())
fd = list.files(tempdir(), pattern = "counts.csv", full.names = TRUE)
d = readr::read_csv(fd)
table(d$road_category)
table(d$road_name)
d$is_a_road = grepl(pattern = "^A", d$road_type)
table(d$road_type) / nrow(d)
da = d %>% 
  group_by(count_point_id) %>% 
  summarise(
    pcu = mean(all_motor_vehicles),
    lon = mean(longitude),
    lat = mean(latitude)
    )
di = d %>% filter(local_authority_name == "Isle of Wight")
dai = di %>% 
  group_by(count_point_id, longitude, latitude) %>% 
  summarise(
    n = n(),
    pcu = mean(all_motor_vehicles),
    lon = mean(longitude),
    lat = mean(latitude)
    ) %>% 
  ungroup()
dai %>% arrange(desc(n)) %>% slice(1:9)
dasf = sf::st_as_sf(dai, coords = c("lon", "lat"), crs = 4326)
dasf_iow = dasf[pct::pct_regions %>% filter(region_name == "isle-of-wight"), ]
tm_shape(dasf_iow) + tm_dots(size = "pcu")
saveRDS(dasf_iow, "dasf_iow.Rds")
# aim: try to estimate PCU values across IoW, if it works, try for all of UK (big data)
```

``` r
dasf_iow = readRDS("dasf_iow.Rds")
tm_shape(dasf_iow) + tm_dots(size = "pcu")
#> Linking to GEOS 3.5.1, GDAL 2.1.2, PROJ 4.9.3
#> Legend for symbol sizes not available in view mode.
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
