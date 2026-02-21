# Remote Sensing Applications: Group Glacial Retreat
Justin Lingg-Laham, Liev Paustian, Lisa Wohlgemuth

## Introduction

Remote Sensing Applications: Group Glacial Retreat Part 2: statistical
Analysis in R

This R script works on the post-processing, Glacial conditions, climate
variables and Hypsometry for the Fox Glacier in New Zealand

Inputs (from GEE): - fox_glacier.csv (Part 1.2) -
Fox_Climate_ERA5Land_1989_2025.csv (Part 1.3) - Fox_hypsometry_bins.csv
(Part 1.4)

## Preparation

To ensure Reproducibility, it is crucial to understand that when the
document is rendered, the working directory is automatically set to the
folder in which it is located. Therefore, it is necessary to maintain
the attached folder structure! The data import takes place via relative
paths to be independent from local storage pathways from different
computersystems.

``` r
getwd()
```

    [1] "C:/Users/justi/OneDrive/Desktop/UGM Semester 1/remote_sensing_applications/storyMap_glacier"

## Librarys

The following Libraries are used to provide access to additional
functions:

``` r
library(tidyverse)
library(lubridate)
library(MASS)   
library(ggrepel)
library(scales)
```

# 1. Data

``` r
#Load Data
glacier_30 <- read.csv("./data/raw/fox_glacier.csv", sep = ",", header = TRUE) # Satellite-Data of Fox-Glacier
clim <- read.csv("./data/raw/Fox_Climate_ERA5Land_1989_2025.csv", sep = ",", header = TRUE) # ERA5-land climate data

#Clean glacier-data ((first steps in GEE))
df2 <- glacier_30 %>%
  mutate(
    year = as.integer(year),
    doSnowline = doSnowline %in% c(TRUE, "TRUE", 1, "1"),
    ok_overall = ok_overall %in% c(TRUE, "TRUE", 1, "1"),
    ok_clear   = ok_clear   %in% c(TRUE, "TRUE", 1, "1")
  )

# Snowline subset (Define a Quality-Threshold)
snow_cf_thr  <- 0.20 #at least 20% clear pixel for snowline
snow_min_km2 <- 5 #at least 5 km2 glacier-area

#produce quality-snowline data
snow_df <- df2 %>%
  filter(
    doSnowline == TRUE,
    !is.na(snowline_edge_m),
    nSnow >= 3,
    clearFraction_snow >= snow_cf_thr,
    !is.na(glacier_km2),
    glacier_km2 >= snow_min_km2
  )

# Area subset (Quality-trheshold)
area_min_km2 <- 10 #at least 10 km2 to be area

area_strict <- df2 %>%
  filter(
    ok_overall == TRUE,
    !is.na(glacier_km2),
    imgCount >= 5,
    clearFraction >= 0.20,
    glacier_km2 >= area_min_km2
  )

#Output Statistics
cat("Years total:", nrow(df2), "\n") 
```

    Years total: 37 

``` r
cat("Snowline clean years:", nrow(snow_df), " (unique:", n_distinct(snow_df$year), ")\n")
```

    Snowline clean years: 26  (unique: 26 )

``` r
cat("Area STRICT years:", nrow(area_strict), " (unique:", n_distinct(area_strict$year), ")\n")
```

    Area STRICT years: 28  (unique: 28 )

``` r
head(glacier_30, 15)
```

       system.index LST_median_C clearFraction clearFraction_snow             dem
    1             0   -2.6655239     0.9110165                 NA NASADEM_HGT_001
    2             1   -1.3791870     0.8288045         0.09540215 NASADEM_HGT_001
    3             2           NA            NA                 NA NASADEM_HGT_001
    4             3    0.1883497     0.5972790                 NA NASADEM_HGT_001
    5             4           NA            NA         0.52652108 NASADEM_HGT_001
    6             5           NA            NA                 NA NASADEM_HGT_001
    7             6           NA            NA                 NA NASADEM_HGT_001
    8             7           NA            NA                 NA NASADEM_HGT_001
    9             8           NA            NA                 NA NASADEM_HGT_001
    10            9           NA            NA                 NA NASADEM_HGT_001
    11           10    0.1141623     0.9926243                 NA NASADEM_HGT_001
    12           11   -1.6247016     0.9947841         0.95314867 NASADEM_HGT_001
    13           12   -0.1235998     0.9948972         0.95734230 NASADEM_HGT_001
    14           13   -0.1190369     0.9930122         0.94625484 NASADEM_HGT_001
    15           14   -3.8138686     0.9948972         0.94847206 NASADEM_HGT_001
       doSnowline  firstDate glacier_km2 imgCount   lastDate minClearFraction
    1           0 1989-11-09    54.59791        6 1990-02-13              0.1
    2           0 1990-12-14    51.72279        2 1990-12-30              0.1
    3           0                     NA        0                         0.1
    4           0 1992-12-19    34.01140        3 1993-02-21              0.1
    5           0                     NA        0                         0.1
    6           0                     NA        0                         0.1
    7           0                     NA        0                         0.1
    8           0                     NA        0                         0.1
    9           0                     NA        0                         0.1
    10          0                     NA        0                         0.1
    11          0 1999-08-01    58.53726       19 2000-05-31              0.1
    12          1 2000-08-03    61.03747       19 2001-05-18              0.1
    13          1 2001-08-22    58.14506       16 2002-05-21              0.1
    14          1 2002-08-09    60.99791       19 2003-05-24              0.1
    15          1 2003-08-04    62.50721       33 2004-05-26              0.1
       minClearFraction_snowline minGlacierKm2_forSnowline minScenes_quality
    1                        0.2                         5                 2
    2                        0.2                         5                 2
    3                        0.2                         5                 2
    4                        0.2                         5                 2
    5                        0.2                         5                 2
    6                        0.2                         5                 2
    7                        0.2                         5                 2
    8                        0.2                         5                 2
    9                        0.2                         5                 2
    10                       0.2                         5                 2
    11                       0.2                         5                 2
    12                       0.2                         5                 2
    13                       0.2                         5                 2
    14                       0.2                         5                 2
    15                       0.2                         5                 2
       minScenes_snowline nSnow ndsiThr ok_clear ok_overall ok_scenes scale_m
    1                   3     0     0.4        1          1         1      30
    2                   3     3     0.4        1          1         1      30
    3                   3     0     0.4    false      false         0      30
    4                   3     0     0.4        1          1         1      30
    5                   3     2     0.4    false      false         0      30
    6                   3     0     0.4    false      false         0      30
    7                   3     0     0.4    false      false         0      30
    8                   3     0     0.4    false      false         0      30
    9                   3     0     0.4    false      false         0      30
    10                  3     0     0.4    false      false         0      30
    11                  3     0     0.4        1          1         1      30
    12                  3     6     0.4        1          1         1      30
    13                  3     6     0.4        1          1         1      30
    14                  3     5     0.4        1          1         1      30
    15                  3     5     0.4        1          1         1      30
         season snowElev_p10_m snowElev_p50_m snowElev_p90_m snowWindow
    1  8/1–5/31             NA             NA             NA   1/1–3/31
    2  8/1–5/31             NA             NA             NA   1/1–3/31
    3  8/1–5/31             NA             NA             NA   1/1–3/31
    4  8/1–5/31             NA             NA             NA   1/1–3/31
    5  8/1–5/31             NA             NA             NA   1/1–3/31
    6  8/1–5/31             NA             NA             NA   1/1–3/31
    7  8/1–5/31             NA             NA             NA   1/1–3/31
    8  8/1–5/31             NA             NA             NA   1/1–3/31
    9  8/1–5/31             NA             NA             NA   1/1–3/31
    10 8/1–5/31             NA             NA             NA   1/1–3/31
    11 8/1–5/31             NA             NA             NA   1/1–3/31
    12 8/1–5/31       1735.667       2279.669       2759.372   1/1–3/31
    13 8/1–5/31       1703.758       2263.553       2743.553   1/1–3/31
    14 8/1–5/31       1752.184       2279.684       2743.477   1/1–3/31
    15 8/1–5/31       1735.636       2279.669       2759.465   1/1–3/31
       snowline_edge_m swir1Max windowUsed year
    1               NA      0.2     SEASON 1989
    2               NA      0.2     SEASON 1990
    3               NA      0.2     SEASON 1991
    4               NA      0.2     SEASON 1992
    5               NA      0.2     SEASON 1993
    6               NA      0.2     SEASON 1994
    7               NA      0.2     SEASON 1995
    8               NA      0.2     SEASON 1996
    9               NA      0.2     SEASON 1997
    10              NA      0.2     SEASON 1998
    11              NA      0.2     SEASON 1999
    12        1799.779      0.2     SEASON 2000
    13        2040.116      0.2     SEASON 2001
    14        2359.085      0.2     SEASON 2002
    15        1847.977      0.2     SEASON 2003
                                         .geo
    1  {"type":"MultiPoint","coordinates":[]}
    2  {"type":"MultiPoint","coordinates":[]}
    3  {"type":"MultiPoint","coordinates":[]}
    4  {"type":"MultiPoint","coordinates":[]}
    5  {"type":"MultiPoint","coordinates":[]}
    6  {"type":"MultiPoint","coordinates":[]}
    7  {"type":"MultiPoint","coordinates":[]}
    8  {"type":"MultiPoint","coordinates":[]}
    9  {"type":"MultiPoint","coordinates":[]}
    10 {"type":"MultiPoint","coordinates":[]}
    11 {"type":"MultiPoint","coordinates":[]}
    12 {"type":"MultiPoint","coordinates":[]}
    13 {"type":"MultiPoint","coordinates":[]}
    14 {"type":"MultiPoint","coordinates":[]}
    15 {"type":"MultiPoint","coordinates":[]}

# 2. Statistics: Trends, correlation (Area vs. Snowline)

``` r
#Trends 
snow_lm <- lm(snowline_edge_m ~ year, data = snow_df) #Snowline (linear model with Snow-height per year)
area_lm_strict <- lm(glacier_km2 ~ year, data = area_strict) #Area (linear model with glacier area per year)

print(summary(snow_lm))
```


    Call:
    lm(formula = snowline_edge_m ~ year, data = snow_df)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -466.81  -59.59   17.82   88.88  346.49 

    Coefficients:
                  Estimate Std. Error t value Pr(>|t|)    
    (Intercept) -38763.865   9155.407  -4.234 0.000291 ***
    year            20.368      4.549   4.477 0.000157 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 174 on 24 degrees of freedom
    Multiple R-squared:  0.4551,    Adjusted R-squared:  0.4324 
    F-statistic: 20.05 on 1 and 24 DF,  p-value: 0.000157

``` r
print(summary(area_lm_strict))
```


    Call:
    lm(formula = glacier_km2 ~ year, data = area_strict)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -18.545  -1.257   1.588   3.241   4.822 

    Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
    (Intercept) 789.4964   225.3860   3.503  0.00168 **
    year         -0.3640     0.1121  -3.248  0.00320 **
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 5.194 on 26 degrees of freedom
    Multiple R-squared:  0.2886,    Adjusted R-squared:  0.2613 
    F-statistic: 10.55 on 1 and 26 DF,  p-value: 0.003198

``` r
# Mechanism: area vs snowline (correlation between both)
rel_lm  <- lm(glacier_km2 ~ snowline_edge_m, data = snow_df)
rel_rob <- rlm(glacier_km2 ~ snowline_edge_m, data = snow_df)

print(summary(rel_lm))
```


    Call:
    lm(formula = glacier_km2 ~ snowline_edge_m, data = snow_df)

    Residuals:
         Min       1Q   Median       3Q      Max 
    -21.9781  -1.3023   0.8466   2.9671   6.2804 

    Coefficients:
                     Estimate Std. Error t value Pr(>|t|)    
    (Intercept)     82.364474  11.263818   7.312 1.49e-07 ***
    snowline_edge_m -0.011162   0.005033  -2.218   0.0363 *  
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 5.811 on 24 degrees of freedom
    Multiple R-squared:  0.1701,    Adjusted R-squared:  0.1355 
    F-statistic: 4.918 on 1 and 24 DF,  p-value: 0.0363

``` r
print(summary(rel_rob))
```


    Call: rlm(formula = glacier_km2 ~ snowline_edge_m, data = snow_df)
    Residuals:
         Min       1Q   Median       3Q      Max 
    -23.6238  -2.0034   0.6167   1.5450   4.5778 

    Coefficients:
                    Value   Std. Error t value
    (Intercept)     75.6886  6.4878    11.6663
    snowline_edge_m -0.0077  0.0029    -2.6419

    Residual standard error: 2.928 on 24 degrees of freedom

``` r
# Correlations
cor_test_p <- cor.test(snow_df$glacier_km2, snow_df$snowline_edge_m, method = "pearson") #Pearson correlation: linear correlation
cor_test_s <- cor.test(snow_df$glacier_km2, snow_df$snowline_edge_m, method = "spearman") #Spearman correlation: monotone correlation

print(cor_test_p)
```


        Pearson's product-moment correlation

    data:  snow_df$glacier_km2 and snow_df$snowline_edge_m
    t = -2.2177, df = 24, p-value = 0.0363
    alternative hypothesis: true correlation is not equal to 0
    95 percent confidence interval:
     -0.68958792 -0.02979931
    sample estimates:
           cor 
    -0.4123911 

``` r
print(cor_test_s)
```


        Spearman's rank correlation rho

    data:  snow_df$glacier_km2 and snow_df$snowline_edge_m
    S = 4630, p-value = 0.00214
    alternative hypothesis: true rho is not equal to 0
    sample estimates:
          rho 
    -0.582906 

# 3. Plots (Glacier)

``` r
# Snowline-TimeSeries (January-March)
snow_fit <- lm(snowline_edge_m ~ year, data = snow_df)
ann <- broom::glance(snow_fit) %>%
  transmute( r2 = round(r.squared, 2), p  = signif(p.value, 2))
slope <- broom::tidy(snow_fit) %>%
  filter(term == "year") %>%
  pull(estimate)
label_txt <- paste0(
  "Trend: ", round(slope, 1), " m/Year",
  " | R² = ", ann$r2,
  " | p = ", ann$p
)
p1 <- ggplot(snow_df, aes(year, snowline_edge_m)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#00B8D9") +
  annotate(
  "label",
  x = max(snow_df$year, na.rm = TRUE) - 1,
  y = min(snow_df$snowline_edge_m, na.rm = TRUE) + 40,
  label = label_txt,
  hjust = 1, vjust = 0,
  label.size = 0,
  fill = "white",
  alpha = 0.9
) +
  labs(
    title = "TimeSeries: Snowline (snowmask-edge) for the Fox Glacier",
    subtitle = "Yearly Snowline-height [m] (Jan–Mar)",
    x = "Year",
    y = "Snowline elevation [m]"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Area-TimeSeries (August-May)
area_fit <- lm(glacier_km2 ~ year, data = area_strict)
gl <- broom::glance(area_fit)
slope <- broom::tidy(area_fit) %>% filter(term == "year") %>% pull(estimate)

label_txt <- paste0(
  "Trend: ", round(slope, 2), " km²/Year",
  " | R² = ", round(gl$r.squared, 2),
  " | p = ", signif(gl$p.value, 2)
)
p2_strict <- ggplot(area_strict, aes(year, glacier_km2)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "#00B8D9",) +
  labs(
    title = "TimeSeries: Glacier-area of the Fox Glacier",
    subtitle = "Yearly Area [km²] (Aug-May)",
    x = "Year",
    y = "Area (km²)"
  )+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank())

# Correlatiion: Area vs. Snowline
# Robust model
rel_rob <- rlm(glacier_km2 ~ snowline_edge_m, data = snow_df)
# Slope
rob_slope <- broom::tidy(rel_rob) %>%
  filter(term == "snowline_edge_m") %>%
  pull(estimate)
# correlation
cor_p <- cor.test(snow_df$glacier_km2, snow_df$snowline_edge_m, method = "pearson")
cor_s <- cor.test(snow_df$glacier_km2, snow_df$snowline_edge_m, method = "spearman")

label_txt <- paste0(
  "Robust slope: ", round(rob_slope, 3), " km² per m\n",
  "Pearson r = ", round(unname(cor_p$estimate), 2), ", p = ", signif(cor_p$p.value, 2), "\n",
  "Spearman ρ = ", round(unname(cor_s$estimate), 2), ", p = ", signif(cor_s$p.value, 2)
)
p3_rob <- ggplot(snow_df, aes(snowline_edge_m, glacier_km2)) + geom_point() +geom_smooth(method = "rlm", se = FALSE, color = "#00B8D9") + annotate(
    "label",
    x = max(snow_df$snowline_edge_m, na.rm = TRUE),
    y = min(snow_df$glacier_km2, na.rm = TRUE) + 1.2,
    label = label_txt,
    hjust = 1, vjust = 0,
    label.size = 0,
    fill = "white",
    alpha = 0.9) +
  labs(
    title = "Correlation: Area vs. Snowline for Fox Glacier",
    subtitle = "Higher snowline years tend to have smaller area ",
    x = "Snowline elevation [m]",
    y = "Area (km²)"
  ) +
  theme_minimal()+theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

# save and export
ggsave("figures/plots_glacier/fox_snowline_timeseries.png", p1, width = 9, height = 5, dpi = 200)      
ggsave("figures/plots_glacier/fox_area_timeseries_strict.png", p2_strict, width = 9, height = 5, dpi = 200) 
ggsave("figures/plots_glacier/fox_area_vs_snowline_robust.png", p3_rob, width = 7, height = 6, dpi = 200)   

p1
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
p2_strict
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-6-1.png)

``` r
p3_rob
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-7-1.png)

# 4. Climate Analysis (Preciptation, Temperature) ERA5

``` r
#prepare the data
clim2 <- clim %>%
  mutate(year = as.integer(year)) %>%
  filter(nMonths_season > 0, nMonths_summer > 0)
#join climate data to glacier
df <- df2 %>%
  dplyr::left_join(
    clim2 %>%
      dplyr::select(
        year,
        t2m_C_mean_season, tp_mm_sum_season,
        t2m_C_mean_summer, tp_mm_sum_summer,
        t2m_C_pos_sum_summer,
        nMonths_season, nMonths_summer
      ),
    by = "year"
  )

cat("Joined years:", nrow(df), "\n")
```

    Joined years: 37 

``` r
cat("Years with climate:", sum(!is.na(df$t2m_C_mean_summer)), "\n")
```

    Years with climate: 34 

``` r
# recreate subsets WITH climate availability
snow_df_c <- df %>%
  filter(
    doSnowline == TRUE,
    !is.na(snowline_edge_m),
    nSnow >= 3,
    clearFraction_snow >= snow_cf_thr,
    !is.na(glacier_km2),
    glacier_km2 >= snow_min_km2,
    !is.na(t2m_C_mean_summer)
  )

area_strict_c <- df %>%
  filter(
    ok_overall == TRUE,
    !is.na(glacier_km2),
    imgCount >= 5,
    clearFraction >= 0.20,
    glacier_km2 >= area_min_km2,
    !is.na(tp_mm_sum_season),
    !is.na(t2m_C_mean_season)
  )
#output
cat("Snow years (strict + climate):", nrow(snow_df_c), "\n")
```

    Snow years (strict + climate): 23 

``` r
cat("Area years (strict + climate):", nrow(area_strict_c), "\n")
```

    Area years (strict + climate): 25 

# 5. Plot: Climate

``` r
clim_year <- df %>%
  filter(!is.na(t2m_C_mean_summer)) %>%
  group_by(year) %>%
  summarise(t2m_C_mean_summer = mean(t2m_C_mean_summer), .groups = "drop")

clim_fit <- lm(t2m_C_mean_summer ~ year, data = clim_year)
gl <- broom::glance(clim_fit)
slope <- broom::tidy(clim_fit) %>% filter(term == "year") %>% pull(estimate)

label_txt <- paste0(
  "Trend: ", round(slope, 3), " °C/Jahr",
  " | R² = ", round(gl$r.squared, 2),
  " | p = ", signif(gl$p.value, 2)
)

p_clim_ts_up <- ggplot(clim_year, aes(year, t2m_C_mean_summer)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(
    method = "lm", se = TRUE,
    color = "#00B8D9",
  ) +
  annotate(
    "label",
    x = max(clim_year$year, na.rm = TRUE) - 1,
    y = min(clim_year$t2m_C_mean_summer, na.rm = TRUE) + 0.15,
    label = label_txt,
    hjust = 1, vjust = 0,
    label.size = 0,
    fill = "white",
    alpha = 0.9
  ) +
  labs(
    title = "TimeSeries: Summer temperature for Fox Glacier",
    subtitle = "Annual mean 2 m temperature (Jan–Mar) from ERA5-Land",
    x = "Year",
    y = "T2m mean [°C]",) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p_clim_ts_up
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-9-1.png)

``` r
ggsave("figures/plots_climate/fox_summerTemp_timeseries.png", p_clim_ts_up, width = 8, height = 5, dpi = 200)
```

# 6. Plot: Preciptation + Temperature = Snow-proxy (cold-weighted preciptation)

``` r
# precip only helps if it's cold enough

T0 <- 1.0  # °C threshold (snow to rain)
kW <- 1.5  # smoothness (bigger = flat translation, low = scharper)

area_snowproxy <- area_strict_c %>%
  mutate(
    snow_weight = 1 / (1 + exp((t2m_C_mean_season - T0) / kW)),   # if t2m_C_mean_season colder then T0 = expression(negative) low = 1, if warmer expression(positive) = 0
    tp_snowproxy = tp_mm_sum_season * snow_weight
  )#snow proxy: much cold preciptation = high tp_snowproxy (high accumulation of snow), much warm preciptation = tp_snowproxy low (rain, no development for glacier)

m_area_snowproxy <- lm(glacier_km2 ~ tp_snowproxy, data = area_snowproxy)
print(summary(m_area_snowproxy))
```


    Call:
    lm(formula = glacier_km2 ~ tp_snowproxy, data = area_snowproxy)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -5.1912 -0.9097  0.6295  1.3423  4.2875 

    Coefficients:
                 Estimate Std. Error t value Pr(>|t|)    
    (Intercept)   56.3218     1.3410   42.00   <2e-16 ***
    tp_snowproxy   0.3128     0.1378    2.27   0.0329 *  
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 2.549 on 23 degrees of freedom
    Multiple R-squared:  0.1831,    Adjusted R-squared:  0.1475 
    F-statistic: 5.154 on 1 and 23 DF,  p-value: 0.03288

``` r
m_area_snowproxy_T <- lm(glacier_km2 ~ tp_snowproxy + t2m_C_mean_season, data = area_snowproxy)
print(summary(m_area_snowproxy_T))
```


    Call:
    lm(formula = glacier_km2 ~ tp_snowproxy + t2m_C_mean_season, 
        data = area_snowproxy)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -5.4437 -0.6447  0.9703  1.5236  3.1459 

    Coefficients:
                      Estimate Std. Error t value Pr(>|t|)    
    (Intercept)        76.4623    11.5703   6.608 1.21e-06 ***
    tp_snowproxy       -0.2441     0.3443  -0.709   0.4858    
    t2m_C_mean_season  -3.3639     1.9205  -1.752   0.0938 .  
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 2.442 on 22 degrees of freedom
    Multiple R-squared:  0.283, Adjusted R-squared:  0.2179 
    F-statistic: 4.342 on 2 and 22 DF,  p-value: 0.02573

``` r
# plot
m_area_snowproxy <- lm(glacier_km2 ~ tp_snowproxy, data = area_snowproxy)
gl <- broom::glance(m_area_snowproxy)
slope <- broom::tidy(m_area_snowproxy) %>% filter(term == "tp_snowproxy") %>% pull(estimate)

label_txt <- paste0(
  "Slope: ", round(slope, 3), " km² per mm",
  " | R² = ", round(gl$r.squared, 2),
  " | p = ", signif(gl$p.value, 2)
)

p_area_snowproxy_col <- ggplot(
  area_snowproxy,
  aes(x = tp_snowproxy, y = glacier_km2, color = snow_weight)
) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_smooth(
    aes(group = 1),               
    method = "lm",
    se = TRUE,
    color = "#00B8D9",
  ) +
  scale_color_viridis_c(
  name = "Snow fraction\n(1=much snow-weight, 0=not much)",
  limits = quantile(area_snowproxy$snow_weight, c(0.05, 0.95), na.rm = TRUE),
  oob = scales::squish
) +
  labs(
  title = "Combination: Area vs snowfall-precipitation for the Fox Glacier (Aug–May)",
  subtitle = "Effectiveness by Preciptation for increasing Glacier-area of cold Seasons",
  x = "Cold-weighted precipitation [mm]",
  y = expression(paste("Area proxy [km"^2, "]")),
) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p_area_snowproxy_col
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
ggsave("figures/plots_climate/fox_area_vs_snowproxy.png", p_area_snowproxy_col, width = 8.5, height = 5.2, dpi = 200)
```

\#7. Plot: Anomaly TimeSeries

``` r
# Define Baseline 
base_start <- 1991
base_end   <- 2020

# Calculate Mean baseline
area_anom <- area_strict_c %>%
  mutate(in_base = year >= base_start & year <= base_end)

base_vals <- area_anom %>%
  filter(in_base) %>%
  summarise(
    tp_base   = mean(tp_mm_sum_season, na.rm = TRUE),
    t_base    = mean(t2m_C_mean_season, na.rm = TRUE),
    area_base = mean(glacier_km2, na.rm = TRUE)
  )

tp_base   <- base_vals$tp_base
t_base    <- base_vals$t_base
area_base <- base_vals$area_base

#calculate the Anomalie
area_anom <- area_anom %>%
  mutate(
    tp_anom   = tp_mm_sum_season - tp_base,
    t_anom    = t2m_C_mean_season - t_base,
    area_anom = glacier_km2 - area_base
  )

cat("Baseline:", base_start, "-", base_end, "\n")
```

    Baseline: 1991 - 2020 

``` r
print(base_vals)
```

       tp_base   t_base area_base
    1 95.25449 4.348737  59.59845

``` r
# Models on anomalies # sensitivity on temperature
m_anom_T   <- lm(area_anom ~ t_anom, data = area_anom)         # temp and precip together
m_anom_TP  <- lm(area_anom ~ t_anom + tp_anom, data = area_anom) 
# interation of both
m_anom_int <- lm(area_anom ~ t_anom * tp_anom, data = area_anom)

print(summary(m_anom_T))
```


    Call:
    lm(formula = area_anom ~ t_anom, data = area_anom)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -5.4982 -0.4850  0.6937  1.6377  3.6435 

    Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
    (Intercept)  -0.1483     0.4950  -0.300  0.76722   
    t_anom       -2.1063     0.7284  -2.892  0.00823 **
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 2.415 on 23 degrees of freedom
    Multiple R-squared:  0.2666,    Adjusted R-squared:  0.2348 
    F-statistic: 8.363 on 1 and 23 DF,  p-value: 0.008226

``` r
print(summary(m_anom_TP))
```


    Call:
    lm(formula = area_anom ~ t_anom + tp_anom, data = area_anom)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -5.6124 -1.9671  0.6853  1.4211  3.5703 

    Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
    (Intercept) -0.17662    0.48642  -0.363  0.72000   
    t_anom      -2.27974    0.72623  -3.139  0.00477 **
    tp_anom     -0.05468    0.04004  -1.366  0.18585   
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 2.371 on 22 degrees of freedom
    Multiple R-squared:  0.324, Adjusted R-squared:  0.2625 
    F-statistic: 5.271 on 2 and 22 DF,  p-value: 0.01348

``` r
print(summary(m_anom_int))
```


    Call:
    lm(formula = area_anom ~ t_anom * tp_anom, data = area_anom)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -5.5997 -1.1448  0.6153  1.3600  3.6002 

    Coefficients:
                   Estimate Std. Error t value Pr(>|t|)   
    (Intercept)    -0.15975    0.49581  -0.322  0.75049   
    t_anom         -2.20477    0.75277  -2.929  0.00802 **
    tp_anom        -0.06451    0.04495  -1.435  0.16599   
    t_anom:tp_anom  0.02446    0.04737   0.516  0.61096   
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 2.412 on 21 degrees of freedom
    Multiple R-squared:  0.3324,    Adjusted R-squared:  0.2371 
    F-statistic: 3.486 on 3 and 21 DF,  p-value: 0.03389

``` r
# Sensitivity automation (ΔArea per +1°C)
bT  <- coef(m_anom_T)[["t_anom"]]
ciT <- confint(m_anom_T)["t_anom", ]
cat("\nSENSITIVITY: per +1°C (Aug–May), area_anom changes by",
    round(bT, 2), "km² (95% CI:", round(ciT[1], 2), "to", round(ciT[2], 2), ")\n")
```


    SENSITIVITY: per +1°C (Aug–May), area_anom changes by -2.11 km² (95% CI: -3.61 to -0.6 )

``` r
# TimeSeries: anomalie
anom_fit <- lm(area_anom ~ year, data = area_anom)
gl <- broom::glance(anom_fit)
slope <- broom::tidy(anom_fit) %>% filter(term == "year") %>% pull(estimate)

label_txt <- paste0(
  "Trend: ", round(slope, 2), " km²/Year",
  " | R² = ", round(gl$r.squared, 2),
  " | p = ", signif(gl$p.value, 2)
)

p_anom_ts_up <- ggplot(area_anom, aes(x = year, y = area_anom)) +
  geom_hline(yintercept = 0, linewidth = 0.4, alpha = 0.7) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE,
    color = "#00B8D9") +
  annotate(
    "label",
    x = max(area_anom$year, na.rm = TRUE) - 1,
    y = min(area_anom$area_anom, na.rm = TRUE) + 0.6,
    label = label_txt,
    hjust = 1, vjust = 0,
    label.size = 0,
    fill = "white",
    alpha = 0.9
  ) +
  labs(
    title = paste0("Area anomaly of the Fox Glacier (1991-2020)"),
    subtitle = "when was the Glacier under the Reference-Niveau (Aug–May)",
    x = "Year",
    y = expression(paste("Area anomaly [km"^2, "]")),
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# save and export
ggsave("figures/plots_climate/fox_area_anomaly_timeseries.png", p_anom_ts_up, width = 8, height = 5, dpi = 200) 
 p_anom_ts_up
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
# Export csv (for the step 8)
write_csv(snow_df, "fox_snowline_clean.csv")
write_csv(area_strict, "fox_area_clean_strict.csv")
write_csv(area_anom, "fox_area_strict_with_anomalies.csv")
```

# 8. Hypsometrie of Snowarea with Snowline Overlay

``` r
# (Measurement and representation of the distribution of elevations and depressions on the Earth's surface relative to sea level)
#Where at the Glacier is Snow-Area missing

# We Need the GEE-Output: glacier_hypo
# 8. Hypsometry of snow/ice area with snowline overlay (robust)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

#data
hypso_path    <- "./data/raw/Fox_hypsometry_bins.csv"
snowline_path <- "fox_snowline_clean.csv"   # from earlier chapter

# lables
# Hypsometry file contain "Early/Late" = "pre-2006/post-2006"
map_to_period <- function(x){
  x <- as.character(x)
  case_when(
    x %in% c("Early", "early", "pre", "pre2006", "pre-2006", "Pre-2006") ~ "pre-2006",
    x %in% c("Late", "late", "post", "post2006", "post-2006", "Post-2006") ~ "post-2006",
    TRUE ~ NA_character_
  )
}

period_levels <- c("pre-2006", "post-2006")

# read data
hyp <- read_csv(hypso_path, show_col_types = FALSE) %>%
  mutate(
    year_group_raw = year_group,
    year_group = map_to_period(year_group),
    year_group = factor(year_group, levels = period_levels),
    elev_bin   = as.numeric(elev_bin),
    area_km2   = as.numeric(area_km2)
  ) %>%
  filter(!is.na(year_group), !is.na(elev_bin), !is.na(area_km2)) %>%
  arrange(year_group, elev_bin)

stopifnot(all(levels(droplevels(hyp$year_group)) %in% period_levels))

# snowline groups
snowline_groups <- NULL

if (file.exists(snowline_path)) {

  snow <- readr::read_csv(snowline_path, show_col_types = FALSE) %>%
    dplyr::transmute(
      year = as.integer(year),
      snowline_edge_m = as.numeric(snowline_edge_m)
    ) %>%
    dplyr::filter(!is.na(year), !is.na(snowline_edge_m))

  summarise_snow <- function(dat, y_min, y_max, label){
    out <- dat %>%
      dplyr::filter(year >= y_min, year <= y_max) %>%
      dplyr::summarise(
        n_years = dplyr::n(),
        snowline_mean_m = mean(snowline_edge_m, na.rm = TRUE),
        snowline_p10_m  = stats::quantile(snowline_edge_m, 0.10, na.rm = TRUE),
        snowline_p90_m  = stats::quantile(snowline_edge_m, 0.90, na.rm = TRUE),
        .groups = "drop"
      )
    out$year_group <- label  
    out
  }

  # default split at 2006
  pre  <- summarise_snow(snow, -Inf, 2005L, "pre-2006")
  post <- summarise_snow(snow, 2006L,  Inf, "post-2006")

  # fallback split if one side is empty
  if (pre$n_years == 0 || post$n_years == 0) {
    yrs <- sort(unique(snow$year))
    mid <- yrs[floor(length(yrs) / 2)]
    pre  <- summarise_snow(snow, min(yrs), mid, "pre-2006")
    post <- summarise_snow(snow, mid + 1L, max(yrs), "post-2006")
  }

  # only build if both groups have data
  if (pre$n_years > 0 && post$n_years > 0) {
    snowline_groups <- dplyr::bind_rows(pre, post) %>%
      dplyr::select(year_group, snowline_mean_m, snowline_p10_m, snowline_p90_m) %>%
      dplyr::mutate(
        year_group = factor(as.character(year_group), levels = period_levels)
      )
  } else {
    snowline_groups <- NULL
  }
}
# median
w_median <- function(x, w){
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[min(which(cw >= 0.5))]
}

# calculate metrics
metrics <- hyp %>%
  group_by(year_group) %>%
  summarise(
    total_area_km2 = sum(area_km2, na.rm = TRUE),
    elev_w_median  = w_median(elev_bin, area_km2),
    elev_w_mean    = sum(elev_bin * area_km2, na.rm = TRUE) / sum(area_km2, na.rm = TRUE),
    .groups = "drop"
  )
print(metrics)
```

    # A tibble: 2 × 4
      year_group total_area_km2 elev_w_median elev_w_mean
      <fct>               <dbl>         <dbl>       <dbl>
    1 pre-2006             59.7          2250       2155.
    2 post-2006            58.5          2200       2147.

``` r
# AAR-Proxy = Accumulation Area Ratio (How much snow stays on the glacier above the snowline)
aar <- NULL
if (!is.null(snowline_groups)) {
  aar <- hyp %>%
    dplyr::left_join(
      dplyr::select(snowline_groups, year_group, snowline_mean_m),
      by = "year_group"
    ) %>%
    dplyr::group_by(year_group) %>%
    dplyr::summarise(
      snowline_mean_m = dplyr::first(snowline_mean_m),
      area_total = sum(area_km2, na.rm = TRUE),
      area_above_snowline = sum(area_km2[elev_bin >= snowline_mean_m], na.rm = TRUE),
      aar_proxy = area_above_snowline / area_total,
      .groups = "drop"
    )

  print(aar)
}
```

    # A tibble: 2 × 5
      year_group snowline_mean_m area_total area_above_snowline aar_proxy
      <fct>                <dbl>      <dbl>               <dbl>     <dbl>
    1 pre-2006             1983.       59.7                41.8     0.700
    2 post-2006            2299.       58.5                26.1     0.446

``` r
# Plot:  Hypsometry distribution
blue_post <- "#00B8D9"   
red_pre   <- "#E76F51"

p_hypso <- ggplot(hyp, aes(elev_bin, area_km2, color = year_group)) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 0.9, alpha = 0.85) +
  scale_color_manual(
    values = c("pre-2006" = red_pre, "post-2006" = blue_post),
    breaks = c("pre-2006", "post-2006")
  ) +
  labs(
    title = "Elevation distribution of Snowarea for the Fox Glacier (Aug–May)",
    subtitle = "Mean Snow-area [km2] per elevation; pre-2006 vs post-2006",
    x = "Elevation [m]",
    y = expression(paste("Mean Snowarea per class (band) [km"^2, "]")),
    color = "Period"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

if (!is.null(snowline_groups)) {
  p_hypso <- p_hypso +
    geom_vline(
      data = snowline_groups,
      aes(xintercept = snowline_mean_m, color = year_group),
      linetype = "dashed",
      linewidth = 0.8,
      alpha = 0.9,
      show.legend = FALSE
    )
}

ggsave("figures/plots_hypsometry/fox_hypso_distribution.png", p_hypso, width = 9, height = 5, dpi = 250)
p_hypso
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-12-1.png)

``` r
# Plot: Cumulative area above elevation 

red_pre   <- "#E76F51"
blue_post <- "#00B8D9" 

hyp_cum <- hyp %>%
  dplyr::group_by(year_group) %>%
  dplyr::arrange(elev_bin) %>%
  dplyr::mutate(cum_above_km2 = rev(cumsum(rev(area_km2)))) %>%
  dplyr::ungroup()

p_cum <- ggplot(hyp_cum, aes(elev_bin, cum_above_km2, color = year_group)) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(
    values = c("pre-2006" = red_pre, "post-2006" = blue_post),
    breaks = c("pre-2006", "post-2006")
  ) +
  labs(
    title = "Cumulative Snowarea above elevation for Fox Glacier",
    subtitle = "Cumulative mean area remaining above a given elevation",
    x = "Elevation [m]", y = expression(paste("Cumulative area above elevation [km"^2, "]")),
    color = "Period"
  )+
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

if (!is.null(snowline_groups)) {

  # 1) Snowline vertical lines
  p_cum <- p_cum +
    geom_vline(
      data = snowline_groups,
      aes(xintercept = snowline_mean_m, color = year_group),
      linetype = "dashed", linewidth = 0.85, alpha = 0.9,
      show.legend = FALSE
    )
  # 2) Compute AAR points from hyp_cum at the snowline elevation
  aar_pts <- hyp_cum %>%
    dplyr::group_by(year_group) %>%
    dplyr::summarise(area_total = max(cum_above_km2, na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(
      dplyr::select(snowline_groups, year_group, snowline_mean_m),
      by = "year_group"
    ) %>%
    dplyr::left_join(hyp_cum, by = "year_group") %>%
    dplyr::group_by(year_group, area_total, snowline_mean_m) %>%
    dplyr::summarise(
      # nearest cum-value at snowline elevation
      area_above_snowline = cum_above_km2[which.min(abs(elev_bin - snowline_mean_m))],
      elev_at = elev_bin[which.min(abs(elev_bin - snowline_mean_m))],
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      aar_proxy = area_above_snowline / area_total,
      label = paste0(
        "AAR ≈ ", scales::percent(aar_proxy, accuracy = 1), "\n",
        "Area above snowline: ", round(area_above_snowline, 1), " km²"
      )
    )
  # 3) Add point + label
  p_cum <- p_cum +
    geom_point(
      data = aar_pts,
      aes(x = elev_at, y = area_above_snowline, color = year_group),
      size = 2.6,
      show.legend = FALSE
    ) +
    geom_label(
      data = aar_pts,
      aes(x = elev_at, y = area_above_snowline, label = label, color = year_group),
      fill = "white",
      alpha = 0.9,
      label.size = 0,
      size = 3.2,
      hjust = -0.05, vjust = -0.6,
      show.legend = FALSE
    )
}

ggsave("figures/plots_hypsometry/fox_hypso_cumulative_above.png", p_cum, width = 9, height = 5, dpi = 250)

# Plot: Difference per class/band (post-2006 minus pre-2006) 
red_pre   <- "#E76F51"
blue_post <- "#00B8D9"  

hyp_wide <- hyp %>%
  dplyr::select(year_group, elev_bin, area_km2) %>%
  tidyr::pivot_wider(names_from = year_group, values_from = area_km2) %>%
  dplyr::mutate(
    diff_km2 = `post-2006` - `pre-2006`,
    change_dir = dplyr::if_else(diff_km2 >= 0, "Gain (post > pre)", "Loss (post < pre)")
  )

p_diff <- ggplot(hyp_wide, aes(elev_bin, diff_km2)) +
  geom_hline(yintercept = 0, linewidth = 0.6, alpha = 0.8) +
  geom_line(color = "black", linewidth = 0.9) +
  geom_point(aes(color = change_dir), size = 1.6, alpha = 0.9) +
  scale_color_manual(
    values = c("Loss (post < pre)" = red_pre, "Gain (post > pre)" = blue_post)
  ) +
  labs(
    title = "Hypsometry change for the Fox Glacier(post-2006 − pre-2006)",
    subtitle = "Positive values = more snowarea in post-2006 within that elevation class/band; negative = loss",
    x = "Elevation [m]",
    y = expression(paste("Change in mean area per class/band [km"^2, "]")),
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave("figures/plots_hypsometry/fox_hypso_diff_post_minus_pre.png", p_diff, width = 9, height = 5, dpi = 250)


# Plot: Accumulation Area Ratio (AAR)
red_pre   <- "#E76F51"
blue_post <- "#00B8D9"

if (!is.null(aar)) {

  aar2 <- aar %>%
    dplyr::mutate(
      period = factor(year_group, levels = c("pre-2006","post-2006"),
                      labels = c("Before 2006", "After 2006")),
      pct = aar_proxy * 100
    )

  delta <- round(aar2$pct[aar2$period=="After 2006"] - aar2$pct[aar2$period=="Before 2006"], 1)

  p_aar_simple <- ggplot(aar2, aes(period, pct, fill = period)) +
    geom_col(width = 0.62, alpha = 0.95) +
    geom_text(aes(label = paste0(round(pct, 0), "%")), vjust = -0.4,
              fontface = "bold", size = 4) +
    scale_fill_manual(values = c("Before 2006" = red_pre, "After 2006" = blue_post)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = "Percentage of snowarea above the snowline",
      subtitle = paste0("After 2006: ", abs(delta), " percentage points ",
                        ifelse(delta < 0, "lower", "higher"),
                        " → less potential accumulation area"),
      x = NULL,
      y = "Area above snowline [%]"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  ggsave("figures/plots_hypsometry/fox_area_above_snowline_share.png", p_aar_simple, width = 7, height = 4.5, dpi = 250)

}

# tables
readr::write_csv(metrics, "fox_hypso_metrics_by_group.csv")
if (!is.null(aar)) readr::write_csv(aar, "fox_hypso_aar_proxy.csv")
if (!is.null(snowline_groups)) readr::write_csv(snowline_groups, "fox_snowline_groups_for_hypso.csv")
```

``` r
p_cum
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-13-1.png)

``` r
p_diff
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-14-1.png)

``` r
p_aar_simple
```

![](glacier_test_files/figure-commonmark/unnamed-chunk-15-1.png)
