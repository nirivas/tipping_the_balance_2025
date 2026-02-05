

# Load Libraries ----------------------------------------------------------



library(tidyverse)
library(glmmTMB)
library(performance)
library(DHARMa)
library(emmeans)
library(multcomp)
library(multcompView)
library(MuMIn)
library(ggeffects)
library(viridis)
library(ggpubr)
library(vegan)
library(cowplot)
library(png)        
library(grid)
library(pairwiseAdonis)




# Load Data ---------------------------------------------------------------
## And data wrangling


df = read_csv("Data/urchinabundance.csv") |> 
  group_by(yearmonth, location, site, plot) |> 
  mutate(inside  = replace_na(inside, 0), outside = replace_na(outside, 0),
         yearmonth = as.character(yearmonth)) |> 
  mutate(abundance = inside+outside) |>
  rename(period = yearmonth) |> 
  ungroup() |> 
  mutate(period = recode(period,  "2021" = "Pre-die-off",
                                  "202206" = "Early-die-off",
                                  "202210" = "Post-die-off",
                                  "202305" = "1-year-post")) |> 
  mutate(period = factor(period, levels = c("Pre-die-off", "Early-die-off", "Post-die-off", "1-year-post")),
         density = abundance/area)

head(df)

community_matrix = read_csv("Data/community.csv") |> 
  mutate(period = factor(period, levels =c("baseline", "dieoff", "after 6mo", "after 12mo")),
         period = recode(period, 
                         "baseline" = "Pre-die-off", 
                         "dieoff" = "Early-die-off",
                         "after 6mo" = "Post-die-off", 
                         "after 12mo" = "1-year-post")) |> 
  dplyr::select(-...1)

pcover = read_csv("Data/pcover.csv") |> 
  mutate(period = factor(period, levels =c("baseline", "dieoff", "after 6mo", "after 12mo")),
         period = recode(period, 
                         "baseline" = "Pre-die-off", 
                         "dieoff" = "Early-die-off",
                         "after 6mo" = "Post-die-off", 
                         "after 12mo" = "1-year-post"),
         location=as.factor(location), 
         site=as.numeric(site), 
         plot=as.numeric(plot))


df_recr = read_csv("Data/recruitment.csv") 



# Urchin Density Analaysis ------------------------------------------------

#### Exploratory Plot
ggplot(df, aes(period, density, fill = location))+
  #geom_boxplot(outliers = F)+
  geom_point(aes(color = location), position = position_dodge(width = .75))+
  geom_boxplot(outliers = F, alpha = 0.5)+
  theme_classic()


## Change in Urchin Density % ----------------------------------------------


urchins_counts_summary = df |> 
  group_by(period) |> 
  summarise(mean_c = mean(density), sd_c = sd(density))

relative_differences <- urchins_counts_summary |>
  pivot_wider(names_from = period, values_from = c(mean_c, sd_c), 
              names_glue = "{.value}_{period}") |>
  mutate(
    relative_difference_early = (`mean_c_Early-die-off` - `mean_c_Pre-die-off`) / `mean_c_Pre-die-off` * 100,
    relative_difference_post = (`mean_c_Post-die-off` - `mean_c_Early-die-off`) / `mean_c_Early-die-off` * 100,
    relative_difference_year = (`mean_c_1-year-post` - `mean_c_Post-die-off`) / `mean_c_Post-die-off` * 100
  ) |>
  dplyr::select(relative_difference_early, `sd_c_Early-die-off`,
                relative_difference_post, `sd_c_Post-die-off`, 
                relative_difference_year, `sd_c_1-year-post`)

print(relative_differences)

vals <- urchins_counts_summary |>
  filter(period %in% c("Pre-die-off", "Early-die-off")) |>
  transmute(tag = if_else(period == "Pre-die-off", "pre", "early"),
            mean_c, sd_c) |>
  pivot_wider(names_from = tag, values_from = c(mean_c, sd_c))

pre_mean   <- 0.0582
pre_sd     <- 0.153
early_mean <- 0.000639
early_sd   <- 0.00213
n <- 95

pre_sem   <- pre_sd / sqrt(n)
early_sem <- early_sd / sqrt(n)

pct_diff <- 100 * ((early_mean / pre_mean) - 1)

se_pct <- 100 * sqrt( (early_sem / early_mean)^2 +
                        (pre_sem   / pre_mean)^2 )

lwr95 <- pct_diff - 1.96 * se_pct
upr95 <- pct_diff + 1.96 * se_pct

pre_sem; early_sem; pct_diff; se_pct; lwr95; upr95

#### 99% decrease in urchins from before (pre) to post die-off



## Urchin Density Models ----------------------------------------------------



### Tweedie Model -----------------------------------------------------------


density_1 = glmmTMB(density ~ location*period + (1|site), 
                    family = tweedie(), data = df)
#Model check

check_model(density_1)

#' DHARMa can also assess over/under dispersion via simulation.
SimulationOutput1 <- simulateResiduals(fittedModel = density_1, 
                                       refit = FALSE, 
                                       seed = 12345)
testDispersion(SimulationOutput1, 
               type = "DHARMa",
               plot = TRUE)



summary(density_1)
car::Anova(density_1)

options(na.action = "na.fail")
dr_invT = dredge(density_1) |> 
  filter(delta < 4)

top = which.min(dr_invT$df)

top_invT = get.models(dr_invT, subset = top)[[1]]

summary(top_invT)
options(na.action = "na.omit")

#Model summary and fit
#Based on the summary, anova test and dredge, the best most parsimonious model 
#includes only period and site


#Final model for counts
tweedie_periods = glmmTMB(density ~ period + (1|site), 
                         family = tweedie(),
                         data = df)

SimulationOutput1 <- simulateResiduals(fittedModel = tweedie_periods, 
                                       refit = FALSE, 
                                       seed = 12345)
testDispersion(SimulationOutput1, 
               type = "DHARMa",
               plot = TRUE)

check_model(tweedie_periods)
car::Anova(tweedie_periods)
summary(tweedie_periods)


### Predictions -------------------------------------------------------------



gg.density_period = ggpredict(tweedie_periods, terms = c("period"), bias_correction = T) |> 
  rename(period = x, density = predicted)

print(as.data.frame(gg.density_period))

#check region comparisons
model_means_tweediedensity = emmeans(object = tweedie_periods,
                               specs = "period")

pairs(model_means_tweediedensity)

contrast(model_means_tweediedensity)

# add letters to each mean
model_means_cld_density = cld(object = model_means_tweediedensity,
                               Letters = letters,
                               alpha = 0.05)
# show output
gg.density_period2 <- merge(gg.density_period, model_means_cld_density, by = c("period"))



### Plots -------------------------------------------------------------------


#plot of fit for location
density_plot = ggplot(gg.density_period2, aes(period, density, colour = period))+
  geom_pointrange(aes(ymax = conf.high, ymin = conf.low), size = 1, linewidth = 2)+
  geom_text(aes(label = .group, y = .1), 
            size = 6,hjust=.65)+
  scale_color_viridis(option = 'turbo', discrete = TRUE, 
                      labels = c(
                        "Pre-die-off" = "Pre-die-off",
                        "Early-die-off" = "Early-die-off",
                        "Post-die-off" = "Post-die-off",
                        "1-year-post" = "1-year-post"
                      ))+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  theme_classic()+
  labs(x = "Periods relative to urchin die-off", y = "Density of urchins (m\u207B\u00B2)")+
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        axis.line = element_line(colour = "black"))
density_plot



img <- readPNG("Figs/Picture1.png")   
image_plot <- ggdraw() + draw_image(img)


# Top row: a) and b)
top_row <- ggarrange(
  density_plot, 
  labels = c("a)"),
  nrow = 1,
  align = "h",
  legend = "none"
)

# Final layout: top row + bottom image
final_plot <- ggarrange(
  top_row, 
  image_plot, 
  labels = c("", "b)"),  
  ncol = 1,
  heights = c(1, 1)  
)


final_plot


## Benthic Cover Analysis --------------------------------------------------

## Data wrangling ###

coral = pcover |> 
  filter(category %in% c("coral"))|> 
  mutate(period = factor(period, levels = c("Pre-die-off", "Early-die-off", "Post-die-off", "1-year-post")))


algae = pcover |> 
  filter(category %in% c("algae")) |> 
  mutate(period = factor(period, levels = c("Pre-die-off", "Early-die-off", "Post-die-off", "1-year-post")))


sponge = pcover |> 
  filter(category %in% c("sponge")) |> 
  mutate(period = factor(period, levels = c("Pre-die-off", "Early-die-off", "Post-die-off", "1-year-post")))


# Univariate Analysis -----------------------------------------------------



## Corals ------------------------------------------------------------------




### Model -------------------------------------------------------------------



coralbeta = glmmTMB(percentcover ~ period + (1|site), 
                     family = beta_family(), data = coral)



check_model(coralbeta)


summary(coralbeta)
performance(coralbeta)
car::Anova(coralbeta)


### Predictions -------------------------------------------------------------



gg.coralperiod <-ggpredict(coralbeta, terms = c("period")) |> 
  rename(period = x, cover = predicted)

model_means_coralperiodbeta = emmeans(coralbeta, specs = ~ period)
emmeans_coral_resp <- emmeans(coralbeta, ~ period, type = "response")
emmeans_coral_resp

pairs(model_means_coralperiodbeta)

model_means_cld_coralperiodbeta = cld(object = model_means_coralperiodbeta,
                                      Letters = letters,
                                      alpha = 0.05)

gg.coralperiod <- merge(gg.coralperiod, model_means_cld_coralperiodbeta, by = c("period"))
print(gg.coralperiod)


### Plot --------------------------------------------------------------------


cor = ggplot(gg.coralperiod, aes(period, cover*100, colour = period))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position = position_dodge(width = 0.5))+
  geom_text(aes(label = .group,y=1.93 
  ), show.legend = F, size = 6,
  position = position_dodge(width = 0.5))+
  labs(x = NULL, y = "% Coral Cover",color= "period")+
  scale_color_viridis_d(
    option = 'turbo'
  )+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),legend.position = "none",
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank())
cor



## Algae -------------------------------------------------------------------

### Model -------------------------------------------------------------------

algaebeta = glmmTMB(percentcover ~ period + (1|site), 
                    family = beta_family(), data = algae)


check_model(algaebeta)
summary(algaebeta)
performance(algaebeta)

car::Anova(algaebeta)

### Predictions -------------------------------------------------------------

gg.algaeperiod = ggpredict(algaebeta, terms = c("period")) |> 
  rename(period = x, cover = predicted)
print(as.data.frame(gg.algaeperiod))
model_means_algaeperiodbeta = emmeans(algaebeta, specs = ~ period)
pairs(model_means_algaeperiodbeta)
emmeans_algae_resp <- emmeans(algaebeta, ~ period, type = "response")
emmeans_algae_resp

model_means_cld_algaeperiodbeta = cld(object = model_means_algaeperiodbeta,
                                      Letters = letters,
                                      alpha = 0.05)

gg.algaeperiod <- merge(gg.algaeperiod, model_means_cld_algaeperiodbeta, by = "period")
print(gg.algaeperiod)



### Plot --------------------------------------------------------------------

alg <- ggplot(gg.algaeperiod, aes(period, cover*100, colour = period))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position =position_dodge(width = 0.5))+
  labs(x = NULL, y = "% Algae Cover")+
  scale_color_viridis_d(
    option = 'turbo'
  )+
  geom_text(aes(label=.group,  y=40,
  ), size=6,position = position_dodge(width = 0.5))+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
   theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"),  
        legend.text = element_text(size = 12),
        ,
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank())

alg


## Sponge ------------------------------------------------------------------

### Model -------------------------------------------------------------------

spongebeta = glmmTMB(percentcover ~ period + (1|site), 
                     family = beta_family(), data = sponge)


check_model(spongebeta)
summary(spongebeta)
car::Anova(spongebeta)

model_means_spongebeta = emmeans(spongebeta, specs = ~ period)
pairs(model_means_spongebeta)
emmeans_sponge_resp <- emmeans(spongebeta, ~ period, type = "response")
emmeans_sponge_resp


### Predictions -------------------------------------------------------------

gg.spongebeta = ggpredict(spongebeta, terms = c("period")) |> 
  rename(period = x, cover = predicted)
print(as.data.frame(gg.spongebeta))

model_means_cld_spongebeta = cld(object = model_means_spongebeta,
                                 Letters = letters,
                                 alpha = 0.05)

gg.spongebeta <- merge(gg.spongebeta, model_means_cld_spongebeta, by = c("period"))
print(gg.spongebeta)

### Plot --------------------------------------------------------------------

spo = ggplot(gg.spongebeta, aes(period, cover*100, colour = period))+
  geom_pointrange(aes(ymax = conf.high*100, ymin = conf.low*100), size = 1, linewidth = 2,
                  position =position_dodge(width = 0.5))+
  geom_text(aes(label=.group, y = 9),
            size=6)+
  labs(x = NULL, y = "% Sponge Cover")+
  scale_color_viridis_d(
    option = 'turbo')+
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 14, face = "bold"), 
        legend.text = element_text(size = 12),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank()) 

spo


#### Arranged plot -----------------------------------------------------------



combined_plot <- ggarrange(cor, alg, spo,
                           labels = c("a)", "b)", "c)"),
                           ncol = 3, nrow = 1, 
                           common.legend = TRUE,
                           legend = "bottom"
)

combined_plot





# Multivariate Analysis ---------------------------------------------------


## NMDS --------------------------------------------------------------------



nmds <- readRDS("rds/nmds_results.rds")

nmds_coords <- as.data.frame(nmds$points) 

colnames(nmds_coords) <- c("NMDS1", "NMDS2")  

nmds_coords$period <- community_matrix$period   

stress_value <- nmds$stress   
nmds_coords$period <- factor(nmds_coords$period, levels = c("Pre-die-off", "Early-die-off", "Post-die-off", "1-year-post"))

hulls <- nmds_coords %>%
  group_by(period) %>%
  slice(chull(NMDS1, NMDS2))

nmdscom = ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = period)) +  
  geom_point(size = 3, alpha = 0.7) +  
  labs(title = "NMDS of Community Data", x = "NMDS 1", y = "NMDS 2", color = 
         "Periods") +
  annotate("text", x = Inf, y = Inf, hjust = 1.0, vjust = 1.1,
           label = paste("Stress =", round(stress_value, 4)), 
           size = 5, fontface = "bold") + 
  theme_classic() +   
  theme(axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_viridis_d(option = "turbo"
  )

nmdscom

custom_colors <- c(
  "Pre-die-off"     = "#0072B2",  # blue
  "Early-die-off"  = "#E69F00",  # orange
  "Post-die-off"    = "#009E73",  # green
  "1-year-post"  = "#CC79A7"   # pink/purple
)


nmdscom <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2)) +  
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = period, group = period),
               alpha = 0.2, color = NA) +  
  geom_point(size = 3, alpha = 0.7, aes(color = period)) +  
  labs(title = "NMDS of Community Data", x = "NMDS 1", y = "NMDS 2", color = "Periods") +
  annotate("text", x = Inf, y = Inf, hjust = 1.0, vjust = 1.1,
           label = paste("Stress =", round(stress_value, 4)), 
           size = 5, fontface = "bold") + 
  geom_segment(data = vec, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  ggrepel::geom_text_repel(data = vec, aes(x = NMDS1, y = NMDS2, label = spp),
                           size = 4) +
  theme_classic() +   
  theme(axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors, guide = "none")  

nmdscom


# multivariate----

# m_b = df_b |> 
#   
#   filter(impact == 'Inside') |> 
#   
#   mutate(id = paste(station, year, sep = '_')) |> 
#   
#   select(id, `Mojarra spp`:`Spanish sardine`) |> 
#   
#   column_to_rownames(var = 'id')
# 
# vars = df_b |> 
#   
#   filter(impact == 'Inside') |> 
#   
#   select(station, year, dieoff, impact)
# 
# # site = tibble(station = c(146, 168, 241, 266))
# 
# pb = adonis2(m_b ~ dieoff, data = vars)
# 
# nmds_b = metaMDS(m_b, distance = "bray", k = 2, try = 100)

# vectors 

m_b = community_matrix |> 
  dplyr::select(-1,-2,-3)

sp.fit = envfit(nmds, m_b, permutations = 999)

vec = tibble(spp = rownames(scores(sp.fit, display = "vectors")),
             
             NMDS1 = scores(sp.fit, display = "vectors")[,1],
             
             NMDS2 = scores(sp.fit, display = "vectors")[,2],
             
              p = sp.fit$vectors$pvals) 
   
  # 
  # filter(p < 0.05)



# df = tibble(vars, data.frame(nmds_b[['points']]))
# 
# hull = df  |> 
#   
#   group_by(dieoff) |> 
#   
#   slice(chull(MDS1, MDS2))
# 
# # nmds plot
# 
# ggplot(df, aes(MDS1, MDS2))+
#
#   geom_polygon(data = hull, aes(color = dieoff, fill = dieoff), alpha = 0.5) +
#
#   geom_point(aes(shape = as.factor(station),
#
#                  color = dieoff, fill = dieoff), size = 3)+
#
#   labs(x = 'NMDS1', y = 'NMDS2',
#
#        color = 'Die-off', fill = 'Die-off', shape = 'Station')+
#
#   scale_color_manual(values = c('darkgreen', 'burlywood2'))+
#
#   scale_fill_manual(values = c('darkgreen', 'burlywood2'))+
#
#   geom_segment(data = vec,
#
#                aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
#
#                arrow = arrow(length = unit(0.25, "cm")), color = "black") +
#
#   ggrepel::geom_text_repel(data = vec, aes(x = MDS1, y = MDS2, label = spp),
#
#                            size = 4)+
#
#   theme_bw()+
#
#   theme(axis.title = element_text(size = 14),
#
#         axis.text.y = element_text(size = 14, colour = "black"),
#
#         axis.text.x = element_text(size = 12, colour = "black"),
#
#         plot.title = element_text(size = 14, hjust=0.5),
#
#         panel.grid.major = element_blank(),
#
#         panel.grid.minor = element_blank(),
#
#         legend.position = 'right',
#
#         legend.title = element_text(size = 14),
#
#         strip.text.x = element_text(size = 14),
#
#         legend.text = element_text(size = 12))











## PERMANOVA ---------------------------------------------------------------



community_matrix <- pcover |>
  mutate(plot = as.numeric(plot), site = as.numeric(site),location = as.factor(location)) |>
  filter(category %in% c("coral", "sponge", "algae","seagrass")) |>
  dplyr::select(site, plot, period, category, percentcover) |>
  pivot_wider(
    names_from = category,
    values_from = percentcover,
    values_fill = 0
  )

distance_matrix <- vegdist(community_matrix[4:7], method = "bray") 


permanova_results <- adonis2(distance_matrix ~ period, 
                             data = community_matrix, 
                             permutations = 999,
                             strata = community_matrix$site, 
                             by ="terms"
) 

permanova_results
summary(permanova_results)

#Test homogeneity of multivariate dispersions
dispersion <- betadisper(distance_matrix, group = community_matrix$period)

# Permutation test for differences in dispersion
permdisp_results <- permutest(dispersion, permutations = 999)

permdisp_results

# Perform pairwise comparisons for 'period'
pairwise_results_period <- pairwiseAdonis::pairwise.adonis2(distance_matrix ~ period, data = community_matrix, nperm = 999, p.method = "bonferroni")

pairwise_results_period




# Log Resonse Ratios ------------------------------------------------------



## Urchins and Algae  -------------------------------------------------


### Data wrangle ------------------------------------------------------------



cover2<-pcover |> 
  dplyr::select( period, category, percentcover, site, plot) |>
  pivot_wider(names_from = category, values_from = percentcover)

algae2 = algae |> 
  dplyr::select(-yearmonth,-category) |> 
  mutate(plot = as.character(plot), site= as.character(site),
         plot = sub("^0+", "", plot),  
         site = sub("^0+", "", site),
         location = as.factor(location),
         plot = as.factor(plot),
         site = as.factor(site))

anyNA(algae2)
colSums(is.na(algae2))

urchins_counts2 = df |> 
  mutate(plot = as.character(plot),site=as.character(site),
         plot = sub("^0+", "", plot),  
         site = sub("^0+", "", site),
         plot = as.factor(plot),
         site = as.factor(site)) |> 
  pivot_wider(names_from = period, values_from = c(density))

head(urchins_counts2)


log_algae = algae2 |> 
  pivot_wider(names_from = period, values_from = percentcover) |> 
  mutate(log_ratio = log(`Post-die-off`/`Pre-die-off`)) |> 
  drop_na() |> 
  dplyr::select(-`Pre-die-off`, -`Early-die-off`, -`Post-die-off`, -`1-year-post`)

head(log_algae)

algaeratiodf = left_join(log_algae, urchins_counts2, by = c("plot", "site", "location")) |> 
  filter(!is.na(`Pre-die-off`))|> 
  rename(pre = `Pre-die-off`)


### Model -------------------------------------------------------------------

lrralgae = glmmTMB(log_ratio ~  pre * location + (1|site),
             family = gaussian, data = algaeratiodf)


options(na.action = "na.fail")
dr_invT = dredge(lrralgae) |> 
  filter(delta < 4)

top = which.min(dr_invT$df)

top_invT = get.models(dr_invT, subset = top)[[1]]

summary(top_invT)
options(na.action = "na.omit")

#Top model

algaeratiomodel = glmmTMB(log_ratio ~ pre + (1|site),
           family = gaussian, data = algaeratiodf)

check_model(algaeratiomodel)
car::Anova(algaeratiomodel)
summary(algaeratiomodel)

### Predictions -------------------------------------------------------------

gg.algaeratiomodel <- ggpredict(algaeratiomodel, terms = "pre [all]") |>
  rename(Density = x, cover = predicted)



print(as.data.frame(gg.algaeratiomodel))

### Plot --------------------------------------------------------------------

ratiocounts = ggplot(gg.algaeratiomodel, aes(Density, cover)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#2A788EFF", alpha = 0.2) +
  scale_fill_viridis_c(option = "D")+
  labs(x = "Urchin Density", y = "Log Response Ratio of Algae Cover \n (Post-die-off / Pre-die-off)") +
  theme_classic() +
  theme(axis.text = element_text(size = 16, face = "bold", colour = "black"),
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        plot.title = element_text(size = 18, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) 

ratiocounts





## Recruitment  ------------------------------------------------------------

### Data wrangle ------------------------------------------------------------

log_recr <- df_recr |> 
  pivot_wider(names_from = period, values_from = recruits) |> 
  mutate(`Summer 2023` = replace_na(`Summer 2023`, 0)) |> 
  mutate(log_ratio = log(`Summer 2023` / `Baseline`)) |> 
  filter(is.finite(log_ratio)) |>  # This removes -Inf, Inf, and NaN
  dplyr::select(-Baseline, -`Summer 2023`) |> 
  mutate(plot = as.factor(plot),
         site = as.factor(site))

ratiodf_recr = left_join(log_recr, urchins_counts2, by = c("plot", "site", "location")) |> 
  filter(!is.na(`Pre-die-off`)) |> 
  mutate(pre = `Pre-die-off`)


head(ratiodf_recr)

### Model -------------------------------------------------------------------

rcr_ratio = glmmTMB(log_ratio ~  pre * location + (1|site),
             family = gaussian, data = ratiodf_recr)


options(na.action = "na.fail")

dr_invT = dredge(rcr_ratio) |> 
  filter(delta < 4)

top = which.min(dr_invT$df)

top_invT = get.models(dr_invT, subset = top)[[1]]

summary(top_invT)
options(na.action = "na.omit")



rcrpre <- glmmTMB(log_ratio ~ pre + (1|site),
                   family = gaussian, data = ratiodf_recr)


rcrloc <- glmmTMB(log_ratio ~ location + (1|site),
                   family = gaussian, data = ratiodf_recr)

compare_performance(rcrpre, rcrloc)

check_model(rcrpre)
summary(rcrpre)
car::Anova(rcrpre)

check_model(rcrloc)
summary(rcrloc)
car::Anova(rcrloc)

##pre values are the better model

### Predictions -------------------------------------------------------------

gg.ratiomodel = ggpredict(rcrpre, terms = c("pre [all]")) |> 
  rename(counts = x, recruits = predicted)


print(as.data.frame(gg.ratiomodel))

### Plot --------------------------------------------------------------------

ratiorecru = ggplot(gg.ratiomodel, aes(counts, recruits)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#2A788EFF", alpha = 0.2) +
  scale_fill_viridis_c(option = "D")+
  labs(x = "Urchin Density", y = "Log Response Ratio of Recruitment \n (1-year-post / Pre-die-off)") +
  theme_classic() +
  theme(axis.text = element_text(size = 16, face = "bold", colour = "black"),
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        plot.title = element_text(size = 18, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  geom_text(x = 24.5, 
            y = 2.1 , 
            label = "p = 0.0045", 
            size = 6, color = "black")

ratiorecru





# Sponges and Urchins -----------------------------------------------------


### Data wrangle ------------------------------------------------------------

sponge2 = sponge |> 
  dplyr::select(-yearmonth,-category) |> 
  mutate(plot = as.character(plot), site= as.character(site),
         plot = sub("^0+", "", plot),  # Remove leading zeros from 'plot'
         site = sub("^0+", "", site),
         location = as.factor(location),
         plot = as.factor(plot),
         site = as.factor(site))

anyNA(sponge2)
colSums(is.na(sponge2))


head(urchins_counts2)


log_sponge = sponge2 |> 
  pivot_wider(names_from = period, values_from = percentcover) |> 
  mutate(log_ratio = log(`Post-die-off`/`Pre-die-off`)) |> 
  drop_na() |> 
  dplyr::select(-`Pre-die-off`, -`Early-die-off`, -`Post-die-off`, -`1-year-post`)

head(log_sponge)

ratiodf2 = left_join(log_sponge, urchins_counts2, by = c("plot", "site", "location")) |> 
  filter(!is.na(`Pre-die-off`)) |> 
  rename(pre = `Pre-die-off`) |> 
  mutate(location = as.factor(location))

head(ratiodf2)


### Model -------------------------------------------------------------------

ms = glmmTMB(log_ratio ~  pre  + (1|site),
             family = gaussian, data = ratiodf2)


options(na.action = "na.fail")
dr_invT = dredge(ms) |>
  filter(delta < 4)
top = which.min(dr_invT$df)
top_invT = get.models(dr_invT, subset = top)[[1]]
summary(top_invT)
options(na.action = "na.omit")

ms = glmmTMB(log_ratio ~  location + (1|site),
             family = gaussian, data = ratiodf2)

check_model(ms)
car::Anova(ms)
summary(ms)

### Predictions -------------------------------------------------------------

gg.ratiomodel4 = ggpredict(ms, terms = c("location")) |> 
  rename(location = x, cover = predicted)

plot(ggpredict(ms, terms = c("location ")))


# Correlation analysis ----------------------------------------------------



## Urchin metrics correlation ----------------------------------------------



urchins = read_csv("./data/urchin_counts.csv") |>   mutate(period = recode(period, 
                                                                           "baseline" = "Pre-die-off", 
                                                                           "dieoff" = "Early-die-off",
                                                                           "after 6mo" = "Post-die-off", 
                                                                           "after 12mo" = "1-year-post"),
                                                           location=as.factor(location), 
                                                           site=as.numeric(site), 
                                                           plot=as.numeric(plot))

urchins

df2 = left_join(df, urchins, by = c("period", "location", "site", "plot"))

correlation_results2 <- df2 %>%
  mutate(
    cluster_abundance_corr = cor(density, counts, use = "complete.obs"),
    cluster_p = cor.test(density, counts)$p.value,
    
    size_abundance_corr = cor(density, mean_size, use = "complete.obs"),
    size_p = cor.test(density, mean_size)$p.value,
    
    density_abundance_corr = cor(density, abundance, use = "complete.obs"),
    density_p = cor.test(density, abundance)$p.value,
    
  )

AC = ggplot(correlation_results2, aes(density,counts))+
  geom_point(aes(color = period))+
  geom_smooth(method="lm")+
  scale_color_viridis(option = 'turbo', discrete = T)+
  #facet_wrap(~location)+
  theme_classic()+
  labs(color = "Periods", y = "Count of Urchin Clusters", x= "Urchin Density")+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16, face = "bold"),  # Legend title size
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 14, face = "bold")
  ) + 
  geom_text(
    aes(
      x = Inf, y = Inf,  # Position in the top right corner of each facet
      label = paste("r =", round(cluster_abundance_corr, 2), "\np =", format.pval(cluster_p, digits = 3, eps = 0.001))
    ),
    hjust = 1.2, vjust = 1.2, size = 6, color = "black",
    data = correlation_results2
  )

AC

AM = ggplot(correlation_results2, aes(density,mean_size))+
  geom_point(aes(color = period))+
  geom_smooth(method="lm")+
  scale_color_viridis(option = 'turbo', discrete = T)+
  #facet_wrap(~location)+
  theme_classic()+
  labs(color = "Periods", y = "Cluster Size", x= "Urchin Density")+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16, face = "bold"),  # Legend title size
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 14, face = "bold")
  ) + 
  geom_text(
    aes(
      x = Inf, y = Inf,  # Position in the top right corner of each facet
      label = paste("r =", round(size_abundance_corr, 2), "\np =", format.pval(size_p, digits = 3, eps = 0.001))
    ),
    hjust = 1.2, vjust = 1.2, size = 6, color = "black",
    data = correlation_results2
  )

AM

DA = ggplot(correlation_results2, aes(density,abundance))+
  geom_point(aes(color = period))+
  geom_smooth(method="lm")+
  scale_color_viridis(option = 'turbo', discrete = T)+
  #facet_wrap(~location)+
  theme_classic()+
  labs(color = "Periods", y = "Urchin Abundance", x= "Urchin Density")+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16, face = "bold"),  # Legend title size
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold", colour = "black"),
        plot.title = element_text(size = 16, face = "bold", colour = "black"),
        strip.text = element_text(size = 14, face = "bold")
  ) + 
  geom_text(
    aes(
      x = Inf, y = Inf,  # Position in the top right corner of each facet
      label = paste("r =", round(density_abundance_corr, 2), "\np =", format.pval(density_p, digits = 3, eps = 0.001))
    ),
    hjust = 1.2, vjust = 1.2, size = 6, color = "black",
    data = correlation_results2
  )

DA

ggarrange(AC,AM,DA,
          nrow=1, ncol=3,common.legend = TRUE, # To have a common legend for all plots
          legend = "bottom",
          labels = c("a)", "b)", "c)"))


