setwd("~/Project/Palaeofire_v4")
path<- "/Users/yichengshen/Project/Palaeofire_v4"
# create a folder to save figures
dir.create(paste0(path, "/figures"))
# create a folder to save some tables in the process of data analysis
dir.create(paste0(path, "/process0311"))

library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(rgdal)
library(raster)
library(sf)
library(latticeExtra)
library(sp)
library(maptools)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
#install.packages('ncdf4',repos='http://cran.us.r-project.org')
library(ncdf4)
library(tmap)
library(locfit)

# Entity list (entities with both pollen and charcoal) -- entitylist_pc.xlsx
# the "used" column indicate whether it has modern charcoal sample (age<=100).
# the "rep" column indicate whether it has micro and macro entities
# we do not consider entities without macro_micro, used and rep information (for some reasons)

entitylist<- read_excel("~/Project/Palaeofire_v4/entitylist_pc.xlsx") %>%
  filter(is.na(used)==F)

mc<- entitylist %>% dplyr::filter(used==1) 

# Load pollen data
# p-- all pollen data
# p_relationship-- pollen data used in deriving fire-veg relationships

p<- read.csv("/Users/yichengshen/Project/Palaeofire_v4/Iberia_data_v3/Iberia_pollen_records_v3_0307.csv")
p<- p %>% 
  dplyr::rename("depth"="avg_depth..cm.", 
                        "median"="INTCAL2020_median", 
                        "Entity.name"="entity_name",
                        "Site.Name"="site_name") 

p$IPE.age..cal.<- as.numeric(p$IPE.age..cal.)
p$median<- as.numeric(p$median)

# Villarquemado use IPE ages
for(i in 1:length(p$Site.Name)){
  if (is.na(p$median[i])==TRUE){
    p$median[i]<- p$IPE.age..cal.[i]
  }
}

## target study region-- iberia peninsula
ip_shp_rgdal<- readOGR("~/Project/Palaeofire/UseForMap/iberianpeninsula/iberianpeninsula_doc.shp")
ip<- st_as_sf(ip_shp_rgdal)
st_crs(ip)<- "+proj=longlat"

## filter the entities outsite the main region of IP

p.entity <- p %>% dplyr::select(Site.Name, Entity.name, latitude, longitude) %>%
  distinct() %>%
  st_as_sf(., coords = c("longitude", "latitude"))

st_crs(p.entity)<- "+proj=longlat"
p.entity$testlocation<- st_within(p.entity$geometry, ip, sparse = FALSE)
p.entity<- p.entity %>% filter(testlocation==T) %>% dplyr::select(-testlocation) %>%
  st_drop_geometry(.)
# 122 to 113 entities

p<- p %>% filter(Entity.name %in% p.entity$Entity.name)

p_relationship<- p %>% dplyr::filter(Entity.name %in% mc$pollen_entity)
# 33 entities

# Pollen bin ####
p_relationship_bin<- pollen_bin_function(p_relationship)
# write.csv(p_relationship_bin, "process0311/p_relationship_bin.csv", row.names = F)



# Load charcoal data
c<- read.csv("/Users/yichengshen/Project/Palaeofire_v4/Iberia_data_v3/Iberia_charcoal_records_v3_0307.csv")
c<- c %>% dplyr::rename("depth"="avg_depth..cm.", 
                        "median"="INTCAL2020_median",
                        "quantity"="charcoal.quantity")
c<- c %>% filter(entity_name %in% entitylist$charcoal_entity)

## update on 2022-01-25
#pollen concentration-- C0P0
#concentration-- CONC 
#raw count-- COUN
#count-- COUN
#influx-- INFL
#other-- OTHE
#per unit weight-- per unit weight
for (i in 1: length(c$TYPE)){
  if (c$TYPE[i]=="concentration"){
    c$TYPE[i]<- "CONC"
  }
  if (c$TYPE[i]=="pollen concentration"){
    c$TYPE[i]<- "C0P0"
  }
  if (c$TYPE[i]=="count"){
    c$TYPE[i]<- "COUN"
  }
  if (c$TYPE[i]=="raw count"){
    c$TYPE[i]<- "COUN"
  }
  if (c$TYPE[i]=="influx"){
    c$TYPE[i]<- "INFL"
  }
  if (c$TYPE[i]=="other"){
    c$TYPE[i]<- "OTHE"
  }
}

# Add ID_SAMPLE, ID_SITE, ID_ENTITY
c <- c %>% mutate(ID_SAMPLE = row_number())
c$ID_ENTITY <- c %>% dplyr::group_indices(entity_name)
c$ID_SITE <- c %>% dplyr::group_indices(site_name)

# Charcoal influx transformation ####
c_influx<- charcoal_fluxtrans_function(c)

# Charcoal bin ####
# c_influx<- read.csv("charcoaltrans/checkdat.csv")
c_influx<- c %>% 
  dplyr::select(-ID_ENTITY, -ID_SITE, -depth) %>%
  dplyr::left_join(c_influx, by="ID_SAMPLE")

c_influx<- c_influx %>% filter(sed_rate>=0) 
charcoal_bin<- charcoal_bin_function(c_influx)
# write.csv(charcoal_bin, "process0311/charcoal_bin.csv", row.names = F)

# Conversion factor####
# Filter entities with modern charcoal
modern_entities<- charcoal_bin%>%
  dplyr::filter(age==100) 

# Extract entities pre_ba
ip_shp_rgdal<- readOGR("/Users/yichengshen/Project/Palaeofire_v2/UseForMap/iberianpeninsula/iberianpeninsula_doc.shp")
ba<- brick("/Users/yichengshen/Project/Palaeofire_v2/UseForMap/2001_2016_bapre_mean.nc")
modern_entities_points<- st_as_sf(modern_entities, coords = c("longitude", "latitude"), crs=4326)
bavalue<- raster::extract(ba, modern_entities_points)
modern_entities$ba_pre<- bavalue*12 # calculate annual mean
# add macro_micro info
entitylistsmodern<- entitylist %>% dplyr::select(charcoal_entity, macro_micro, rep) %>% dplyr::rename(entity_name=charcoal_entity)
modern_entities<- modern_entities %>% dplyr::left_join(entitylistsmodern, by="entity_name")
# write.csv(modern_entities, "process0311/modern_entities.csv", row.names = F, fileEncoding = "UTF-8")

# Conversion factor conversion
charcoalentityhasmodern<-data.frame()
for (i in 1:length(modern_entities$entity_name)){
  modern_target<- modern_entities$entity_name[i]
  targetentity<- charcoal_bin %>% 
    filter(entity_name %in% modern_target)
  factor<- modern_entities$ba_pre[i]/modern_entities$max[i] #try square-root (or log) ba_pre#not good
  modern_entities$factor[i]<- factor
  targetentity$weighteddata<- targetentity$max * factor
  charcoalentityhasmodern<- rbind(charcoalentityhasmodern, targetentity)
}

# Ignore entity that weighting factor = inf
inf_entity<- modern_entities$entity_name[which(modern_entities$factor==Inf)]
charcoalentityhasmodern<- charcoalentityhasmodern %>% dplyr::filter(entity_name %notin% inf_entity)
## 41entities

# Merge charcoal and pollen ####
pollen_bin<- p_relationship_bin
entities<- unique(pollen_bin$Entity.name)
df_merge<- data.frame()
for (entity in entities) {
  # entity<- entities[2] 
  targetsite_p<- pollen_bin %>% filter(Entity.name==entity)
  c_name_id<- which(entitylist$pollen_entity==entity)
  if (length(c_name_id)==2){
    c_name<- entitylist %>% dplyr::filter(pollen_entity == entity) %>% dplyr::filter(macro_micro=="macro")
  }else{
    c_name<- entitylist %>% dplyr::filter(pollen_entity == entity)
  }
  c_name<- c_name$charcoal_entity
  targetsite_c<- charcoalentityhasmodern %>% dplyr::filter(entity_name %in% c_name)
  targetmerge<- inner_join(targetsite_c, targetsite_p, by = "age")
  df_merge<- rbind(df_merge, targetmerge)
}
pc_merge_modern<- df_merge
# 31 entities
# write.csv(pc_merge_modern, "process0311/pc_merge_modern.csv", row.names = F, fileEncoding = "UTF-8")


# Check rare taxa
startafter<- which(colnames(pc_merge_modern) == "reference")
df_taxa<- pc_merge_modern[,(startafter+1):length(colnames(pc_merge_modern))] 
countrare<- as.data.frame(colSums(df_taxa != 0))
colnames(countrare)<- "sum0"
raretaxa<- countrare %>% dplyr::filter(sum0<=5)

# Remove rare taxa (occurence <=5) ####
taxa_start<- which(colnames(p) == "Abies")
taxa_end<- length(colnames(p))
p_taxa<- p[,taxa_start:taxa_end] #205
p_info<- p[,1:taxa_start-1]
p_taxa<- dplyr::select(p_taxa, -(rownames(raretaxa)))
# p_relationship<- p %>% dplyr::filter(Entity.name %in% mc$pollen_entity)
p_rmrare<- cbind(p_info, p_taxa)
p_rmrare_bin<- pollen_bin_function(p_rmrare)
# write.csv(p_rmrare_bin, "process0311/p_rmrare_bin.csv", row.names = F, fileEncoding = "UTF-8")

# Merge pollen and charcoal after remove rare taxa ####
pollen_bin<- p_rmrare_bin
entities<- unique(pollen_bin$Entity.name)
df_merge<- data.frame()
for (entity in entities) {
  # entity<- entities[6] 
  targetsite_p<- pollen_bin %>% filter(Entity.name==entity)
  c_name_id<- which(entitylist$pollen_entity==entity)
  if (length(c_name_id)==2){
    c_name<- entitylist %>% dplyr::filter(pollen_entity == entity) %>% dplyr::filter(macro_micro=="macro")
  }else{
    c_name<- entitylist %>% dplyr::filter(pollen_entity == entity)
  }
  if (length(c_name$charcoal_entity)>0){
    c_name<- c_name$charcoal_entity
    targetsite_c<- charcoalentityhasmodern %>% dplyr::filter(entity_name == c_name)
    targetmerge<- inner_join(targetsite_c, targetsite_p, by = "age")
    df_merge<- rbind(df_merge, targetmerge)
  }
}
# write.csv(df_merge, "process0311/pc_merge_modern_rm.csv", row.names = F, fileEncoding = "UTF-8")
pc_merge_modern<- df_merge

# fxTWAPLS ####
# Use Box-Cox transformation to reduce skewness
BCTransform <- function(y, lambda=0) {
  if (lambda == 0L) { log(y) }
  else { (y^lambda - 1) / lambda }
}
BCTransformInverse <- function(yt, lambda=0) {
  if (lambda == 0L) { exp(yt) }
  else { exp(log(1 + lambda * yt)/lambda) }
}

pc_merge_modern$BC<- BCTransform(pc_merge_modern$weighteddata, lambda=0.25)
w_charcoal<- pc_merge_modern$BC
ggplot(pc_merge_modern, aes(x=BC))+
  geom_histogram(binwidth=0.1, fill="grey", color="grey", alpha=0.9)+
  xlab("Burnt area fraction (Box-Cox transformation)")+
  theme_bw()
ggsave("figures/baf_bc_distribution.pdf", width = 4, height = 4)

ggplot(pc_merge_modern, aes(x=weighteddata))+
  geom_histogram(binwidth=0.01, fill="grey", color="grey", alpha=0.9)+
  xlab("Burnt area fraction")+
  theme_bw()
ggsave("figures/baf_distribution.pdf", width = 4, height = 4)

start<- which(colnames(pc_merge_modern) == "Abies")
end<- which(colnames(pc_merge_modern) == "Ziziphus")
df_taxa<- pc_merge_modern[,(start):(end)] 
df_taxa<- df_taxa[, colSums(df_taxa != 0) > 0] ## Just to make sure we removed all the taxa that do not appear

library(fxTWAPLS)
fx_charcoal <- fxTWAPLS::fx(w_charcoal, bin = 0.002, show_plot = TRUE)
#wapls
fit_charcoal <- fxTWAPLS::WAPLS.w(modern_taxa = df_taxa, modern_climate = w_charcoal, nPLS = 8, usefx = FALSE, fx = NA)
#twapls
fit_t_charcoal <- fxTWAPLS::TWAPLS.w(modern_taxa = df_taxa, modern_climate = w_charcoal, nPLS = 8,usefx = FALSE, fx = NA)
#fxwapls
fit_f_charcoal <- fxTWAPLS::WAPLS.w(modern_taxa = df_taxa, modern_climate = w_charcoal, nPLS = 8, usefx = TRUE, fx = fx_charcoal)
#fxtwapls
fit_tf_charcoal <- fxTWAPLS::TWAPLS.w(modern_taxa = df_taxa, modern_climate = w_charcoal, nPLS = 8, usefx = TRUE, fx = fx_charcoal)

cv_charcoal<-fxTWAPLS::cv.w(df_taxa,w_charcoal,nPLS=8, fxTWAPLS::WAPLS.w, fxTWAPLS::WAPLS.predict.w) #, cpus = 2, test_mode = test_mode
cv_t_charcoal<-fxTWAPLS::cv.w(df_taxa,w_charcoal,nPLS=8, fxTWAPLS::TWAPLS.w, fxTWAPLS::TWAPLS.predict.w)
cv_f_charcoal<-fxTWAPLS::cv.w(df_taxa,w_charcoal,nPLS=8, fxTWAPLS::WAPLS.w, fxTWAPLS::WAPLS.predict.w,usefx=TRUE,fx=fx_charcoal) #, cpus = 2, test_mode = test_mode
cv_tf_charcoal<-fxTWAPLS::cv.w(df_taxa,w_charcoal,nPLS=8, fxTWAPLS::TWAPLS.w, fxTWAPLS::TWAPLS.predict.w,usefx=TRUE,fx=fx_charcoal)

rand_charcoal<-fxTWAPLS::rand.t.test.w(cv_charcoal,n.perm=999)
rand_t_charcoal<-fxTWAPLS::rand.t.test.w(cv_t_charcoal,n.perm=999)
rand_f_charcoal<-fxTWAPLS::rand.t.test.w(cv_f_charcoal,n.perm=999)
rand_tf_charcoal<-fxTWAPLS::rand.t.test.w(cv_tf_charcoal,n.perm=999)

methods_results<- rbind(rand_charcoal,rand_t_charcoal,rand_f_charcoal,rand_tf_charcoal)
write.csv(methods_results, "methods_results_BC.csv")

# Plot fitted reconstructed burnt area (SI)
f1<- fxTWAPLS::plot_train(fit_charcoal, 4)
f2<- fxTWAPLS::plot_train(fit_t_charcoal, 5) 
f3<- fxTWAPLS::plot_train(fit_f_charcoal, 2) 
f4<- fxTWAPLS::plot_train(fit_tf_charcoal, 4)


library(ggpubr)
ggarrange(f1, f2, f3, f4,
          #labels = c("WAPLS", "TWAPLS", "fxWAPLS", "fxTWAPLS"),
          ncol = 2, nrow = 2)
ggsave("figures/fitted_combined.pdf", width = 8, height = 8)
#size 8*8

# Plot residuals (SI)
r1<- fxTWAPLS::plot_residuals(fit_charcoal, 4)
r2<- fxTWAPLS::plot_residuals(fit_t_charcoal, 5)
r3<- fxTWAPLS::plot_residuals(fit_f_charcoal, 2) 
r4<- fxTWAPLS::plot_residuals(fit_tf_charcoal, 4)

library(ggpubr)
ggarrange(r1, r2, r3, r4,
          # labels = c("WAPLS", "TWAPLS", "fxWAPLS", "fxTWAPLS"),
          ncol = 2, nrow = 2)
ggsave("figures/residual_combined.pdf", width = 8, height = 8)

# Plot fitted and residual plots for fxTWAPLS
# fitted plot
train_output<- fit_tf_charcoal
col<- 4
x <- train_output[["x"]]
fitted <- train_output[["fit"]][, col]
plotdata <- cbind.data.frame(x, fitted)

max <- max(fitted, x)
min <- min(fitted, x)

# plot the fitted curve, the black line is the 1:1 line, the red line is the 
# linear regression line to fitted and x, which shows the overall compression
ggplot2::ggplot(plotdata, ggplot2::aes(x, fitted)) + 
  ggplot2::geom_point(size = 0.4) + ggplot2::theme_bw() +
  ggplot2::geom_abline(slope = 1, intercept = 0) + 
  # ggplot2::xlim(min, 0.35) + ggplot2::ylim(min, 0.35) +
  ggplot2::geom_smooth(method = 'lm',
                       formula = y ~ x,
                       color = 'red')+
  ylab("fxTWA-PLS fitted burnt area fraction") + xlab("Burnt area fraction")
ggsave("figures/fxtwapls_fitted.pdf", width = 5, height = 5)
#size 5*5

# residual plot
x <- train_output[["x"]]
residuals <- train_output[["fit"]][, col] - train_output[["x"]]
plotdata <- cbind.data.frame(x, residuals)

maxr <- max(abs(residuals))

# plot the residuals, the black line is 0 line, the red line is the locally 
# estimated scatterplot smoothing, which shows the degree of local compression
ggplot2::ggplot(plotdata, ggplot2::aes(x, residuals)) + 
  geom_rect(aes(xmin=-3.25,xmax=-2.5,ymin=-1.5,ymax=1.5),alpha=0.5,fill="#D8D8D8")+
  ggplot2::geom_point(size = 0.4) + ggplot2::theme_bw() +
  ggplot2::geom_abline(slope = 0, intercept = 0) + 
  ggplot2::xlim(min(x), max(x)) + ggplot2::ylim(-1.5, 1.5) +
  ggplot2::geom_smooth(method = 'loess',
                       color = 'red', 
                       formula = "y ~ x",
                       se = FALSE)+
  ylab("Residuals of fxTWA-PLS Burnt area fraction") + xlab("Burnt area fraction")
ggsave("figures/fxtwapls_residual.pdf", width = 5, height = 5)

# Plot modern sites and fossil sites ####
p_location<- pollen_bin %>% dplyr::select(Site.Name, Entity.name, latitude, longitude, elevation) %>% 
  dplyr::distinct(Entity.name, .keep_all = T) 

p_location<-  st_as_sf(p_location, coords = c("longitude","latitude"),crs=4326)
mc_location<- charcoalentityhasmodern %>% 
  dplyr::filter(entity_name %in% mc$charcoal_entity) %>%
  dplyr::select(site_name, entity_name, latitude, longitude, elevation) %>% 
  dplyr::distinct(entity_name, .keep_all = T)
mc_location<-  st_as_sf(mc_location, coords = c("longitude","latitude"),crs=4326)

iplow_topo<- raster("~/Project/Palaeofire_v2/UseForMap/iplow_topo.nc")
iplow_topo[iplow_topo < 0] <- 0

map<-tm_shape(iplow_topo) +
  tm_graticules(lines = FALSE,labels.size = 1)+
  tm_raster(title = "", palette=c("#00A600FF", "#24B300FF", "#4CBF00FF", "#7ACC00FF", "#ADD900FF", "#E6E600FF", "#E8C727FF",
                                  "#EAB64EFF", "#ECB176FF"),breaks = c(0, 100,500,1000, 2000, 3500))+
  tm_shape(p_location) +
  tm_symbols(col = "#7491F2", border.col = "white",size = 0.5, shape=21)+
  tm_shape(mc_location) +
  tm_symbols(col = "#F27474", border.col = "white",size = 0.5, shape=25)+
  tm_layout(legend.show=TRUE,
            legend.position = c("right","bottom"),
            inner.margins = 0.09)
tmap_save(map, "figures/map_sites.pdf", width = 6, height = 8)
# Size: 6*8
# Predicted burnt area fraction-- reconstructions ----

## fxTWAPLS

taxa_start<- which(colnames(pollen_bin) == "reference")+1
taxa_end<- length(colnames(pollen_bin))
core<- pollen_bin[,taxa_start:taxa_end]

fossil_charcoal <- fxTWAPLS::WAPLS.predict.w(fit_charcoal, core) #WAPLS
fossil_t_charcoal <- fxTWAPLS::TWAPLS.predict.w(fit_t_charcoal, core) #TWAPLS
fossil_f_charcoal <- fxTWAPLS::WAPLS.predict.w(fit_f_charcoal, core) #fxWAPLS
fossil_tf_charcoal <- fxTWAPLS::TWAPLS.predict.w(fit_tf_charcoal, core) #fxTWAPLS


data_info<- pollen_bin[,1:8]
data_info<-st_as_sf(data_info, coords = c("longitude","latitude"),crs=4326)

pre_core<- cbind(data_info,
                 (as.data.frame(BCTransformInverse(fossil_charcoal$fit[,4], lambda=0.25))), 
                 (as.data.frame(BCTransformInverse(fossil_t_charcoal$fit[,5], lambda=0.25))),
                 (as.data.frame(BCTransformInverse(fossil_f_charcoal$fit[,2], lambda=0.25))),
                 (as.data.frame(BCTransformInverse(fossil_tf_charcoal$fit[,4], lambda=0.25))))

start<- which(colnames(pre_core) == "reference")+1
end<- length(colnames(pre_core))-1
colnames(pre_core)[start:end]<- c("pre_wapls", "pre_twapls", "pre_fxwapls", "pre_fxtwapls")

pre_core<- pre_core %>% drop_na(c("pre_wapls", "pre_twapls", "pre_fxwapls", "pre_fxtwapls"))
# if value<-4 than na will appear. remove them

# write.csv(pre_core, "pre_core_rm.csv", row.names = T, fileEncoding = "UTF-8")
write.csv(pre_core, "pre_core_rm_BC.csv", row.names = T, fileEncoding = "UTF-8")

# Create a folder named rm_rare
dir.create("rm_rare")
# Composite plot: show the patterns in the Iberia region ####
Composite_plot(data=pre_core, method="pre_fxtwapls", hw=300, nreps=1000)
# SI: Compare 4 methods
Composite_plot(data=pre_core, method="pre_wapls", hw=300, nreps=1000)
Composite_plot(data=pre_core, method="pre_twapls", hw=300, nreps=1000)
Composite_plot(data=pre_core, method="pre_fxwapls", hw=300, nreps=1000)

curveout_wapls<-read.csv("rm_rare/pre_wapls_curveout.csv")
curveout_fxwapls<-read.csv("rm_rare/pre_fxwapls_curveout.csv")
curveout_twapls<-read.csv("rm_rare/pre_twapls_curveout.csv")
curveout_fxtwapls<-read.csv("rm_rare/pre_fxtwapls_curveout.csv")

# plot composite curves
curveout_plot(curveout_fxtwapls)
curveout_plot(curveout_wapls)
curveout_plot(curveout_twapls)
curveout_plot(curveout_fxwapls)
# *10^-3

# Palaeo spatial patterns of burnt area ####
target = pre_core %>%
  filter(age %in% c(100, 600, 3000, 7000, 9000, 11300, 12500, 13300, 14300))

colors <- brewer.pal(9, "Reds")
# setbreak <- c(0, 0.0001,0.0002,0.0003,0.0004,0.0005,0.0010,0.0015,0.0030,0.0060,0.0500)
setbreak <- c(0, 0.001,0.002,0.003,0.004,0.005,0.010,0.015,0.030,0.060,0.500)

time_map<- tm_shape(ip_shp_rgdal)+
  tm_polygons()+
  tm_shape(target) +
  tm_symbols(col = "pre_fxtwapls", 
             style = "cont",
             breaks= setbreak,   
             palette = "YlOrRd",
             border.col = "white", size=0.8, alpha = 0.9) + 
  tm_facets(by = "age", nrow = 3, free.coords = FALSE)+
  tm_layout(frame = FALSE, frame.lwd = NA, panel.label.bg.color = NA)
tmap_save(time_map, "figures/map_time.pdf", width = 8, height = 6)

# Compare micro and macro entities (SI) #### 
macro_all<- entitylist %>% dplyr::filter(rep == 1 & macro_micro=="macro")
macro_all_data<- charcoal_bin %>% dplyr::filter(entity_name %in% macro_all$charcoal_entity)
micro_all<- entitylist %>% dplyr::filter(rep == 1 & macro_micro=="micro")
micro_all_data<- charcoal_bin %>% dplyr::filter(entity_name %in% micro_all$charcoal_entity)
Composite_plot_char(data=micro_all_data, method="micro_all_data", hw=300, nreps=1000)
Composite_plot_char(data=macro_all_data, method="macro_all_data", hw=300, nreps=1000)

curveout_micro<-read.csv("rm_rare/micro_all_data_curveout.csv")
curveout_macro<-read.csv("rm_rare/macro_all_data_curveout.csv")
curveout_char_plot(curveout_macro)
curveout_char_plot(curveout_micro)

# Compare charcoal and reconstructed burnt area trends ####
# all charcoal(modern and fossil)--allc
# consider samples from shared age bins
allc<- entitylist %>% dplyr::filter(!is.na(pollen_entity))
rmallc<- allc %>% dplyr::filter(rep == 1 & macro_micro=="micro")
allc<- allc %>% dplyr::filter(charcoal_entity %notin% rmallc$charcoal_entity)

allc_data<- charcoal_bin %>% dplyr::filter(entity_name %in% allc$charcoal_entity)
pre_core<- pre_core %>% mutate(Entity.name=gsub("Abi 05/07", "Abi 05_07", Entity.name))
# so weird. abi 05/07 in allc cannot be recognize
allc<- allc %>% mutate(pollen_entity=gsub(allc$pollen_entity[1], "Abi 05_07", pollen_entity))
allc_recon<- pre_core %>% dplyr::filter(Entity.name %in% allc$pollen_entity)



giveid<- function(allc_data,allc_recon){
  namelist<- cbind(sort(unique(allc_data$site_name)),sort(unique(allc_recon$Site.Name)))
  colnames(namelist)<- c("char","ba")
  namelist<- as.data.frame(namelist) %>% mutate(siteid = row_number())
  allc_data<- allc_data %>% left_join(namelist,by = c("site_name"="char"))
  allc_recon<- allc_recon %>% left_join(namelist,by = c("Site.Name"="ba"))
  allc_merge<-  allc_data %>% inner_join(allc_recon, by=c("age", "siteid"))
  return(allc_merge)
}
allc_merge<- giveid(allc_data,allc_recon)
avg_allc_merge_recon <- allc_merge %>%filter(!is.na(pre_fxtwapls)) %>% group_by(age) %>%
  summarize(mean_recon = mean(pre_fxtwapls),
            median_recon= median(pre_fxtwapls))
avg_allc_merge_char <- allc_merge %>%filter(!is.na(max)) %>% group_by(age) %>%
  summarize(mean_charcoal = mean(max),
            median_charcoal= median(max))

allc_merge %>% filter(!is.na(pre_fxtwapls)) %>%filter(!is.na(max))%>%
  ggplot() +
  geom_smooth(method = "loess", span = 0.04, aes(age, median_recon*1000), data=avg_allc_merge_recon,col= "#E74C3C", se=F, size = 1, alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.04, aes(age, median_charcoal*10), data=avg_allc_merge_char,col= "#1E90FF", se=F, size = 1, alpha = 0.5) +
  xlab("Age (BP 1950)") + 
  scale_x_reverse(limits=c(15000,0), breaks = seq(15000,0, by = -2000))+
  theme_bw()+
  scale_y_continuous(#limits=c(0,1),
    # Features of the first axis
    name = "Max transformed charcoal*10",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~., name="Burnt area fraction*1000"))+
  theme(axis.title.y = element_text(color = "#1E90FF"),
        axis.title.y.right = element_text(color = "#E74C3C"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("figures/compare_results.pdf", width = 8, height = 6)

# CCA ####
# Use modern pollen records, environmental variables used in GLM and burnt area fraction from GFED4

modern_pollen<- p %>% dplyr::filter(median<=150)
## Transform to relative abundance
modern_pollen_ra<- modern_pollen[, which(colnames(modern_pollen) == "Abies"):length(modern_pollen)] 
modern_pollen_ra<- modern_pollen_ra/rowSums(modern_pollen_ra,na.rm = T)
modern_pollen_ra<- cbind(modern_pollen[, 1: which(colnames(modern_pollen) == "Abies")-1], modern_pollen_ra)

## Load GLM data
glmdata<- read.csv("~/Project/Palaeofire_v2/glmdataset.csv")
glmdata_mean<- glmdata %>% group_by(id) %>% summarise_each(mean)
## Load burnt area data
ba<- brick("~/Project/Palaeofire_v2/ip_0.25_burnedareafraction_annually_2001_2016.nc")
ba_mean <- calc(ba, fun = mean, na.rm = T)

mp_sf <- st_as_sf(modern_pollen_ra, coords = c("longitude", "latitude"),crs=4326)
mp_sf$cellno <- raster::extract(ba_mean, mp_sf,cellnumbers=T)[,"cells"]

ccadata<- glmdata_mean %>% filter(id %in% mp_sf$cellno)
mp_glm<- mp_sf %>% left_join(ccadata, by = c("cellno" = "id"))%>% 
  dplyr::select(-c(3:13))%>% 
  st_set_geometry(NULL)

mp_glm<- mp_glm %>% drop_na(Abies,burntareafraction,dtr,
                            drydays,wind,gpp,cropland,grazingland,urb_pop,nontreecover) 
mp_glm<- mp_glm[rowSums(mp_glm[, which(colnames(mp_glm) == "Abies"):which(colnames(mp_glm) == "Zygophyllaceae")])>0, ]

cca_pollen<- mp_glm[,which(colnames(mp_glm) == "Abies"):which(colnames(mp_glm) == "Zygophyllaceae")]
cca_glm<- mp_glm[,which(colnames(mp_glm) == "burntareafraction"):which(colnames(mp_glm) == "treecover")]
colnames(cca_glm)

## Apply transformation methods to reduce skewness
cca_glm[,13:14]<- cca_glm[,13:14]/100
cca_glm[,2:6]<- log(cca_glm[,2:6])
cca_glm[,2:6] <- cca_glm[,2:6] %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
cca_glm[,10:12] <- sqrt(cca_glm[,10:12])

## Start CCA
library(vegan)
mp.cca <- cca(cca_pollen ~ burntareafraction+dtr+drydays+wind+gpp+cropland+grazingland+urb_pop+nontreecover, data=cca_glm)
mp.cca2 <- cca(cca_pollen ~ burntareafraction, data=cca_glm)
mp.cca.clim <- cca(cca_pollen ~ burntareafraction+dtr+drydays+wind, data=cca_glm)
mp.cca.pop <- cca(cca_pollen ~ burntareafraction+cropland+grazingland+urb_pop, data=cca_glm)
mp.cca.veg <- cca(cca_pollen ~ burntareafraction+gpp+nontreecover, data=cca_glm)

summary(mp.cca)
summary(mp.cca2)


plot(mp.cca,type="n", xlim=c(-4, 4), ylim=c(-3, 2)) #
points(mp.cca, pch=16,scaling="species", dis="species", col="red", cex=0.8)
text(mp.cca, dis="bp",scaling="species", col="blue", cex=.8)
# text(mp.cca, dis="species",scaling="species", col="red", cex=0.8, pos=1)



summary(mp.cca.clim)
summary(mp.cca.pop)
summary(mp.cca.veg)

# Supplementary-- choose lambda####
pc_merge_modern$BC_0.1<- BCTransform(pc_merge_modern$weighteddata, lambda=0.1)
pc_merge_modern$BC_0.25<- BCTransform(pc_merge_modern$weighteddata, lambda=0.25)
pc_merge_modern$BC_0.33<- BCTransform(pc_merge_modern$weighteddata, lambda=0.33)
pc_merge_modern$BC_0.5<- BCTransform(pc_merge_modern$weighteddata, lambda=0.5)

lambda_test<- pc_merge_modern %>% dplyr::select(BC_0.1,BC_0.25,BC_0.33,BC_0.5)
lambda_test<- lambda_test %>% tidyr::gather("lambda","value")

library(hrbrthemes)
library(viridis)
lambda.labs <- c("λ = 0.1", "λ = 0.25", "λ = 0.33","λ = 0.5")
names(lambda.labs) <- c("BC_0.1", "BC_0.25","BC_0.33","BC_0.5")
lambda_test %>%
  ggplot( aes(x=value, color=lambda, fill=lambda)) +
  geom_histogram(alpha=0.6, binwidth = 0.1) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
    strip.background = element_rect(
      color="black", fill="white"
    )
  ) +
  xlab("Burnt area fraction (Box-Cox transformation)") +
  ylab("Count") +
  facet_wrap(~lambda, scales="free", labeller = labeller(lambda = lambda.labs))
ggsave("figures/compare_bc_lambda.png", width = 5, height = 5)

ggplot(pc_merge_modern, aes(x=weighteddata))+
  geom_histogram(binwidth=0.01, fill="grey", color="black", alpha=0.9)+
  xlab("Burnt area fraction")+
  ylab("Count")+
  theme_bw()
ggsave("figures/compare_raw.png", width = 3, height = 3)

# Get results using 0.1, 0.33, 0.5
w_charcoal<- pc_merge_modern$BC_0.5
start<- which(colnames(pc_merge_modern) == "Abies")
end<- which(colnames(pc_merge_modern) == "Ziziphus")
df_taxa<- pc_merge_modern[,(start):(end)] 
df_taxa<- df_taxa[, colSums(df_taxa != 0) > 0] ## Just to make sure we removed all the taxa that do not appear

library(fxTWAPLS)
fx_charcoal <- fxTWAPLS::fx(w_charcoal, bin = 0.002, show_plot = TRUE)
#fxtwapls
fit_tf_charcoal <- fxTWAPLS::TWAPLS.w(modern_taxa = df_taxa, modern_climate = w_charcoal, nPLS = 8, usefx = TRUE, fx = fx_charcoal)
cv_tf_charcoal<-fxTWAPLS::cv.w(df_taxa,w_charcoal,nPLS=8, fxTWAPLS::TWAPLS.w, fxTWAPLS::TWAPLS.predict.w,usefx=TRUE,fx=fx_charcoal)
rand_tf_charcoal<-fxTWAPLS::rand.t.test.w(cv_tf_charcoal,n.perm=999)
write.csv(rand_tf_charcoal, "process0311/methods_results_BC0.5_lambda.csv")
