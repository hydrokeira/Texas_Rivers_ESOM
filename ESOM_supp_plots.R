setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/GGM_Paper2Materials/Codes and Data")

lulc<-read.csv("LULC_filled_interpolated_joinkey_simplified.csv")

param <- c("Bicarbonate","Calcium", "Chloride","Potassium","Magnesium", "Sodium","Sulfate", "Silica")
Rtypes <- c("Sdep", "Conglomerate", "Carbonates", "Water", "Evaporites",
            "Mudstone", "Igmet", "Sandstone")
LCtypes <- colnames(lulc)[4:10]
attributes <- c("date", "siteno","River.x.x","joinkey","lat.x","long.x",
                "area","TDS","Flow", "pH", "temp",
                param, Rtypes)
id.vars <- c("date", "siteno", "lat.x", "long.x", "River.x.x","joinkey", "area",
             "TDS", "Flow", "pH", "temp")

dat <- read.csv("CDAT_ESOM.csv")
dat$Year<-year(dat$date)
dat<-dat[,-c(35:41)]

lulc<-read.csv("LULC_filled_interpolated_joinkey_simplified.csv")
dat<-left_join(dat, lulc[2:10])

dat_unique<-dat[!duplicated(dat$joinkey),]

dat_melt<-melt(dat_unique[,c(4, 13, 28:34, 38:44)], id.vars = "River.x.x")

dat_melt_big<-subset(dat_melt, dat_melt$value > 0.05)

pal=c("Brazos"="#0071FE", "Red"="#C500FF", "Pecos"="#FF5500", "Colorado"="#55FF00")

pdf("Feature_Density_both_newLULC.pdf", width = 20, height = 7.5)

p1<-ggplot(dat_melt, aes(x=value, fill=River.x.x))+geom_density(alpha=0.5)+facet_wrap(~variable, scales = "free")+
  theme_classic()+theme(text = element_text(size = 12), legend.position = "null")+labs(x="Feature Value", y="Denisty", fill="", tag="a")+
  scale_fill_manual(values = pal)

p2<-ggplot(dat_melt_big, aes(x=value, fill=River.x.x))+geom_density(alpha=0.5)+facet_wrap(~variable, scales = "free")+
  theme_classic()+theme(text = element_text(size = 12))+labs(x="Feature Value", y="Denisty", fill="River Basin", tag="b")+
  scale_fill_manual(values = pal)

ggarrange(p1, p2, widths = c(0.45,0.55))

dev.off()

unique(dat$joinkey)

dat_count<-dat %>%
  group_by(River.x.x, year(date)) %>%
  count()

colnames(dat_count)<-c("River", "year", "num")

pdf("POR_Texas_ESOM.pdf", width = 10, height = 4)

ggplot(dat_count, aes(year, River, col=num))+geom_line(size=15)+theme_classic()+
  scale_color_gradientn(colors = c("lightblue", "dodgerblue", "midnightblue"))+
  labs(x="Year", y="River Basin", col="Number of Observations")+
  theme(text = element_text(size = 20), legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust=1, size = 12))

dev.off()

dat <- readRDS("allriv_compwts04212022.RDS")

count_df<-dat %>%
  dplyr::group_by(joinkey) %>%
  dplyr::count()

pdf("ObservationsPerSite.pdf", width = 6, height = 5)

ggplot(count_df, aes(n))+geom_histogram(fill="black")+theme_classic()+theme(text = element_text(size = 20))+
  labs(x="Number of Observations", y="Number of Sites")

dev.off()

count_df<-dat %>%
  group_by(joinkey) %>%
  count() %>%
  filter(n > 49)

colnames(count_df)[2]<-"tot"

count_year_df<-dat %>%
  filter(joinkey %in% count_df$joinkey) %>%
  group_by(joinkey, year(date)) %>%
  count() %>%
  filter(n > 5) %>%
  ungroup() %>%
  group_by(joinkey) %>%
  count() %>%
  filter(n > 9)

long_term_sites<-count_year_df$joinkey

count_all<-inner_join(count_df, count_year_df)

count_all$obs_year<-count_all$tot/count_all$n

cdat<-read.csv("CDAT_ESOM.csv")

cdat_shifting<-subset(cdat, cdat$joinkey %in% long_term_sites)

clust_sum<-cdat_shifting %>%
  group_by(joinkey, clust) %>%
  count() 

clust_sum2<-clust_sum %>%
  group_by(joinkey) %>%
  mutate(modal_clust=max(n), total_sum=sum(n), prop_mode=modal_clust/total_sum)

clust_unique<-clust_sum2[!duplicated(clust_sum2$joinkey),]

length(which(clust_unique$prop_mode > .5))

table(clust_sum$num_clust)

pdf("Cluster_Shifting.pdf", width = 5, height = 8)

ggplot(cdat_shifting, aes(date, joinkey, col=as.character(clust)))+geom_point()+theme_classic()+
  scale_color_manual(values=colors)+labs(x="Date", y="", col="Cluster")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 20))

dev.off()

res <- readRDS("ESOM_04242022.RDS")

#add unique ID for BMU
BestMatches <- as.data.frame(res$BestMatches)
#row of BM -1 times number of columns (82) plus column of BM (not sure why its like this..)
ids <- ((res$BestMatches[,2]-1)*length(unique(res$BestMatches[,3])))+res$BestMatches[,3]
#add unique ID to BMU
BestMatches<-cbind(BestMatches,ids)
names(BestMatches) <- c("obs", "row", "col", "id")

cdat_shifting_BM<-left_join(cdat_shifting, BestMatches)
cdat_shifting_BM$node<-paste0(cdat_shifting_BM$row, ",", cdat_shifting_BM$col)

length(unique(cdat_shifting_BM$node))/4100

length(unique(cdat_shifting_BM$joinkey))
