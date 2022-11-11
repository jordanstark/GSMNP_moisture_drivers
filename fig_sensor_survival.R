#### script to generate figures showing data availability across sites

library(ggplot2)
library(lubridate)
library(tidyr)
library(patchwork)

dat <- read.csv(paste0(intermediate_path,"model_data.csv"))
dat$date <- ymd(dat$date)
dat <- dat[,c("date","SiteID","vmc_Surf","vmc_Deep","elev")]
dat <- dat[which(!is.na(dat$vmc_Surf)|!is.na(dat$vmc_Deep)),]

wide_dat_Surf <- pivot_wider(dat,
                        id_cols=date,
                        names_from=SiteID,
                        values_from=vmc_Surf)
wide_dat_Surf <- wide_dat_Surf[order(wide_dat_Surf$date),]
wide_dat_Surf$nsensors_Surf <- rowSums(!is.na(wide_dat_Surf[,names(wide_dat_Surf) != "date"]))

wide_dat_Deep <- pivot_wider(dat,
                             id_cols=date,
                             names_from=SiteID,
                             values_from=vmc_Deep) 
wide_dat_Deep <- wide_dat_Deep[order(wide_dat_Deep$date),]
wide_dat_Deep$nsensors_Deep <- rowSums(!is.na(wide_dat_Deep[,names(wide_dat_Deep) != "date"]))

plotdat <- merge(wide_dat_Surf[,c("date","nsensors_Surf")],
                 wide_dat_Deep[,c("date","nsensors_Deep")])

plotdat <- pivot_longer(plotdat,
                        starts_with("nsensors"),
                        names_to="type",
                        names_prefix="nsensors_",
                        values_to="n")

# figure of number of functioning sensors over time
ggplot(plotdat, aes(x=date,y=n,color=type)) +
  geom_line() +
  theme_bw() +
  scale_x_date(date_labels="%b - %Y") +
  labs(y="Number of functioning sensors",
       x="Date",
       color="Sensor Depth") +
  theme(text=element_text(size=10))
#ggsave(paste0(fig_path,"Number_of_sensors.png"),width=6,height=4)


# figure of survival of each sensor
plot2dat <- dat
plot2dat$`Data 0-5 cm` <- ifelse(is.na(plot2dat$vmc_Surf),"no data","data")
plot2dat$`Data 10-15 cm` <- ifelse(is.na(plot2dat$vmc_Deep),"no data","data")


plot2dat <- pivot_longer(plot2dat,
                         c(`Data 0-5 cm`,`Data 10-15 cm`),
                         names_to="depth",
                         values_to="data")
plot2dat <- plot2dat[which(plot2dat$data=="data"),]
plot2dat <- plot2dat[order(plot2dat$elev,decreasing=T),]
plot2dat$SiteID <- factor(plot2dat$SiteID,levels=unique(plot2dat$SiteID),
                          ordered=T)

ggplot(plot2dat, aes(x=date,y=SiteID)) +
  geom_point() +
  theme_bw() +
  scale_x_date(date_labels="%b - %Y") +
  facet_wrap(~depth)
ggsave(paste0(fig_path,"Sensor_survival.png"),width=5,height=9)

