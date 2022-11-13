#### script to generate figures showing data availability across sites

library(ggplot2)
library(lubridate)
library(tidyr)
library(patchwork)

dat <- read.csv(paste0(intermediate_path,"model_data.csv"))
dat$date <- ymd(dat$date)

val_IDs <- unique(dat$SiteID[which(dat$val==T)])

dat <- dat[,c("date","SiteID","vmc_Surf","vmc_Deep","elev")]
dat <- dat[which(!is.na(dat$vmc_Surf)|!is.na(dat$vmc_Deep)),]

wide_dat_Surf <- pivot_wider(dat,
                        id_cols=date,
                        names_from=SiteID,
                        values_from=vmc_Surf)
wide_dat_Surf <- wide_dat_Surf[order(wide_dat_Surf$date),]
wide_dat_Surf$`All Sensors_0-5 cm` <- rowSums(!is.na(wide_dat_Surf[,names(wide_dat_Surf) != "date"]))
wide_dat_Surf$`Validation Sensors_0-5 cm` <- rowSums(!is.na(wide_dat_Surf[,names(wide_dat_Surf) %in% val_IDs]))


wide_dat_Deep <- pivot_wider(dat,
                             id_cols=date,
                             names_from=SiteID,
                             values_from=vmc_Deep) 
wide_dat_Deep <- wide_dat_Deep[order(wide_dat_Deep$date),]
wide_dat_Deep$`All Sensors_10-15 cm` <- rowSums(!is.na(wide_dat_Deep[,names(wide_dat_Deep) != "date"]))
wide_dat_Deep$`Validation Sensors_10-15 cm` <- rowSums(!is.na(wide_dat_Deep[,names(wide_dat_Deep) %in% val_IDs]))


plotdat <- merge(wide_dat_Surf[,c("date","All Sensors_0-5 cm","Validation Sensors_0-5 cm")],
                 wide_dat_Deep[,c("date","All Sensors_10-15 cm","Validation Sensors_10-15 cm")])

plotdat <- pivot_longer(plotdat,
                        -date,
                        names_to=c("summary","type"),
                        names_sep="_",
                        values_to="n")


# figure of number of functioning sensors over time
ggplot(plotdat[which(plotdat$summary=="All Sensors"),], aes(x=date,y=n,color=type)) + #linetype=summary,
  geom_line() +
  theme_classic() +
  scale_x_date(date_labels="%b - %Y") +
  scale_color_manual(values=c("10-15 cm"="black","0-5 cm"="grey60")) +
  #scale_linetype_manual(values=c("All Sensors"="solid","Validation Sensors"="dashed")) +
  scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,0.05))) +
  labs(y="Number of functioning sensors",
       x="Date",
       color="Sensor Depth") +
  theme(text=element_text(size=10),
        legend.title=element_blank(),
        legend.position=c(0.15,0.85),
        legend.margin=margin(0,1,0,1,"pt"),
        legend.spacing=unit(0,"pt"))
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

