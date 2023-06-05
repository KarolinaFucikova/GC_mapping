library("ggplot2")
library("phytools")

data.name = "data/GC_content.csv"
gc_data <- read.csv(file=data.name,header=TRUE)
## changing column names as needed
colnames(gc_data)
## if some columns start with numbers, they should be renamed like so:
##names(gc_data)[names(gc_data) == "X18S_GC"] <- "SSU_GC"

#excluding Trochiscia and M. homosphaera which don't have any environmental data
gc_clean <- gc_data[-c(13,30),]

# simple scatter plot of gc as a function of average temp in the warmest month
plot(gc_data$max_temp,gc_data$SSU_GC,na.rm=TRUE)
cor.test(gc_data$max_temp, gc_data$SSU_GC, method=c("pearson"))

# other exploratory plots are possible
#plot(gc_data$cp_GC_total,gc_data$cp_GC_coding,na.rm=TRUE)
#plot(gc_data$max_temp,gc_data$cp_GC_coding,na.rm=TRUE)
#plot(gc_data$max_temp,gc_data$cp_GC_total,na.rm=TRUE)
#plot(gc_data$cp_percent_coding,gc_data$cp_GC_coding,na.rm=TRUE)
#plot(gc_data$cp_size,gc_data$cp_GC_coding,na.rm=TRUE)
## not a lot of obvious correlation, except for total GC to coding GC

#chloroplast rrs vs. temperature
plot(gc_data$max_temp,gc_data$SSU_cp_GC)
cor.test(gc_data$max_temp, gc_data$SSU_cp_GC, method=c("pearson"))
# weak positive correlation, barely significant 0.04993, r = 0.2733

## 18S GC as function of temperature, divided by habitat type
fancyplot<-ggplot(gc_clean, aes(x=max_temp, y=SSU_GC, color=habitat, shape=habitat)) + 
  geom_point(size=5) + geom_smooth(method=lm, aes(fill=habitat))
fancyplot 
# the above adds a regression line separately for each habitat

# cleaner plotting
fancyplot1<-ggplot(gc_clean, aes(x=max_temp, y=SSU_GC)) + 
  geom_smooth(method=lm, se=FALSE, color="gray") +          
  geom_point(aes(color=habitat, shape=habitat), size=5) +
  geom_text(aes(label=ifelse(SSU_GC>54,as.character(species),''),hjust=1.2,vjust=1.5), size=5) +
  scale_color_manual(values=c('#56B4E9','#999999','#E69F00')) +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
  labs(x = "max temperature at site of origin (°C)",
       y = "% GC content in 18S gene")

fancyplot1

# grayscale plotting
fancyplot2<-ggplot(gc_clean, aes(x=max_temp, y=SSU_GC)) + 
  geom_smooth(method=lm, se=FALSE, color="gray") +          
  geom_point(aes(color=habitat, shape=habitat), size=5) +
  geom_text(aes(label=ifelse(SSU_GC>54,as.character(species),''),hjust=1.2,vjust=1.5), size=5) +
  scale_color_manual(values=c('#666666','#999999','#333333')) +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
  labs(x = "max temperature at site of origin (°C)",
       y = "% GC content in 18S gene")

fancyplot2

# correlation excluding the outlier strain (WJT) and also without taxa with no temp data
gc_noWJT <- gc_data[-c(1,13,30),]
# simple scatter plot of gc as a function of average temp in the warmest month
plot(gc_noWJT$max_temp,gc_noWJT$SSU_GC)
cor.test(gc_noWJT$max_temp, gc_noWJT$SSU_GC, method=c("pearson"))
fancyplot3<-ggplot(gc_noWJT, aes(x=max_temp, y=SSU_GC)) + 
  geom_smooth(method=lm, se=FALSE, color="gray") +          
  geom_point(aes(color=habitat, shape=habitat), size=5) +
  geom_text(aes(label=ifelse(SSU_GC>54,as.character(species),''),hjust=1.2,vjust=1.5), size=5) +
  scale_color_manual(values=c('#666666','#999999','#333333')) +
  theme(legend.key.size = unit(2, 'cm'), legend.title = element_text(size=18),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20)) + 
  labs(x = "max temperature at site of origin (°C)",
       y = "% GC content in 18S gene")

fancyplot3

## mapping GC content (or other traits) onto a tree
WJTtree <- "data/cp_nt_concatenated.con.tre"
WJT<-read.nexus(WJTtree)
#node number subtending OCC group happens to be 103 - this is not particularly easy to figure out
#I am rerooting the tree to have OCC as outgroup
# the current WJT can be viewed at anytime by simply running the line 
# plot(WJT)

WJT<-reroot(WJT,103)
plot(WJT)

# mapping traits and inferring ancestral states
gc_data <- read.csv(file=data.name,header=TRUE, row.names=1)
# pick the column to map and store in a new object
gc_ssu_tomap <- subset(gc_data,select=SSU_GC)
gc_cp_tomap <- subset(gc_data,select=cp_GC_total)
gc_cpcoding_tomap <- subset(gc_data,select=cp_GC_coding)

#turn the GC data into a vector
gc_ssu<-as.matrix(gc_ssu_tomap)[,1]
gc_cp<-as.matrix(gc_cp_tomap)[,1]
gc_cpcoding<-as.matrix(gc_cpcoding_tomap)[,1]

#estimate ancestral states, substitute objects from above
fit_ssu<-fastAnc(WJT,gc_ssu,vars=TRUE,CI=TRUE)
# anc.ML is the other method to infer ancestral states
# columns with missing data have to be handled separately
fit_cp<-fastAnc(WJT,gc_cp,vars=TRUE,CI=TRUE)
fit_cpcoding<-fastAnc(WJT,gc_cpcoding,vars=TRUE,CI=TRUE)
fit_cpcoding
# here we are mapping the GC content in chloroplast coding regions (fit_cpcoding) and overall chloroplast GC (fit_cp)

#pick variable to map, draw the tree
map <- contMap(WJT,gc_cpcoding,plot=FALSE)
plot(map,legend=0.5*max(nodeHeights(WJT)),fsize=c(0.7,0.7))

## the default rainbow palette has high as blue and low as red
## this function will swap the colors to be more intuitive
setMap<-function(x,...){
  if(hasArg(invert)) invert<-list(...)$invert
  else invert<-FALSE
  n<-length(x$cols)
  if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
  else x$cols[1:n]<-colorRampPalette(...)(n)
  x
}

plot(setMap(map,invert=TRUE))
## more making pretty tree - ladderize
map$tree<-ladderize.simmap(map$tree)

plot(setMap(map,invert=TRUE),fsize=c(0.7,0.7))

#picking ribosomal GC to map, draw the tree
map <- contMap(WJT,gc_ssu,plot=FALSE)
plot(map,legend=0.5*max(nodeHeights(WJT)),fsize=c(0.7,0.7))

# recolorize using the invert function above
plot(setMap(map,invert=TRUE))
## more making pretty tree - ladderize
map$tree<-ladderize.simmap(map$tree)

plot(setMap(map,invert=TRUE),fsize=c(0.7,0.7))
# in this case, WJT stands out as extremely high in ribosomal GC content

# plotting in black and white
bw.contMap<-setMap(map,c("white","black"))
plot(bw.contMap, lwd=2, fsize=c(0.6,0.6))
##########################

## a couple of ways to correct for the effect of phylogeny in the temp-GC regression
## PIC is one
## first drop taxa with missing temp data from the tree
gc_data <- read.csv(file="data/GC_content.csv",header=TRUE, row.names=1)
WJTtree <- "data/cp_nt_concatenated.con.tre"
WJT<-read.nexus(WJTtree)
WJT<-reroot(WJT,103)

miss<-which(is.na(gc_data[,2]))
WJT_clean <- drop.tip(WJT,miss)
plot(WJT_clean)
# and also drop them from the table
data_clean <- gc_data[-miss,]

pic.temp<-pic(data_clean$max_temp,WJT_clean)
pic.gc<-pic(data_clean$SSU_GC,WJT_clean)

fit.pic<-lm(pic.gc~pic.temp -1)
# the -1 1 specifies that the regression is through the origin 
# (the intercept is set to zero) as recommended by Garland et al., 1992.
# https://www.r-phylo.org/wiki/HowTo/Phylogenetic_Independent_Contrasts
fit.pic


summary(fit.pic)
plot(pic.temp, pic.gc,
     xlab="Temperature PIC",
     ylab="18S GC Content PIC",bg="grey",
     cex=1.8,pch=21, cex.lab=1.3, cex.axis=1.2)
abline(fit.pic,lwd=2,lty="dashed",col="grey")

### same thing for cp GC coding
gc_data <- read.csv(file="data/GC_content.csv",header=TRUE, row.names=1)
WJTtree <- "data/cp_nt_concatenated.con.tre"
WJT<-read.nexus(WJTtree)
WJT<-reroot(WJT,103)

miss<-which(is.na(gc_data[,2]))
WJT_clean <- drop.tip(WJT,miss)
plot(WJT_clean)
# and also drop them from the table
data_clean <- gc_data[-miss,]

pic.temp<-pic(data_clean$max_temp,WJT_clean)
pic.gc<-pic(data_clean$cp_GC_coding,WJT_clean)

fit.pic<-lm(pic.gc~pic.temp+0)
fit.pic

summary(fit.pic)
plot(pic.temp,pic.gc,xlab="PICs for average temperature in warmest month",
     ylab="PICs for chloroplast coding GC content",bg="grey",
     cex=1.4,pch=21)
abline(fit.pic,lwd=2,lty="dashed",col="red")

