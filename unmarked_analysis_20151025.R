#################################################
##### Set working directory, load libraries #####
##### and clean data                        #####
#################################################

setwd("~/Documents/*Science/*Research/Marianas/Point Counts/Data")

# Define functions

'%!in%' <- function(x,y)!('%in%'(x,y))

# load library
library(unmarked)
library(reshape)
library(plyr)
library(FD)
library(sciplot)
library(vegan)


# Read in data. Each row = one species at each point. 
birdcounts <- read.table("fulldata.txt", header=T) #load your data.
str(birdcounts)

# change format of date
birdcounts$date <- as.Date(birdcounts$date, "%m/%d/%y")

# make species levels lowercase
birdcounts$spp <- as.factor(as.character(tolower(birdcounts$spp)))

# make unique identifier for each count (point plus date) and each site date (site plus date)
birdcounts$pointdate <- as.factor(paste(birdcounts$point, birdcounts$date, sep=" "))
birdcounts$sitedate <- as.factor(paste(birdcounts$site, birdcounts$date, sep=" "))

# pull out point-level information, and put in alphabetical order by point 
# to be used later, after we use the formatDataDist function
point.data <- birdcounts[!duplicated(birdcounts[2:4]),]
point.data <- point.data[order(point.data$point),]
point.data <- point.data[,-c(5:9)]

# change from wide to long format
birdcounts.long <- untable(birdcounts[,c(1:6,8:11)], num=birdcounts[,7])

# remove data for certain species that we have to few estimates for (<10)
spp.included <- levels(birdcounts.long$spp)[which(levels(birdcounts.long$spp) %!in% c("coku","mela","gyal"))]
birdcounts.long <- birdcounts.long[which(birdcounts.long$spp %in% spp.included),]
birdcounts.long$spp <- as.factor(as.character(birdcounts.long$spp)) # makes r forget those other levels of this factor

# Remove the 6th occasion to achieve equal samples for each site (Tinian just has 5 occasions, rather than 6)
# Otherwise, the code interprets it as though Tinians just has super low densities
birdcounts.long <- subset(birdcounts.long, occasion <6)

birdcounts.long$occasion <- as.factor(birdcounts.long$occasion)

#################################################
##### Use unmarked package to get abundance #####
##### estimates for each point              #####
#################################################

# make empty output dataframes that can handle the data from predict()
lambda.out <- data.frame(spp=character(),Predicted=numeric(),SE=numeric(),lower=numeric(),upper=numeric(),point=character())
phi.out <- lambda.out #do the same for phi
det.out <- lambda.out #do the same for det

# Determine which spp were recorded on which site
present <- as.data.frame(ifelse(table(list(birdcounts.long$island,birdcounts.long$spp))>0,T,F))

# Make blank umdata to fill in
distance.cuttoffs <- c(0, 5, 10, 15, 20, 30, 40, 50, 100, 200)

blank.umdata <- formatDistData(birdcounts.long, 
                               distCol = "dist",
                               transectNameCol = "point",
                               dist.breaks = distance.cuttoffs,
                               occasionCol = "occasion")
blank.umdata <- ifelse(blank.umdata>-1,0,0)



# # Now for the models
# time1<-Sys.time()
# 
# for(i in 1:length(levels(birdcounts.long$spp))){
#   # create subset of data for one focal spp
#   sp.set <- subset(birdcounts.long, spp==levels(birdcounts.long$spp)[i])
#   
#   #### don't get estimates for islands where we didn't see the spp
#   sp.set$point <- as.factor(as.character(sp.set$point))
#   sp.set$island <- as.factor(as.character(sp.set$island))
#   sp.set$site <- as.factor(as.character(sp.set$site))
#   
#   #get into binned format, with each point seen as a replicate by the number of dates sampled (OccasioncCol="date")
#   sp.set.umdata <- formatDistData(sp.set, distCol="dist",
#                                   transectNameCol="point", 
#                                   dist.breaks=distance.cuttoffs,
#                                   occasionCol="occasion")
#   sp.set_rows <- row.names(sp.set.umdata)
#   
#   # Add in observed data into correct format with true zeros present
#   sp.blank.umdata <- blank.umdata
#   sp.blank.umdata[sp.set_rows,1:((length(distance.cuttoffs)-1)*5)] <- sp.set.umdata
#   sp.set.umdata <- sp.blank.umdata[,1:((length(distance.cuttoffs)-1)*5)]
#   # Now overwrite the row names
#   sp.set_rows <- row.names(sp.set.umdata)
#   
#   point.data.set <- point.data
#   point.data.set <- point.data.set[which(as.numeric(point.data.set$island) %in% which(present[,i]==T)),] # Will subset to islands where it's present later for the umdata file as well
#   point.data.set$point <- factor(point.data.set$point)
#   point.data.set$island <- factor(point.data.set$island)
#   point.data.set$site <- factor(point.data.set$site)
#   sp.set.umdata <- sp.set.umdata[which(as.character(point.data$island) %in% rownames(present)[which(present[,i]==T)]),] # This is where we subset to islands present
#   # Now overwrite the row names again!
#   sp.set_rows <- row.names(sp.set.umdata)
#   
#   sp.set.umfull <- cbind(sp.set.umdata,sp.set_rows, point.data.set)
#   
#   # Need to make it forget all the factor levels not represented in sp.set.umfull
#   sp.set.umfull$point <- factor(sp.set.umfull$point)
#   sp.set.umfull$island <- factor(sp.set.umfull$island)
#   sp.set.umfull$site <- factor(sp.set.umfull$site)
#   
#   #create dataframe in UMF format # island=sp.set.umfull$island, 
#   sp.set.umf <- unmarkedFrameGDS(y = sp.set.umdata,  #columns should be labeled dc1,dc2, dc3, dc4 and so on, for each distance class and values within are counts of that bird species at that point for that distance class. 
#                                  siteCovs = data.frame(site=sp.set.umfull$site, point=sp.set.umfull$point), #site covariates are probably island, maybe site, and maybe point (not sure how they deal with repeated measures)
#                                  dist.breaks = c(0, 5, 10, 15, 20, 30, 40, 50, 100, 200),  #these are the distances at which you enter a new distance class
#                                  numPrimary = 5,#numPrimary=ifelse(i==5,5,6), #number of surveys at each point
#                                  survey = "point", unitsIn = "m")  #survey style is point, not transect.  unitsIn are meters, I suspect
#   mod <- gdistsamp(~point, ~1, ~1,data=sp.set.umf, keyfun="halfnorm", output="density", unitsOut="ha",
#                    starts=c(1,rep(0,dim(sp.set.umfull)[1]-1),3,3),
#                    method = ifelse(i==7, "SANN", "BFGS")
#                    )
#   
#   # Make dataframe for predictions
#   newdf<-data.frame(point=factor(point.data.set$point))
#   
#   # Make the predictions
#   lambda <- predict(mod, newdata=newdf, type="lambda", appendData=TRUE)
#   phi <- predict(mod, newdata=newdf, type="phi", appendData=TRUE)
#   det <- predict(mod, newdata=newdf, type="det", appendData=TRUE)
#   
#   # Get ready to paste these on to the end of the output dataframes
#   spp <- rep(levels(birdcounts.long$spp)[i],dim(lambda)[1])
#   lambda <- cbind(spp,lambda)
#   phi <- cbind(spp,phi)
#   det <- cbind(spp,det)
#   
#   # Attach to ouput
#   lambda.out <- rbind(lambda.out,lambda)
#   phi.out <- rbind(phi.out,phi)
#   det.out <- rbind(det.out,det)
#   
#   print(warnings())
#   print(mod)
#   print(paste(i,"out of",length(levels(birdcounts.long$spp)),"done at", Sys.time(),sep=" "))
# }
# time2<-Sys.time()
# 
# 
# 
# write.csv(lambda.out,"lambda.out.csv")
# write.csv(phi.out,"phi.out.csv")
# write.csv(det.out,"det.out.csv")


lambda.out <- read.csv("lambda.out.csv",row.names=1)
phi.out <- read.csv("phi.out.csv",row.names=1)
det.out <- read.csv("det.out.csv",row.names=1)



# Convert to common name 4-letter codes
spp.converter<-c("apop"="mist",
                 "clma"="gowe",
                 "dima"="bldr",
                 "gaxa"="wtgd",
                 "gaga"="chic",
                 "mota"="timo",
                 "myru"="mihe",
                 "ptro"="mafd",
                 "rhru"="rufa",
                 "stbi"="iscd",
                 "toch"="cokf",
                 "zoco"="brwe")

lambda.out$spp<-revalue(lambda.out$spp,spp.converter)
phi.out$spp<-revalue(phi.out$spp,spp.converter)
det.out$spp<-revalue(det.out$spp,spp.converter)

#

fd.mat <- read.csv("FD_matrix_20141210_fewer_col.csv",row.names=1)
fd.mat <- read.csv("FD_matrix_20151026_fewer_col.csv",row.names=1)
fd.mat <- fd.mat[,1:10]# Remove the detailed foraging location data 
#fd.mat <- fd.mat[,-6] #No seed eaters
head(fd.mat)
#fd.mat <- read.csv("FD_matrix_20141210_many_col.csv",row.names=1)



abund <- lambda.out[,c(1,2,6)]
library(reshape)
abund.wide <- reshape(abund,direction="wide",v.names="Predicted",idvar="point",timevar="spp")
colnames(abund.wide) <- gsub("Predicted.","",colnames(abund.wide))

head(abund.wide)
abund.wide[is.na(abund.wide)] <- 0 # These are true zeros

# Make each point rowname and so the dataframe only has numeric abundance values
rownames(abund.wide) <- abund.wide[,1]
abund.wide <- abund.wide[,c(2:dim(abund.wide)[2])]

# Correct order
abund.wide <- abund.wide[,order(colnames(abund.wide))]
fd.mat <- fd.mat[order(rownames(fd.mat)),]
fd.mat <- fd.mat[,3:dim(fd.mat)[2]]# Get rid of "Scientific" and "English" names
dfvals <- dbFD(x=fd.mat,a=abund.wide)
dfvals$FDiv # I think these warnings are okay given the absence of some species

# Put together the dataframe for analysis

fd.df <- birdcounts[!duplicated(paste(birdcounts$island,birdcounts$point)),c("island","site","point")]
rownames(fd.df) <- fd.df$point

fd.df$total.abund <- rowSums(abund.wide)[rownames(fd.df)]
fd.df$shannon.div <- diversity(abund.wide)[rownames(fd.df)]
fd.df$FDiv <- dfvals$FDiv[rownames(fd.df)]
fd.df$FEve <- dfvals$FEve[rownames(fd.df)]
fd.df$FDis <- dfvals$FDis[rownames(fd.df)]
fd.df[,(dim(fd.df)[2]:(dim(fd.df)[2]+dim(dfvals$CWM)[2]))] <- dfvals$CWM[rownames(fd.df),] # This adds a bunch of columns for CWM values


# Relevel so Saipan is reference
fd.df$island <- relevel(fd.df$island, ref="Tinian")
fd.df$island <- relevel(fd.df$island, ref="Saipan")




head(fd.df)

library(lme4)
fd.mod0 <- lmer(FDiv~1+(1|site),dat=fd.df)
fd.mod <- lmer(FDiv~island+(1|site),dat=fd.df)
summary(fd.mod)
anova(fd.mod0,fd.mod) # Likelihood ratio test shows that island is a strong predictor of FD


abund.mod0 <- lmer(total.abund~1+(1|site),dat=fd.df)
abund.mod <- lmer(total.abund~island+(1|site),dat=fd.df)
summary(fd.mod)
anova(abund.mod0,abund.mod) # LRT shows that island is a strong predictor of total abundance

tapply(fd.df$FDiv,fd.df$island,mean)
tapply(fd.df$total.abund,fd.df$island,mean)
colnames(fd.df)
#tapply(fd.df$Fruit,fd.df$island,mean)
tapply(fd.df$Diet.Fruit,fd.df$island,mean)

colnames(fd.df)

#cols.to.test <- c(1:9,11)
cols.to.test <- c(4:13,15:16)
analysis.out <- as.data.frame(colnames(fd.df)[cols.to.test],strings.as.factors=F)
colnames(analysis.out)[1] <- "response"

analysis.out$saipan <- c(NA)
analysis.out$tinian <- c(NA)
analysis.out$rota <- c(NA)
analysis.out$saipan.se <- c(NA)
analysis.out$tinian.se <- c(NA)
analysis.out$rota.se <- c(NA)
analysis.out$lrt <- c(NA)

# 
for(i in cols.to.test){
  mod0 <- lmer(fd.df[,i]~1+(1|site),dat=fd.df)
  mod1 <- lmer(fd.df[,i]~island+(1|site),dat=fd.df)
  analysis.out[which(cols.to.test==i),2:4] <- fixef(mod1)
  analysis.out[which(cols.to.test==i),5:7] <- summary(mod1)$coefficients[,2]
  analysis.out[which(cols.to.test==i),8] <- anova(mod0,mod1)$"Pr(>Chisq)"[2]
}

# Species richness. This doesn't use the abundance data, so we're going back to the bidcounts data

#subset(birdcounts.long,birdcounts.long$pointdate=="Ltb 2011-07-19")
#head(birdcounts.long)
rich.dat <- birdcounts.long[!duplicated(birdcounts.long[,c(1:6)]),]
rich.dat$pointdate <- as.factor(as.character(rich.dat$pointdate))
richness.pointdate <- table(rich.dat$pointdate)
rich.no.duplicates <- rich.dat[!duplicated(rich.dat$pointdate),]
rich.no.duplicates <- rich.no.duplicates[order(rich.no.duplicates$pointdate),]
rich.no.duplicates <- cbind(rich.no.duplicates,richness.pointdate)
# Relevel so Saipan is reference
rich.no.duplicates$island <- relevel(rich.no.duplicates$island, ref="Tinian")
rich.no.duplicates$island <- relevel(rich.no.duplicates$island, ref="Saipan")

mod.rich0 <- lmer(Freq ~ 1 + (1|site),data=rich.no.duplicates)
mod.rich1 <- lmer(Freq ~ island + (1|site),data=rich.no.duplicates)
analysis.out[(dim(analysis.out)[1]+1),2:4] <- fixef(mod.rich1)
analysis.out[(dim(analysis.out)[1]),5:7] <- summary(mod.rich1)$coefficients[,2]
analysis.out[(dim(analysis.out)[1]),8] <- anova(mod.rich0,mod.rich1)$"Pr(>Chisq)"[2]
analysis.out$response <- as.character(analysis.out$response)
analysis.out[(dim(analysis.out)[1]),1] <- "richness"

# Make these absolute numbers rather than differences from the reference level
analysis.out$tinian <- analysis.out$saipan + analysis.out$tinian
analysis.out$rota <- analysis.out$saipan + analysis.out$rota
tail(analysis.out)
analysis.long <- reshape(analysis.out,varying=list(c(2:4),c(5:7)),idvar="response",direction="long")
colnames(analysis.long)[c(3:5)] <- c("island","estimate","se")
analysis.long$island <- c("Saipan","Tinian","Rota")[analysis.long$island]
analysis.long$island <- factor(analysis.long$island,levels=c("Saipan","Tinian","Rota"))


# Make a Community Weighted Mean plot

range1 <- function(x) range(x)[1]
range2 <- function(x) range(x)[2]

cwm <- c("Diet.Fruit","Diet.Inv","Diet.Nect","Diet.Vert","Diet.Seed")
islands <- c("Saipan","Tinian","Rota")

interval <- c(1,1,1,1,10)
x.vals <- c(1,1+cumsum(rep(interval,3)))
x.vals <- x.vals[1:(length(x.vals)-1)]
x.vals <- matrix(x.vals,nrow=5,byrow=F)

plot(x.vals,
     y=rep(-100,length(x.vals)),
     ylim=c(0,80),
     frame.plot=F,
     ylab="Community Weighted Mean",
     xlab="",
     las=1,
     xaxt="n",
     cex=2,
     pch=16)

for(i in 1:3) {
  fd.set <- subset(fd.df,island==islands[i])
  med <- apply(fd.set[,c("Diet.Fruit","Diet.Inv","Diet.Nect","Diet.Vert","Diet.Seed")],2,median)
  low <- apply(fd.set[,c("Diet.Fruit","Diet.Inv","Diet.Nect","Diet.Vert","Diet.Seed")],2,range1)
  high <- apply(fd.set[,c("Diet.Fruit","Diet.Inv","Diet.Nect","Diet.Vert","Diet.Seed")],2,range2)
  
  segments(x0=x.vals[,i],
           y0=low,
           y1=high)
  
  points(x.vals[,i],med,
         pch=c(21:25),
         bg="white")

}
axis(1, at = x.vals[3,], labels=islands)

legend("topleft",pch=c(21,25),
       legend=cwm)




###
### Five panel with Total Abundance, Shannon Diversity, Functional Diversity, Func Evenness
###
par(mfrow=c(3,2),pin=c(2,2.5))
abund.div.cols <- c("total.abund","richness","shannon.div","FEve","FDiv")
for(i in 1:length(abund.div.cols)){
  
  if(i==2){
    plot.new()
  }
  
  col.set <- subset(analysis.long,analysis.long$response==abund.div.cols[i])
  plot(estimate~c(1,2,3),data=col.set,
       xlim=c(0.75,3.25),
       ylim=c((min(col.set$estimate)-max(col.set$se)*1.5),(max(col.set$estimate)+max(col.set$se)*1.5)),
       frame.plot=F,
       ylab="",
       xlab="",
       las=1,
       xaxt="n",
       cex=2,
       pch=16)
  mtext(c("Total Abundance (birds/ha)","Richness","Shannon Diversity","Functional Evenness","Functional Diversity")[i],
        side=2,
        line=4,
        cex=0.75)
  segments(x0=c(1,2,3),y0=(col.set$estimate-col.set$se),y1=(col.set$estimate+col.set$se),
           lwd=10,
           col=rgb(0,0,0,alpha=0.4),
           lend="butt")
  if(i==length(abund.div.cols)){
    axis(1,at=c(1,2,3),labels=c("Saipan","Tinian","Rota"),cex.lab=1)
  }
  if(i==length(abund.div.cols)-1){
    axis(1,at=c(1,2,3),labels=c("Saipan","Tinian","Rota"),cex.lab=1)
  }
}






#################################################
##### Use unmarked package to get abundance #####
##### estimates for each island             #####
#################################################

lambda.island.out <- data.frame(spp=character(),Predicted=numeric(),SE=numeric(),lower=numeric(),upper=numeric(),point=character())
phi.island.out <- lambda.island.out #do the same for phi
det.island.out <- lambda.island.out #do the same for det

# # Now for the models
# time1<-Sys.time()
# 
# for(i in 1:length(levels(birdcounts.long$spp))){
#   # create subset of data for one focal spp
#   sp.set <- subset(birdcounts.long, spp==levels(birdcounts.long$spp)[i])
#   
#   #### don't get estimates for islands where we didn't see the spp
#   sp.set$point <- as.factor(as.character(sp.set$point))
#   sp.set$island <- as.factor(as.character(sp.set$island))
#   sp.set$site <- as.factor(as.character(sp.set$site))
#   
#   #get into binned format, with each point seen as a replicate by the number of dates sampled (OccasioncCol="date")
#   sp.set.umdata <- formatDistData(sp.set, distCol="dist",
#                                   transectNameCol="point", 
#                                   dist.breaks=distance.cuttoffs,
#                                   occasionCol="occasion")
#   sp.set_rows <- row.names(sp.set.umdata)
#   
#   # Add in observed data into correct format with true zeros present
#   sp.blank.umdata <- blank.umdata
#   sp.blank.umdata[sp.set_rows,1:((length(distance.cuttoffs)-1)*5)] <- sp.set.umdata
#   sp.set.umdata <- sp.blank.umdata[,1:((length(distance.cuttoffs)-1)*5)]
#   # Now overwrite the row names
#   sp.set_rows <- row.names(sp.set.umdata)
#   
#   point.data.set <- point.data
#   point.data.set <- point.data.set[which(as.numeric(point.data.set$island) %in% which(present[,i]==T)),] # Will subset to islands where it's present later for the umdata file as well
#   point.data.set$point <- factor(point.data.set$point)
#   point.data.set$island <- factor(point.data.set$island)
#   point.data.set$site <- factor(point.data.set$site)
#   sp.set.umdata <- sp.set.umdata[which(as.character(point.data$island) %in% rownames(present)[which(present[,i]==T)]),] # This is where we subset to islands present
#   # Now overwrite the row names again!
#   sp.set_rows <- row.names(sp.set.umdata)
#   
#   sp.set.umfull <- cbind(sp.set.umdata,sp.set_rows, point.data.set)
#   
#   # Need to make it forget all the factor levels not represented in sp.set.umfull
#   sp.set.umfull$point <- factor(sp.set.umfull$point)
#   sp.set.umfull$island <- factor(sp.set.umfull$island)
#   sp.set.umfull$site <- factor(sp.set.umfull$site)
#   
#   #create dataframe in UMF format # island=sp.set.umfull$island, 
#   sp.set.umf <- unmarkedFrameGDS(y = sp.set.umdata,  #columns should be labeled dc1,dc2, dc3, dc4 and so on, for each distance class and values within are counts of that bird species at that point for that distance class. 
#                                  siteCovs = data.frame(island=sp.set.umfull$island, point=sp.set.umfull$point), #site covariates are probably island, maybe site, and maybe point (not sure how they deal with repeated measures)
#                                  dist.breaks = c(0, 5, 10, 15, 20, 30, 40, 50, 100, 200),  #these are the distances at which you enter a new distance class
#                                  numPrimary = 5,#numPrimary=ifelse(i==5,5,6), #number of surveys at each point
#                                  survey = "point", unitsIn = "m")  #survey style is point, not transect.  unitsIn are meters, I suspect
#   if(length(levels(point.data.set$island))>1) {
#     mod <- gdistsamp(~island, ~1, ~1,data=sp.set.umf, keyfun="halfnorm", output="density", unitsOut="ha")#,
#     #starts=c(1,rep(0,dim(sp.set.umfull)[1]-1),3,3),
#     #method = ifelse(i==7, "SANN", "BFGS")
#     #)
#   } else {
#     mod <- gdistsamp(~1, ~1, ~1,data=sp.set.umf, keyfun="halfnorm", output="density", unitsOut="ha")
#   }
#   
#   
#   # Make dataframe for predictions
#   newdf<-data.frame(island=levels(point.data.set$island))
#   
#   # Make the predictions
#   lambda <- predict(mod, newdata=newdf, type="lambda", appendData=TRUE)
#   phi <- predict(mod, newdata=newdf, type="phi", appendData=TRUE)
#   det <- predict(mod, newdata=newdf, type="det", appendData=TRUE)
#   
#   # Get ready to paste these on to the end of the output dataframes
#   spp <- rep(levels(birdcounts.long$spp)[i],dim(lambda)[1])
#   lambda <- cbind(spp,lambda)
#   phi <- cbind(spp,phi)
#   det <- cbind(spp,det)
#   
#   # Attach to ouput
#   lambda.island.out <- rbind(lambda.island.out,lambda)
#   phi.island.out <- rbind(phi.island.out,phi)
#   det.island.out <- rbind(det.island.out,det)
#   
#   print(warnings())
#   print(mod)
#   print(paste(i,"out of",length(levels(birdcounts.long$spp)),"done at", Sys.time(),sep=" "))
# }
# time2<-Sys.time()
# 
# 
# 
# write.csv(lambda.island.out,"lambda.island.out.csv")
# write.csv(phi.island.out,"phi.island.out.csv")
# write.csv(det.island.out,"det.island.out.csv")


lambda.island.out <- read.csv("lambda.island.out.csv",row.names=1)
phi.island.out <- read.csv("phi.island.out.csv",row.names=1)
det.island.out <- read.csv("det.island.out.csv",row.names=1)

lambda.island.out$spp <- as.character(lambda.island.out$spp)
lambda.island.out$island <- as.character(lambda.island.out$island)

interval <- c(1,1,4)
x.vals <- c(1,1+cumsum(rep(interval,length(levels(birdcounts.long$spp)))))
x.vals <- x.vals[1:(length(x.vals)-1)]
x.vals <- matrix(x.vals,ncol=3,byrow=T)

par(mfrow=c(1,1),pin=c(10,6))
plot(-10,
     xlim=c(min(x.vals)-1,max(x.vals)+1),
     ylim=c(min(lambda.island.out$lower),max(lambda.island.out$upper)),
     frame.plot=F,
     ylab="",
     xlab="",
     las=1,
     xaxt="n",
     cex=2,
     pch=16)

levels(birdcounts.long$spp)[]
orderbird <- order(tapply(lambda.island.out$Predicted,lambda.island.out$spp,max),decreasing=T)
for(i in 1:length(levels(birdcounts.long$spp))) {
  sp.set <- subset(lambda.island.out,lambda.island.out$spp==levels(birdcounts.long$spp)[orderbird[i]])
  rownames(sp.set) <- sp.set$island
  if(all(is.na(sp.set["Rota",]))) sp.set["Rota",] <- NA
  if(all(is.na(sp.set["Tinian",]))) sp.set["Tinian",] <- NA
  if(all(is.na(sp.set["Saipan",]))) sp.set["Saipan",] <- NA
  sp.set <- sp.set[c("Saipan","Tinian","Rota"),]
  sp.set[is.na(sp.set)] <- 0 # These are true zeros
  
  segments(x0=x.vals[i,],y0=sp.set$lower,y1=sp.set$upper,
           lwd=10,
           lend="butt",
           col=c("grey30","grey60","grey90"))
  points(x.vals[i,],sp.set[,"Predicted"],
         pch=c(21,22,23),
         bg="white")
  
}
legend("topleft",pch=c(21,22,23),
       legend=c("Saipan","Tinian","Rota"),
       bty="n")
axis(1, at= x.vals[,2], labels = levels(birdcounts.long$spp)[orderbird])


#####
##### CCA
#####

abund.long <- lambda.out[,c(1,2,6)]

head(abund.long)

abund.wide <- reshape(abund.long, v.names="Predicted",idvar="spp",timevar="point",direction="wide")
rownames(abund.wide) <- abund.wide$spp
abund.wide <- abund.wide[,c(2:dim(abund.wide)[2])]
colnames(abund.wide) <- levels(abund.long$point)
abund.wide[is.na(abund.wide)] <- 0

abund.wide <- abund.wide[row.names(fd.mat),]

abund.wide

#data(fd.mat) # fd.mat


## Common but bad way: use all variables you happen to have in your
## environmental data matrix
fd.mat$num.foraging <- as.numeric(as.factor(fd.mat$ForagingLocation))
vare.cca <- cca(abund.wide,fd.mat[,c(1:6,8)])#[,c(1:5)])
vare.cca
par(pin=c(8,6))
plot(vare.cca)
plot(vare.cca, col=2,
     xlim=c(-1,1),
     ylim=c(-0.5,0.5))



