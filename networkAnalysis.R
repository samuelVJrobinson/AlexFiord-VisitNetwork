#### Script written by Gianalberto Losapio and Samuel Robinson

# Load everything ---------------------------------------------------------

library(bipartite)
library(igraph)
library(betalink)
library(vegan)
library(reshape2)
library(car)
library(MASS)

load("./networkData.Rdata") 

netsDay$Date=as.numeric(format(as.Date(netsDay$Date),'%j')) #Changes date to Day of Year (DOY, 0-355)
bowlsDay$Date=as.numeric(format(as.Date(bowlsDay$Date),'%j')) 
flowersDay$Date=as.numeric(format(as.Date(flowersDay$Date),'%j')) 

#Strips out flowers and insects not involved in the network
visitors=as.character(unique(netsDay$InsectSpecies)) #Visiting insects
flowers=as.character(unique(netsDay$PlantSpecies)) #Visited flowers
flowers=flowers[-which(flowers=='Cerastium alpinum')] #Removes C. alpinum (visits, but no counts in control plots.)
bowlsDay=droplevels(subset(bowlsDay,bowlsDay$InsectSpecies%in%visitors))  #Drops from 748 to 418
flowersDay=droplevels(subset(flowersDay,flowersDay$PlantSpecies%in%flowers)) #276 to 184

### bowlsDay
### flowersDay
### netsDay
### weather

# global network
gnet = netsDay[which(netsDay$Count!=0),c(2,3,1)]
gnet$freq = netsDay[which(netsDay$Count!=0),4]/netsDay[which(netsDay$Count!=0),5]

# create a list with plant-poll nets for each sampling day
gnet2 = frame2webs(gnet, varnames = colnames(gnet), type.out = "list", emptylist = TRUE)

# example for visualizing (very basic)
# plotweb(gnet2[[3]])
# visweb(gnet2[[3]],labsize=0.5)

# create outcome dataframe
numdays=which(table(gnet$Date)>0,arr.ind=T) # 20 days over 23 with pollinators 
length(numdays)#SR: 19 ONCE REDUCED

#Uses DOY instead of just day of survey
DOY=as.numeric(names(which(table(gnet$Date)>0)))

task = data.frame(date=DOY,nicheoverpl=NA,nicheoverin= NA, cscorepl=NA, cscorein = NA)

#Generate networks for each day

for(i in 1:length(numdays)){

  task$nicheoverpl[i] = networklevel(gnet2[[numdays[i]]], index="niche overlap")[2] # among plants for pollinators, i.e. LL = lower level
  
  task$nicheoverin[i] = networklevel(gnet2[[numdays[i]]], index="niche overlap")[1] # among pollinators, i.e. HL = higher level
  #NOTE: I BELIEVE THESE WERE SWITCHED ORIGINALLY (I.E. HL<->LL). I CHANGED 1<->2, BUT G.L, CAN YOU VERIFY THIS?
  
  task$cscorepl[i] = networklevel(gnet2[[numdays[i]]], index="C score")[2] # among plants for pollinators, i.e. LL = lover level
  
  task$cscorein[i] = networklevel(gnet2[[numdays[i]]], index="C score")[1] # among pollinators, i.e. HL = higher level
  
  task$complx[i] = networklevel(gnet2[[numdays[i]]], index="linkage density")
  
  task$shdiv[i] = networklevel(gnet2[[numdays[i]]], index="Shannon diversity") 
  
  task$partdivpl[i] = networklevel(gnet2[[numdays[i]]], index= "partner diversity")[2]
  
  task$partdivin[i] = networklevel(gnet2[[numdays[i]]], index= "partner diversity")[1]
  
  # task$connect[i] = networklevel(gnet2[[numdays[i]]], index="connectance")
  # 
  # task$wconnect[i] = networklevel(gnet2[[numdays[i]]], index="weighted connectance")
}


bowls = acast(bowlsDay, Date ~ InsectSpecies)
bowlsdiv=vegan::diversity(bowls, "shannon") #Shannon diversity of bowl catches for each day

#Matches up bowl collection days with netting days (bowls were collected on the day after, so DOY is +1)
matchdays=which(as.numeric(names(bowlsdiv))-1<max(task$date+1)&as.numeric(names(bowlsdiv))-1>min(task$date-1))
task$bowlsdiv = bowlsdiv[matchdays]

flowers = acast(flowersDay, Date ~ PlantSpecies)
flowersdiv=vegan::diversity(flowers, "shannon")

matchdays=match(task$date,as.numeric(names(flowersdiv)))
task$flowersdiv = flowersdiv[matchdays]

sum.bowlsDay=with(bowlsDay,tapply(Count,Date,sum)) #Abundance from bowls for each day
sum.flowersDay=with(flowersDay,tapply(Density,Date,sum)) #Abundance of flowers for each day
sum.bowlsDay=data.frame(date=as.numeric(names(sum.bowlsDay)),count=unname(sum.bowlsDay)) #Convert to dataframe
sum.flowersDay=data.frame(date=as.numeric(names(sum.flowersDay)),count=unname(sum.flowersDay)) 

sum.flowersDay=sum.flowersDay[sum.flowersDay$date!=209,] #Strip out day 209 from flower data (no bowl data on that day)
abund=data.frame(date=sum.flowersDay$date,bowlcount=sum.bowlsDay$count,flwcount=sum.flowersDay$count) #Bowl and flower abundance data (all flowers, all insects)

matchdays=match(task$date,as.numeric(format(as.Date(weather$Date),'%j')))
task$temp = weather$AvTemp[matchdays]
task$sol = weather$AvSol[matchdays]
task$wind = weather$AvWind[matchdays]


# SEMs --------------------------------------------------------------------

#First, plots of model predictors:
library(ggplot2)
library(dplyr)
library(tidyr)

prestheme=theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15), 
                axis.title=element_text(size=20),
                title=element_text(size=20),
                panel.grid.major=element_line(size=0.5,colour='black',linetype='dotted'),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(size=1,colour='black'),
                strip.text=element_text(size=15))
theme_set(theme_bw()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

p=select(task,date,flowersdiv,bowlsdiv,complx,shdiv,temp) %>%
  gather('Index','Value',2:6) %>%
  mutate(Index=factor(Index,levels=c('flowersdiv','bowlsdiv','temp','complx','shdiv'))) %>%
  mutate(Index=factor(Index,labels=c('Flower Diversity','Insect Diversity','Air Temperature (?C)','Network Complexity','Interaction Diversity'))) %>%
  ggplot(aes(date,Value))+geom_line()+facet_wrap(~ Index,strip.position='left',scales='free_y',ncol=1)+
  geom_point()+
  geom_text(data=data.frame(x=175,y=c(1.323166,1.887920,17.709014,3.248587,2.206528),
          Index=c('Flower Diversity','Insect Diversity','Air Temperature (?C)','Network Complexity','Interaction Diversity'),
          Text=c('a','b','c','d','e')),
          aes(y=y,x=x,label=Text))+
  labs(x='Day of Year')+
  theme(axis.title.y = element_blank(),
        strip.background = element_rect(fill = 'transparent',colour='transparent'),
        strip.placement = 'outside',
        strip.text=element_text(size=14),
        axis.text=element_text(size=10), 
        axis.title=element_text(size=14),
        panel.grid.major=element_line(size=0.25,colour='white',linetype='solid'))
ggsave('network_metrics.eps',p,width=4,height=10)

detach("package:ggplot2", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:tidyr", unload=TRUE)

# Longitudinal SEM --------------------------------------------------------

#Idea: measures of Flower diversity, Insect Diversity, Network Complexity, & Network Diversity
tempData <- subset(task,select=c(date,flowersdiv,bowlsdiv,shdiv,complx,temp))
tempData <- cbind(tempData,poly(tempData$date,2)) #Poly term for date
names(tempData)[7:8] <- c('date1','date2')

# #Should be installed from github:
# library(devtools)
# install_github('jslefche/piecewiseSEM@2.0',build_vignette=T)
library(piecewiseSEM)
library(nlme)

#Model 1, temperature only affects network, no indirect path through insect/flower diversity
modlist1 <- psem(
    gls(bowlsdiv~flowersdiv,data=tempData,correlation=corCAR1(form=~date)),
    gls(shdiv~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date)),
    gls(complx~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date)),
    shdiv %~~% complx
)

summary(modlist1,conserve=T)

#Only tests separation of bowlsdiv and temp. Should add flDiv ~ temp test:
pvals <- dSep(modlist1,.progressBar=F,conditioning=T)$P.Value #Claims
pvals <- c(pvals,summary(gls(flowersdiv~temp,data=tempData,correlation=corCAR1(form=~date)))$tTable['temp','p-value'])
Cstat <- -2 * sum(log(pvals)) #4.69
DF <- 2 * length(pvals) #4 df
1 - pchisq(Cstat, DF) #p=0.319, model fits well
K <- do.call(sum, lapply(removeData(modlist1,formulas=1), function(i) attr(logLik(i), "df"))) #Log-lik params
Cstat+2*K #AIC = 4.69 + 2 * K = 30.69

#Calculate standardized coefficients
coefs <- stdCoefs(modlist1,tempData) #Extract coefficients
tempSD <- apply(tempData[,-1],2,sd) #Get SD for each column
coefs$Std.Estimate <- coefs$Estimate* #Standardized coefficients (coef * sd(x)/sd(y))
  (tempSD[match(coefs$Predictor,names(tempSD))]/tempSD[match(coefs$Response,names(tempSD))])


# # Pseudo-r2 values deal with deviance, so meaning is different than regular r2 values. Can be greater than 1 or less than 0.
# # Calculate pseudo-R2 for models (McFadden)
# pseudoR2 <- function(mod,adjust=F) {
#   nullLL <- logLik(update(mod,.~1))
#   if(adjust) K <- log(length(coef(mod))) else K <- 0 #Not sure this is the correct penalty term...
#   fullLL <- logLik(mod)
#   return(unlist(1-((fullLL-K)/nullLL)))
# }
# #fl div - LL is >1 for this,
# pseudoR2(gls(flowersdiv~date,data=tempData,correlation=corCAR1(form=~date),method='ML'))
# #Bowl div
# pseudoR2(gls(bowlsdiv~flowersdiv,data=tempData,correlation=corCAR1(form=~date),method='ML'))
# #shDiv
# pseudoR2(gls(shdiv~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date),method='ML'))
# #complx
# pseudoR2(gls(complx~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date),method='ML'))

#Alternative way of calculating R2 - how much variance explained by predicted value
# #Flower diversity, using a poly date term
# summary(lm(tempData$flowersdiv~predict(gls(flowersdiv~date1+date2,data=tempData,correlation=corCAR1(form=~date)))))

#Insect diversity
summary(lm(tempData$bowlsdiv~predict(gls(bowlsdiv~flowersdiv,data=tempData,correlation=corCAR1(form=~date)))))
#Interaction diversity
summary(lm(tempData$shdiv~predict(gls(shdiv~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date)))))
#Network complexity
summary(lm(tempData$complx~predict(gls(complx~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date)))))


#p-val for test bw beta coefficients for shdiv equation
1-pnorm((0.9218080-0.7388489)/sqrt((0.24673113^2)+(0.32225679^2)))

#p-val for test bw beta coefficients for complx equation
1-pnorm((0.9348146-0.8896561)/sqrt((0.23908793^2)+(0.30050616^2)))

#Model 2, temperature has indirect effect on network
modlist2 <- psem(
  gls(flowersdiv~temp,data=tempData,correlation=corCAR1(form=~date)),
  gls(bowlsdiv~flowersdiv+temp,data=tempData,correlation=corCAR1(form=~date)),
  gls(shdiv~flowersdiv+bowlsdiv,data=tempData,correlation=corCAR1(form=~date)),
  gls(complx~flowersdiv+bowlsdiv,data=tempData,correlation=corCAR1(form=~date)),
  shdiv %~~% complx
)
summary(modlist2,conserve=T)
#Contains all necessary claims, so no extra necessary.
#AIC = 58.34, C = 20.35, p < 0.001

#Since temperature is such an important force in controlling network structure and complexity, it seems like other models without it wouldn't really be that meaningful. Additionally, calculating AIC requires that the models be nested (use the same set of data), so dropping entire variables makes it hard to compare between them. See Section 2.5 of piecewiseSEM vignette.

#End
