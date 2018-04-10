#### Script written by Gianalberto Losapio and Samuel Robinson

# Load everything ---------------------------------------------------------

library(bipartite)
library(igraph)
library(betalink)
library(vegan)
library(reshape2)
library(car)
library(MASS)

setwd("~/Projects/UBC/Alexandra Fiord (All)/Network analysis with GL/Files from GL")
load("networkData.Rdata") 

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
  mutate(Index=factor(Index,labels=c('Flower Diversity','Insect Diversity','Air Temperature (°C)','Network Complexity','Interaction Diversity'))) %>%
  ggplot(aes(date,Value))+geom_line()+facet_wrap(~ Index,strip.position='left',scales='free_y',ncol=1)+
  geom_point()+
  geom_text(data=data.frame(x=175,y=c(1.323166,1.887920,17.709014,3.248587,2.206528),
          Index=c('Flower Diversity','Insect Diversity','Air Temperature (°C)','Network Complexity','Interaction Diversity'),
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
setwd("~/Projects/UBC/Alexandra Fiord (All)/Network analysis with GL")
ggsave('network_metrics.eps',p,width=4,height=10)

detach("package:ggplot2", unload=TRUE)
detach("package:dplyr", unload=TRUE)
detach("package:tidyr", unload=TRUE)

###### try with species distribution models
library(lavaan)
library(semPlot)
library(nlsem)

mod.complx2

#Should the diversity~date relationships be hump-shaped? (Polynomial or spline fit? Not sure if this is possible using this type of SEM model)
M1 <- '
# measurement model, latent variable

###regressions
flowersdiv ~ date
bowlsdiv ~ date
#nicheoverpl ~ flowersdiv + bowlsdiv
shdiv ~ flowersdiv + bowlsdiv
complx ~ shdiv

# residual correlations
'

mo1 <- sem(M1, std.ov=T, missing="ml", estimator="ML", data=task)
summary(mo1, standardized=T, rsquare=T, fit.measures=T)

coef(mo1)
vcov(mo1)
parTable(mo1)
semPlotModel(mo1)

lab=c('Flower\nDiversity','Insect\nDiversity','Plant\nniche\noverlap','Insect\nniche\noverlap','Network\nDiversity','Network\nComplexity','Date')
semPaths(mo1,what='stand',whatLabels='stand',intercepts=F,layout='spring',shapeMan='rectangle',sizeMan=12,sizeMan2=8,
         nodeLabels=lab,
         label.cex=0.75,label.scale=F,posCol='darkblue',esize=10,
         edge.label.cex=1,edge.label.bg=T)

#Trying alternative SEMs ...

# #"Normal" quadratic calcuation
# task$date1=(task$date-min(task$date))+1
# task$date2=task$date1^2 #Quadratic form

#Orthogonal polynomial version - avoids correlation b/w Date and Date^2
task$date1=(task$date-min(task$date))+1
temp=poly(task$date1,2)
task$date1=temp[,1]
task$date2=temp[,2]
rm(temp)

#Potential model 1
M1.pot <- '
# measurement model, latent variable

###regressions
flowersdiv ~ date1 + date2
bowlsdiv ~ date1 + flowersdiv + temp
shdiv ~ flowersdiv + bowlsdiv + temp
complx ~ flowersdiv + bowlsdiv + temp

#orthogonal factors
temp~~0*date1
temp~~0*date2

# residual correlations
shdiv ~~ complx
'
#Potential model 2 - tried using clim =~ temp + sol as a latent variable, but this didn't improve the AIC, BIC, or other fit measures
#Trying model without temperature at all

M2.pot <- '
# measurement model, latent variable

###regressions
flowersdiv ~ date1 + date2
bowlsdiv ~ date1 + flowersdiv
shdiv ~ flowersdiv + bowlsdiv 
complx ~ flowersdiv + bowlsdiv 

# residual correlations
shdiv ~~ complx
'

M3.pot <-'
# measurement model, latent variable

###regressions
flowersdiv ~ date1 + date2
bowlsdiv ~ date1 + flowersdiv + clim
shdiv ~ flowersdiv + bowlsdiv + clim
complx ~ flowersdiv + bowlsdiv + clim

#latent variable
clim =~ temp + sol

# residual correlations
shdiv ~~ complx
'

#Fit models
mo1.pot <- sem(M1.pot, std.ov=T, missing="ml", estimator="ML", 
                data=task)
mo2.pot <- sem(M2.pot, std.ov=T, missing="ml", estimator="ML",
               data=task)
mo3.pot <- sem(M3.pot, std.ov=T, missing="ml", estimator="ML",
               data=task)

summary(mo1.pot, standardized=T, rsquare=T, fit.measures=T)
summary(mo2.pot, standardized=T, rsquare=T, fit.measures=T)
summary(mo3.pot, standardized=T, rsquare=T, fit.measures=T)


#Plot of mo1.pot
lab=c('Flower\nDiversity','Insect\nDiversity','Interaction\nDiversity','Network\nComplexity','Date','Date^2','Air\nTemp')

#define the layout
ly=matrix(c(-0.2,0.5, #Fldiv
            -0.2,-0.5, #bowldiv
            0.2,-0.5, #int div
            0.2,0.5, #net complx
            -0.5,0.2, #date
            -0.5,-0.2, #date2
            0.5,0), #temp
          ncol=2,byrow=TRUE)

semPaths(mo1.pot,what='std',whatLabels='std',intercepts=F,residuals=F,layout=ly,#ThreshAtSide=T,
         label.cex=0.75,label.scale=F,posCol='darkblue',esize=10,edge.label.cex=1,
         nodeLabels=lab,shapeMan='rectangle',sizeMan=12,sizeMan2=8,edge.label.position=0.4)

#Plot of mo3.pot

#define the layout
ly=matrix(c(0.5,0.5,# temp
            0.5,-0.5, #sol
            0.2,-0.5, #net div
            0.2,0.5, #net complx
            -0.2,0.5, #Fldiv
            -0.2,-0.5, #bowldiv
            -0.5,0.2, #date
            -0.5,-0.2, #date2
            0.5,0), #clim
          ncol=2,byrow=TRUE)

lab=c('Air\nTemp','Sun','Network\nDiversity','Network\nComplexity',
       'Flower\nDiversity','Insect\nDiversity',
       'Date','Date^2','Climate')

semPaths(mo3.pot,what='stand',whatLabels='stand',intercepts=F,residuals=F,nodeLabels=lab,
         layout=ly,
         label.cex=0.75,label.scale=F,posCol='darkblue',esize=10,edge.label.cex=1,
         shapeMan='rectangle',sizeMan=12,sizeMan2=8,edge.label.position=0.4)

# Longitudinal SEM --------------------------------------------------------

#Idea: measures of Flower diversity, Insect Diversity, Network Complexity, & Network Diversity
tempData <- subset(task,select=c(date,flowersdiv,bowlsdiv,shdiv,complx,temp))
tempData <- cbind(tempData,poly(tempData$date,2)) #Poly term for date
names(tempData)[7:8] <- c('date1','date2')

# #Should be installed from github:
# library(devtools)
# install_github('jslefche/piecewiseSEM@2.0',build_vignette=T)
library(piecewiseSEM)

#Model 1, temperature only affects network, no indirect path through insect/flower diversity
modlist1 <- psem(
    gls(bowlsdiv~flowersdiv,data=tempData,correlation=corCAR1(form=~date)),
    gls(shdiv~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date)),
    gls(complx~flowersdiv+bowlsdiv+temp,data=tempData,correlation=corCAR1(form=~date)),
    shdiv %~~% complx
)

summary(modlist1,conserve=T)
#Only tests separation of bowlsdiv and temp. Should add flDiv ~ temp test
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


# # Pseudo-r2 values seem unreliable. Deals with deviance, so meaning is different than regular r2 values
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

#Alternative way of calculating R2
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


# #Old method using lavaan 
# dateDist=as.matrix(dist(task$date)/2) #Distance matrix for dates, to show how far apart measurements are. (Done in 2-day increments)
# 
# #Reshape into wide format - looks like there's no easy way to do this in long form :(
# taskWide=task[,c('date','bowlsdiv','flowersdiv','complx','shdiv','temp')]
# taskWide=melt(taskWide,id='date')
# taskWide$meas=paste(taskWide$variable,taskWide$date,sep='.')
# taskWide=taskWide[,c('meas','value')]
# rownames(taskWide)=taskWide$meas #Rename
# taskWide$meas=NULL
# taskWide$value=as.numeric(taskWide$value)
# taskWide=data.frame(t(taskWide)) #Transpose
# names(taskWide)=as.character(taskWide[1,]) 
# 
# template='
# #Template
# i.175 ~~ 1*i.175 + i.c1*i.177 + i.c2*i.179 + i.c3*i.181 + i.c4*i.183 + 0*i.185 + 0*i.187 + 0*i.189 + 0*i.191 + 0*i.193 + 0*i.195 + 0*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.177 ~~ 1*i.177 + i.c1*i.179 + i.c2*i.181 + i.c3*i.183 + i.c4*i.185 + 0*i.187 + 0*i.189 + 0*i.191 + 0*i.193 + 0*i.195 + 0*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.179 ~~ 1*i.179 + i.c1*i.181 + i.c2*i.183 + i.c3*i.185 + i.c4*i.187 + 0*i.189 + 0*i.191 + 0*i.193 + 0*i.195 + 0*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.181 ~~ 1*i.181 + i.c1*i.183 + i.c2*i.185 + i.c3*i.187 + i.c4*i.189 + 0*i.191 + 0*i.193 + 0*i.195 + 0*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.183 ~~ 1*i.183 + i.c1*i.185 + i.c2*i.187 + i.c3*i.189 + i.c4*i.191 + 0*i.193 + 0*i.195 + 0*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.185 ~~ 1*i.185 + i.c1*i.187 + i.c2*i.189 + i.c3*i.191 + i.c4*i.193 + 0*i.195 + 0*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.187 ~~ 1*i.187 + i.c1*i.189 + i.c2*i.191 + i.c3*i.193 + i.c4*i.195 + 0*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.189 ~~ 1*i.189 + i.c1*i.191 + i.c2*i.193 + i.c3*i.195 + i.c4*i.197 + 0*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.191 ~~ 1*i.191 + i.c1*i.193 + i.c2*i.195 + i.c3*i.197 + i.c4*i.199 + 0*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.193 ~~ 1*i.193 + i.c1*i.195 + i.c2*i.197 + i.c3*i.199 + i.c4*i.201 + 0*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.195 ~~ 1*i.195 + i.c1*i.197 + i.c2*i.199 + i.c3*i.201 + i.c4*i.203 + 0*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.197 ~~ 1*i.197 + i.c1*i.199 + i.c2*i.201 + i.c3*i.203 + i.c6*i.209 + 0*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.199 ~~ 1*i.199 + i.c1*i.201 + i.c2*i.203 + i.c5*i.209 + i.c6*i.211 + 0*i.213 + 0*i.216 + 0*i.219
# i.201 ~~ 1*i.201 + i.c1*i.203 + i.c4*i.209 + i.c5*i.211 + i.c6*i.213 + 0*i.216 + 0*i.219
# i.203 ~~ 1*i.203 + i.c3*i.209 + i.c4*i.211 + i.c5*i.213 + i.c65*i.216 + 0*i.219
# i.209 ~~ 1*i.209 + i.c1*i.211 + i.c2*i.213 + i.c35*i.216 + i.c5*i.219
# i.211 ~~ 1*i.211 + i.c1*i.213 + i.c25*i.216 + i.c4*i.219
# i.213 ~~ 1*i.213 + i.c15*i.216 + i.c3*i.219
# i.216 ~~ 1*i.216 + i.c15*i.219
# i.219 ~~ 1*i.219
# 
# #Covariance terms using sigma(variance for measurement), rho(decay with distance: big rho = quick decay), and Distance in 2-day increments.
# i.c1==i.sigma*exp(-(i.rho^2)*(1^2))
# i.c15==i.sigma*exp(-(i.rho^2)*(1.5^2))
# i.c2==i.sigma*exp(-(i.rho^2)*(2^2))
# i.c25==i.sigma*exp(-(i.rho^2)*(2.5^2))
# i.c3==i.sigma*exp(-(i.rho^2)*(3^2))
# i.c35==i.sigma*exp(-(i.rho^2)*(3.5^2))
# i.c4==i.sigma*exp(-(i.rho^2)*(4^2))
# i.c5==i.sigma*exp(-(i.rho^2)*(5^2))
# i.c6==i.sigma*exp(-(i.rho^2)*(6^2))
# i.c65==i.sigma*exp(-(i.rho^2)*(6.5^2))
# '
# 
# Lmod1 <- '
# #Specified covariance matrix by hand
# bowlsdiv.175 ~~ 1*bowlsdiv.175 + bowlsdiv.c1*bowlsdiv.177 + bowlsdiv.c2*bowlsdiv.179 + bowlsdiv.c3*bowlsdiv.181 + bowlsdiv.c4*bowlsdiv.183 + 0*bowlsdiv.185 + 0*bowlsdiv.187 + 0*bowlsdiv.189 + 0*bowlsdiv.191 + 0*bowlsdiv.193 + 0*bowlsdiv.195 + 0*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.177 ~~ 1*bowlsdiv.177 + bowlsdiv.c1*bowlsdiv.179 + bowlsdiv.c2*bowlsdiv.181 + bowlsdiv.c3*bowlsdiv.183 + bowlsdiv.c4*bowlsdiv.185 + 0*bowlsdiv.187 + 0*bowlsdiv.189 + 0*bowlsdiv.191 + 0*bowlsdiv.193 + 0*bowlsdiv.195 + 0*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.179 ~~ 1*bowlsdiv.179 + bowlsdiv.c1*bowlsdiv.181 + bowlsdiv.c2*bowlsdiv.183 + bowlsdiv.c3*bowlsdiv.185 + bowlsdiv.c4*bowlsdiv.187 + 0*bowlsdiv.189 + 0*bowlsdiv.191 + 0*bowlsdiv.193 + 0*bowlsdiv.195 + 0*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.181 ~~ 1*bowlsdiv.181 + bowlsdiv.c1*bowlsdiv.183 + bowlsdiv.c2*bowlsdiv.185 + bowlsdiv.c3*bowlsdiv.187 + bowlsdiv.c4*bowlsdiv.189 + 0*bowlsdiv.191 + 0*bowlsdiv.193 + 0*bowlsdiv.195 + 0*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.183 ~~ 1*bowlsdiv.183 + bowlsdiv.c1*bowlsdiv.185 + bowlsdiv.c2*bowlsdiv.187 + bowlsdiv.c3*bowlsdiv.189 + bowlsdiv.c4*bowlsdiv.191 + 0*bowlsdiv.193 + 0*bowlsdiv.195 + 0*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.185 ~~ 1*bowlsdiv.185 + bowlsdiv.c1*bowlsdiv.187 + bowlsdiv.c2*bowlsdiv.189 + bowlsdiv.c3*bowlsdiv.191 + bowlsdiv.c4*bowlsdiv.193 + 0*bowlsdiv.195 + 0*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.187 ~~ 1*bowlsdiv.187 + bowlsdiv.c1*bowlsdiv.189 + bowlsdiv.c2*bowlsdiv.191 + bowlsdiv.c3*bowlsdiv.193 + bowlsdiv.c4*bowlsdiv.195 + 0*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.189 ~~ 1*bowlsdiv.189 + bowlsdiv.c1*bowlsdiv.191 + bowlsdiv.c2*bowlsdiv.193 + bowlsdiv.c3*bowlsdiv.195 + bowlsdiv.c4*bowlsdiv.197 + 0*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.191 ~~ 1*bowlsdiv.191 + bowlsdiv.c1*bowlsdiv.193 + bowlsdiv.c2*bowlsdiv.195 + bowlsdiv.c3*bowlsdiv.197 + bowlsdiv.c4*bowlsdiv.199 + 0*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.193 ~~ 1*bowlsdiv.193 + bowlsdiv.c1*bowlsdiv.195 + bowlsdiv.c2*bowlsdiv.197 + bowlsdiv.c3*bowlsdiv.199 + bowlsdiv.c4*bowlsdiv.201 + 0*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.195 ~~ 1*bowlsdiv.195 + bowlsdiv.c1*bowlsdiv.197 + bowlsdiv.c2*bowlsdiv.199 + bowlsdiv.c3*bowlsdiv.201 + bowlsdiv.c4*bowlsdiv.203 + 0*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.197 ~~ 1*bowlsdiv.197 + bowlsdiv.c1*bowlsdiv.199 + bowlsdiv.c2*bowlsdiv.201 + bowlsdiv.c3*bowlsdiv.203 + bowlsdiv.c6*bowlsdiv.209 + 0*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.199 ~~ 1*bowlsdiv.199 + bowlsdiv.c1*bowlsdiv.201 + bowlsdiv.c2*bowlsdiv.203 + bowlsdiv.c5*bowlsdiv.209 + bowlsdiv.c6*bowlsdiv.211 + 0*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.201 ~~ 1*bowlsdiv.201 + bowlsdiv.c1*bowlsdiv.203 + bowlsdiv.c4*bowlsdiv.209 + bowlsdiv.c5*bowlsdiv.211 + bowlsdiv.c6*bowlsdiv.213 + 0*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.203 ~~ 1*bowlsdiv.203 + bowlsdiv.c3*bowlsdiv.209 + bowlsdiv.c4*bowlsdiv.211 + bowlsdiv.c5*bowlsdiv.213 + bowlsdiv.c65*bowlsdiv.216 + 0*bowlsdiv.219
# bowlsdiv.209 ~~ 1*bowlsdiv.209 + bowlsdiv.c1*bowlsdiv.211 + bowlsdiv.c2*bowlsdiv.213 + bowlsdiv.c35*bowlsdiv.216 + bowlsdiv.c5*bowlsdiv.219
# bowlsdiv.211 ~~ 1*bowlsdiv.211 + bowlsdiv.c1*bowlsdiv.213 + bowlsdiv.c25*bowlsdiv.216 + bowlsdiv.c4*bowlsdiv.219
# bowlsdiv.213 ~~ 1*bowlsdiv.213 + bowlsdiv.c15*bowlsdiv.216 + bowlsdiv.c3*bowlsdiv.219
# bowlsdiv.216 ~~ 1*bowlsdiv.216 + bowlsdiv.c15*bowlsdiv.219
# bowlsdiv.219 ~~ 1*bowlsdiv.219
# 
# #Covariance terms using sigma(variance for measurement), rho(decay with distance: big rho = quick decay), and Distance in 2-day increments.
# bowlsdiv.c1==i.sigma*exp(-(bowlsdiv.rho^2)*(1^2))
# bowlsdiv.c15==i.sigma*exp(-(bowlsdiv.rho^2)*(1.5^2))
# bowlsdiv.c2==i.sigma*exp(-(bowlsdiv.rho^2)*(2^2))
# bowlsdiv.c25==i.sigma*exp(-(bowlsdiv.rho^2)*(2.5^2))
# bowlsdiv.c3==i.sigma*exp(-(bowlsdiv.rho^2)*(3^2))
# bowlsdiv.c35==i.sigma*exp(-(bowlsdiv.rho^2)*(3.5^2))
# bowlsdiv.c4==i.sigma*exp(-(bowlsdiv.rho^2)*(4^2))
# bowlsdiv.c5==i.sigma*exp(-(bowlsdiv.rho^2)*(5^2))
# bowlsdiv.c6==i.sigma*exp(-(bowlsdiv.rho^2)*(6^2))
# bowlsdiv.c65==i.sigma*exp(-(bowlsdiv.rho^2)*(6.5^2))
# '
# lsem1=lavaan(Lmod1,data=taskWide) 
# #PROBLEM: Error in lav_data_full(data = data, group = group, cluster = cluster,  : 
# #lavaan ERROR: some variables have only 1 observation or no finite variance



# Network turnover --------------------------------------------------------

## measure the difference between two adjacent networks (i.e. sampled in two consecutive days)

# S. dissimilarity in species composition
# OS dimmisimilarity in interactions between common species
# WN dissimilarity of all interactions
# ST dissimilarity due to species turnover

## all comparisons
allbeta=network_betadiversity(betanet, complete=FALSE)
# but we want only 2 consecutive days
betad = data.frame(di=names(betanet)[-length(betanet)], dj=names(betanet)[-1], S=NA, OS=NA, WN=NA, ST=NA)

riga = 1:(length(betanet)-1)
for(i in 1:(length(betanet)-1)){if(i!=1) riga[i]= (length(betanet))-(i-1)+riga[i-1]}

for(i in 1:nrow(allbeta)){if(i%in%riga) betad[which(i==riga),3:6] = allbeta[i,3:6]}

# contribution of species dissimilarity to network dissimilarity
betad$stwn = betad$ST/betad$WN
#Mean date
betad$date=rowMeans(cbind(as.numeric(as.character(betad$di)),as.numeric(as.character(betad$dj))))

## let's try to understand sth...

plot(betad$date,betad$S) # S. dissimilarity in species composition
plot(betad$date,betad$OS) # OS dimmisimilarity in interactions between common species
plot(betad$date,betad$WN) # WN -dissimilarity of all interactions (TOTAL DISSIMILARITY)
plot(betad$date,betad$ST) # ST dissimilarity due to species turnover (FROM MANUAL: "ST STRONGLY CONSTRAINED BY VALUES OF S, AND IS ONLY MEANINGFUL WHEN VALUES OF S ARE 'INTERMEDIATE'")
plot(betad$date,betad$stwn)

## definitively non linear trends....
with(betad,plot(date,S)) # S - dissimilarity in species composition
mod.s = lm(S~poly(date,2), data=betad)
Anova(mod.s)
summary(mod.s) # ok cool!!
## at the begginning the dissim in sp composition is high, then it decreased and by the end of the season it increased again
#SR: THIS MAKES SENSE TO ME. LARGE CHANGES ARE GOING ON DURING THE BEGINNING AND END AS INSECT AND PLANT SPECIES ARE COMING "ONLINE" OR "OFFLINE", SO THERE WOULD BE LARGE DAY-TO-DAY DIFFERENCES DURING THOSE PERIODS.

# OS - dissimilarity in interactions between common species
with(betad,plot(date,OS))
mod.os = lm(OS~poly(date,2), data=betad)
Anova(mod.os)
summary(mod.os) # this doesnt vary
#SR: QUESTION = DOES THIS MEAN THAT INTERACTIONS ARE DRIVEN BY UNCOMMON SPECIES AT BEGINNING AND END OF SEASON?

# WN - dissimilarity of all interactions
with(betad,plot(date,WN))
mod.wn = lm(WN~poly(date,2), data=betad)
Anova(mod.wn) #Date (p=0.002)
summary(mod.wn) # cool, similar pattern to S

# ST - dissimilarity due to species turnover
with(betad,plot(date,ST))
mod.st = lm(ST~poly(date,2), data=betad)
Anova(mod.st) #Date (p=0.05)
summary(mod.st) # cool, similar pattern to S
#SR: QUESTION = IT SEEMS LIKE THIS MIGHT BE MADE UNRELIABLE BY THE SPARSE NETWORK OF INSECTS. E.G. IF SPP1 WAS CAUGHT ON DAY 1 AND 3, BUT NOT 2, IT WOULD APPEAR THERE WAS MUCH MORE TURNOVER FROM DAY-TO-DAY. IS THIS SOMETHING I SHOULD WORRY ABOUT?

#stwn = ST/WN, i.e. proportion of total dissimilarity caused by spp turnover
#contribution of species dissimilarity to network dissimilarity
with(betad,plot(date,stwn))
mod.stwn = lm(asin(sqrt(stwn)) ~ poly(date,2), data=betad)
Anova(mod.stwn) #Date (p=0.11)
summary(mod.stwn) # nothing
#SR: QUESTION = THIS IMPLIES THAT SPP TURNOVER IS NOT VERY IMPORTANT SEASONALLY IN DAY-TO-DAY DIFFERENCES, CORRECT?

#Dissimilarity terms from betalink:
#S = spp.composition
#OS = interactions between common species,
#WN = all interactions
#ST = from species turnover
#From help files: "In the situations where S is either really high or really low, the values of ST are constrained and should no be given importance." #SR: QUESTION = I'M NOT SURE IF THIS IS IMPORTANT OR NOT
colMeans(betad[3:6],na.rm=T)

with(betad,plot(date,S,col='blue',pch=19,ylab='Dissimilarity',xlab='Day of season',ylim=c(0,1))) #Spp composition
with(betad,lines(date,S,col='blue'))
# with(betad,points(date,OS,col='red',pch=19)) #Common species
# with(betad,lines(date,OS,col='red'))
with(betad,points(date,WN,col='black',pch=19)) #All interactions
with(betad,lines(date,WN,col='black'))
with(betad,points(date,ST)) #Species turnover
with(betad,lines(date,ST,lty=2))
# legend('top',legend=c('S=Spp','OS=Com Spp','WN=All Int','ST=Spp Turn'),pch=c(19,19,19,1),col=c('blue','red','black','black'),title='Diss.Type')
legend('top',legend=c('S=Spp','WN=All Int','ST=Spp Turn'),pch=c(19,19,1),col=c('blue','black','black'),title='Diss.Type')

## thus, S WN and ST are quite interesting!
#SR: I THINK ALL OF THEM ARE FAIRLY INTERESTING! I'LL HAVE TO DO A BIT OF THINKING ABOUT WHY WE WOULD SEE SOME OF THESE PATTERNS RATHER THAN OTHERS, BUT THIS LOOKS LIKE AN EXCELLENT START!

#save.image("/g/pubblicazioni/samrobinson/networkData_analysis.Rdata")

#NEXT STEP: MATCH UP BOWL/FLOWER DIVERSITY, AND BOWL/FLOWER DISSIMILARITY, AND SHOW HOW/IF NETWORK DISSIMILARITY IS CONTROLLED BY COMMUNITY
#Reshape bowl data, remove days with no netting, get dissimilarity b/w days
bowlMat=reshape(bowlsDay,v.names='Count',idvar='Date',timevar='InsectSpecies',direction='wide')
names(bowlMat)=gsub('Count.','',names(bowlMat))
rownames(bowlMat)=bowlMat$Date
bowlMat=subset(bowlMat,(bowlMat$Date-1)>=min(as.numeric(as.character(betad$di)))&(bowlMat$Date-1)<=max(as.numeric(as.character(betad$dj)))) #Strips out days where insect netting didn't occur
bowlMat$Date=NULL
#Selects off-diagonal elements of distance matrix (pairwise distances b/w sites)
betad$betabowl=as.matrix(vegdist(bowlMat))[as.logical(rbind(0,cbind(diag(nrow(bowlMat)-1),0)))]

flwMat=reshape(flowersDay,v.names='Density',idvar='Date',timevar='PlantSpecies',direction='wide')
names(flwMat)=gsub('Density.','',names(flwMat))
rownames(flwMat)=flwMat$Date
flwMat=subset(flwMat,(flwMat$Date-1)>=min(as.numeric(as.character(betad$di)))&(flwMat$Date-1)<=max(as.numeric(as.character(betad$dj)))) #Strips out days where insect netting didn't occur
flwMat$Date=NULL
#Selects off-diagonal elements of distance matrix (pairwise distances b/w sites)
betad$betaflw=as.matrix(vegdist(flwMat))[as.logical(rbind(0,cbind(diag(nrow(flwMat)-1),0)))]
#Many days near end of season with 0 flowers, so NaN values are present in betaflw. Since we know there WERE spp around, replacing these NaN values with a small number seems prudent
betad$betaflw[is.nan(betad$betaflw)]=0.01

with(betad,plot(date,betaflw,col='blue',pch=19,ylab='Dissimilarity',xlab='Day of season',ylim=c(0,1))) #Flower betadiv
with(betad,lines(date,betaflw,col='blue'))
with(betad,points(date,betabowl,col='red',pch=19)) #Insect betadiv
with(betad,lines(date,betabowl,col='red'))
with(betad,points(date,WN,col='black',pch=19)) #Network betadiv
with(betad,lines(date,WN,col='black'))

mod.wn.betabowl=lm(WN~betabowl*betaflw,data=betad) #Nothing
mod.st.betabowl=lm(ST~betabowl*betaflw,data=betad) #Nothing
mod.s.betabowl=lm(S~betabowl*betaflw,data=betad) #Nothing
mod.os.betabowl=lm(OS~betabowl*betaflw,data=betad) #Not really. I think the sampling is too sparse to make any concrete statements. however...

#TRYING OUT SEM TO MODEL NETWORK TURNOVER

#Orthogonal polynomial version - avoids correlation b/w Date and Date^2
betad$date1=(betad$date-min(betad$date))+1
temp=poly(betad$date1,2)
betad$date1=temp[,1]
betad$date2=temp[,2]
rm(temp)

#Av Temperature b/w days
betad$temp=rowMeans(cbind(task$temp[1:(nrow(task)-1)],task$temp[2:(nrow(task))]))

#If we hypothesize:

M1.betadiv <- '
# measurement model, latent variable

bowls =~ betabowl
flws =~ betaflw

#regressions
WN ~ temp + bowls + flws #Network betadiv

#correlations

bowls~~flws

#orthogonal factors
'
mo1.betadiv <- sem(M1.betadiv, std.ov=T, missing="ml", estimator="ML", data=betad)

summary(mo1.betadiv, standardized=T, rsquare=T, fit.measures=T) #Chisq test is bad

semPaths(mo1.betadiv,what='std',whatLabels='std',intercepts=F,residuals=F,ThreshAtSide=T,#layout=ly,
         label.cex=0.75,label.scale=F,posCol='darkblue',esize=10,edge.label.cex=1,#nodeLabels=lab,
         shapeMan='rectangle',sizeMan=12,sizeMan2=8,edge.label.position=0.4)

#Summary: relationship b/w insect/flower turnover and network turnover isn't strongly controlled by anything except temperature. I don't think there's enough data to examine this in great detail.