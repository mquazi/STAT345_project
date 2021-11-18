
#######################################################################
########################### packages needed ###########################
#######################################################################
install.packages("stargazer")
install.packages("texreg")
install.packages("leaps")
install.packages("car")
install.packages("MASS")
install.packages("lmerTest")
install.packages("lmtest")
install.packages("lme4")
library(lme4)
library(lmtest)
library(leaps)
library(car)
library(leaps)
library(ggplot2)
library(stargazer)
library(texreg)
library(MASS)
library(car)
library(dae)
library(lmerTest)

#######################################################################
########################### input data ################################
#######################################################################
homedata<-read.csv(file.choose(),header=T) 

#################################################################################
########################### preliminary data analysis ###########################
#################################################################################
head(homedata)
summary(homedata)
table(homedata$Zip.Code)
table(homedata$Beds)
attach(homedata)
tapply(Price,Zip.Code,max)
which.max(Price)
homedata[93,]
summary(Price,Zip.Code)
table(homedata$Baths)
names(homedata)
boxplot(Price~Zip.Code,ylab="Price", xlab="Zip Code")
boxplot(Price~Region,ylab="Price", xlab="Region")

############################################################################
########################### categorial variables ###########################
############################################################################

x16<- Region=="NE"
x16<-as.numeric(Region=="NE")
x17<- Region=="NW"
x17<-as.numeric(Region=="NW")
x18<- Region=="SE"
x18<-as.numeric(Region=="SE")
table(Region)
table(x16)
table(x17)
table(x18)
############################################################################
########################### continuous variables ###########################
############################################################################
contdata<-cbind(Price,SqFt,Beds,Baths,Lat,Long)
head(contdata)
pairs(contdata)
cor(contdata)
#############################################################################
########################### checking selected correlations ##################
#############################################################################
plot(Price,x16)
plot(Price,x17)
plot(Price,x18)
cor(Lat,x16)
cor(Lat,x17)
cor(Lat,x18)
cor(Long,x16)
cor(Long,x17)
cor(Long,x18)
y<-Price
y
x19<-SqFt
x20<-Beds
x21<-Baths
x22<-Lat
x23<-Long


#######################################################################
########################### Model selection ###########################
#######################################################################
maindata<-cbind(y,x16,x17,x18,x19,x20,x21,x22,x23)
head(maindata)


### Backward elimination

upper<-formula(~x16+x17+x18+x19+x20+x21+x22+x23)
lower<-formula(~1)

mod1full<-lm(y~x16+x17+x18+x19+x20+x21+x22+x23)
mod1step<-step(mod1full,scope=list(lower=lower,upper=upper),direction="backward")
summary(mod1step)
backmodel<-lm(y~x19+x20)
vif(backmodel)

#### Forward
upper<-formula(~x16+x17+x18+x19+x20+x21+x22+x23)
lower<-formula(~1)
maindata<-as.data.frame(maindata)
mod1null<-lm(y~1,data=maindata)
mod1forward<-step(mod1null,scope=list(lower=lower,upper=upper),direction="forward")
summary(mod1forward)

### Stepwise
upper<-formula(~x16+x17+x18+x19+x20+x21+x22+x23)
lower<-formula(~1)

mod1<-lm(y~1,data=maindata)
mod1stepwise<-step(mod1full,scope=list(lower=lower,upper=upper),direction="both")
summary(mod1stepwise)


### BEST SUBSET

library(leaps)
### Cp
X<-data.frame(cbind(x16,x17,x18,x19,x20,x21,x22,x23))
cpsubset<-leaps(X,y,method="Cp",nbest=2)
cpsubset$which[order(cpsubset$Cp)[1:5],]
cpsubset$Cp[order(cpsubset$Cp)[1:5]]
cpsubset
cpsubset$Cp

### r2, adjr2, cp, bic
maindata<-data.frame(maindata)
leaps=regsubsets(y~x16+x17+x18+x19+x20+x21+x22+x23,data=maindata, nbest=2)
plot(leaps, scale="adjr2")
plot(leaps, scale="Cp")
plot(leaps, scale="bic")
plot(leaps, scale="r2")
summary(leaps)


### Final Model selected is
finalmodel<-lm(y~x19+x20)
summary(finalmodel)




#######################################################################
############################# Diagnostics #############################
#######################################################################
plot(finalmodel)
##major problem
resid<-residuals(finalmodel)
sum(resid)

# Non const var
bptest(y~x19+x20)
##p<0.05, p=3.048*10^-10
## we have non const error var

# Normality
par(mfrow=c(1,1))
hist(resid)
boxplot(resid,main="Boxplot residuals")
shapiro.test(resid)
##p<alpha, p=0.0004073

# Outliers
leverage<-hatvalues(finalmodel)
sum(leverage)
xoutliers<-which(leverage > 3*9/101)
plot(leverage)

rstud<-rstudent(finalmodel)
outlierTest(finalmodel)
rstuden<-abs(rstud)
youtliers<-which(rstuden >= qt(1-0.01/(2*101),91))
plot(cooks.distance(finalmodel),pch=23,main="Cook's distance")

# VIF
vif(finalmodel)

# transformation

library(MASS)
boxcox(y~x19+x20,lambda=seq(-2,2,length=10))
gmean<-exp(mean(log(y)))
sse<-c()
lambda<-c()
i<-1
for (lam in seq(-2,2,0.1)){
             if (lam !=0){
             k1<-(y^lam - 1)/(lam*gmean^(lam-1))
                    } else {
                             k1<-log(y)*gmean
                             }
             test<-anova(lm(k1~x19+x20))
             sse[i]<-test['Residuals','Sum Sq']
             lambda[i]<-lam
             i<-i+1
         }
cbind(lambda,sse)

abline(0,mean(sse),lty=2)
abline(mean(sse),0,lty=2)


ynew<-y^(0.5)
transformedmodel<-lm(ynew~x19+x20)
summary(transformedmodel)

# Diagnostics for transformed

plot(transformedmodel)
transresid<-residuals(transformedmodel)
plot(transresid,x19,main="Transformed Model Residuals vs SqFt") ##major problem

# Non const var
bptest(ynew~x19+x20)
##p<0.05, p=0.00735
##non const error var but much better than before

# Normality
par(mfrow=c(1,1))
hist(transresid)
boxplot(transresid,main="Boxplot Transformed model residuals")
shapiro.test(transresid)
##p<alpha, p=0.006 ##Not normal but major improvement

## Outliers
leverage<-hatvalues(transformedmodel)
sum(leverage)
plot(leverage)
xoutliers<-which(leverage > 3*9/101)
plot(leverage)

rstud<-rstudent(transformedmodel)
outlierTest(transformedmodel)
rstuden<-abs(rstud)
youtliers<-which(rstuden >= qt(1-0.01/(2*101),91))
plot(cooks.distance(transformedmodel),pch=23,main="Cook's distance")


# Dffits
dffits(transformedmodel)
maindata[which(dffits(transformedmodel)>0.5970223),]
maindata[which(dffits(finalmodel)>0.5970223),]
plot(dffits(transformedmodel),pch=2,ylab="DFFITS", xlab="Obs no.", main="DFFITS")
plot(dffits(finalmodel),pch=2,ylab="DFFITS", xlab="Obs no.", main="DFFITS")


# VIF
vif(transformedmodel)



###############################################################################
############################# Interval Estimation #############################
###############################################################################
newdata<-data.frame(x19=c(2554,1475),x20=c(3,3))
predict(finalmodel,newdata=newdata, interval="prediction",level=1-0.05/2)





