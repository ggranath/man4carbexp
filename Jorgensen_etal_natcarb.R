##############################################################################
# R script for manuscript                                                    #
# "Forest management as a natural climate solution"                          #
# Jörgensen K, Granath G, Lindahl B, Strengbom J                             #
#                                                                            #
# Contact about R script: gustaf.granath@gmail.com                           #
##############################################################################

# This has been tested under:
# R version 3.5.3 (2019-03-11)
# Platform: x86_64-pc-linux-gnu (64-bit), Ubuntu 16.04.6 LTS

# required packages
# Statistical analyses
library(lme4) # ver 1.1-20
library(nlme) # ver 3.1-137
library(lmerTest) # ver 3.1-0
library(car) # ver 3.0-2
library(merTools) # ver 0.4.1

# Plotting and data handling
library(ggplot2) # ver 3.2
library(gridExtra) # ver 2.3
library(sjPlot) # ver 2.6.2
library(tidyr) # ver 0.8.2
library(dplyr) # ver 0.8.0.1

#  Perform imputations
library(Amelia) # # ver 1.7.5

# Load prepared data sets####
dat <- read.csv("jorgensen_etal_data.csv")
resp <- read.csv("jorgensen_etal_soil_resp.csv")

# Processing raw soil respiration data (not uploaded at the moment)
#source("Jorgensen_etal_prep_data.R")

# In these calcualtions we assume 50% carbon content of the biomass (hence, '*0.5'). 

# Figure 1 - C stock tree-soil,removed C, net ecosystem C seq ####

# Calculate mean N application and
# standardize so this is the N effect in models corrsponds
# to 88.1 g N per square meter
mean.N <- mean(dat[dat$N_g_m2>0,"N_g_m2"])
dat$N_mean_g_m2 <- dat$N_g_m2/mean.N
data.t <- dat
data.t$thin <- relevel(data.t$thin, ref = "thin")

#__C stock trees####
# C - stand tree carbon
# Divide by 1000 to get the units right
tree.stand <- lmerTest::lmer((tot.stand.end*0.5/1000) ~ 
                               N_mean_g_m2 * thin + P + (1|no_exp), 
                             data=data.t)

NewDat.treeStand <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                               P = c("no_P","no_P","no_P","no_P","Pxtra"))
NewDat.treeStand$thin <- relevel(NewDat.treeStand$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  NewDat.treeStand)
NewDat.treeStand$eff <- c(fixef(tree.stand)[-1] %*% t(X[,-c(1)]))
X[,1] <-c(1,0,0,0,0)
cis <- cbind(NewDat.treeStand$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.stand)  %*% x)*1.96)),
             NewDat.treeStand$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.stand)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
NewDat.treeStand <- cbind(NewDat.treeStand, cis)
NewDat.treeStand$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                           "Thinning + N + P")
NewDat.treeStand


#__C stock organic soil####
# *10/1000 to get the units right
soilC <- lmerTest::lmer((C_m2/1000)*10 ~  N_mean_g_m2 * thin + P + 
                          (1|no_exp), 
                        data=data.t)
summary(soilC)

newDat.soil <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                          P = c("no_P","no_P","no_P","no_P","Pxtra"))
newDat.soil$thin <- relevel(newDat.soil$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  newDat.soil)
newDat.soil$eff <- c(fixef(soilC)[-1] %*% t(X[,-c(1)]))
X[,1] <-c(1,0,0,0,0)
cis <- cbind(newDat.soil$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilC)  %*% x)*1.96)),
             newDat.soil$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilC)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
newDat.soil <- cbind(newDat.soil, cis)
newDat.soil$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                      "Thinning + N + P")
newDat.soil

#__Change in removed C####
# Changes in carbon removed by thinning. Testing thinnged vs thinned + N 
tree.removed <- lmerTest::lmer((pine_thin_stem + pine_thin_liv.branch + pine_thin_dead.branch + 
                                  spruce_thin_stem + spruce_thin_liv.branch + spruce_thin_dead.branch)*0.5/1000 ~
                                 N_mean_g_m2 + P + (1|no_exp), data=data.t[data.t$thin == "thin",])
summary(tree.removed)
tab_model(tree.removed)
newDat.remTree <- data.frame(N_mean_g_m2 = c(0,1,1), thin = c("thin", "thin", "thin"),
                          P = c("no_P","no_P", "Pxtra"))
newDat.remTree$thin <- relevel(newDat.remTree$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2  + P,
                  newDat.remTree)
newDat.remTree$eff <- c(fixef(tree.removed)[-1] %*% t(X[,-c(1)]))
X[,1] <-c(1,0,0)
cis <- cbind(newDat.remTree$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.removed)  %*% x)*1.96)),
             newDat.remTree$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.removed)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
newDat.remTree <- cbind(newDat.remTree, cis)
newDat.remTree$treatment <- c("Thinning", "Thinning + N", "Thinning + N + P")
newDat.remTree

#__Change in total net C uptake - excluding removed####
# Changes in total C stock. I.e. excluding C removed from the site. 
C.stand.stock <- lmerTest::lmer(((tot.stand.end*0.5/1000)+((C_m2/1000)*10)) ~  
                                  N_mean_g_m2 * thin + P + (1|no_exp), 
                                data=data.t)
summary(C.stand.stock)

NewDat.totstand <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                    P = c("no_P","no_P","no_P","no_P","Pxtra"))
NewDat.totstand$thin <- relevel(NewDat.totstand$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  NewDat.totstand)
NewDat.totstand$eff <- c(fixef(C.stand.stock)[-1] %*% t(X[,-c(1)]))
X[,1] <-c(1,0,0,0,0)
cis <- cbind(NewDat.totstand$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.stand)  %*% x)*1.96)),
             NewDat.totstand$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.stand)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
NewDat.totstand <- cbind(NewDat.totstand, cis)
NewDat.totstand$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                "Thinning + N + P")
NewDat.totstand

#__Change in total net C uptake-including removed####
c.tot <- lmerTest::lmer(((tot.prod.bio*0.5/1000)+((C_m2/1000)*10)) ~  
                          N_mean_g_m2 * thin + P + (1|no_exp), 
                        data=data.t)

newDat.ctot <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                          P = c("no_P","no_P","no_P","no_P","Pxtra"))
newDat.ctot$thin <- relevel(newDat.ctot$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  newDat.ctot)
newDat.ctot$eff <- c(fixef(c.tot)[-1] %*% t(X[,-c(1)]))
X[,1] <-c(1,0,0,0,0)
cis <- cbind(newDat.ctot$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(c.tot)  %*% x)*1.96)),
             newDat.ctot$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(c.tot)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
newDat.ctot <- cbind(newDat.ctot, cis)
newDat.ctot$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                           "Thinning + N + P")
newDat.ctot

#__Stat tables fig 1####
#____stats_models####
tab_model(list(tree.stand,soilC, tree.removed, C.stand.stock, c.tot), 
          pred.labels = c("Intercept (thinned)", "Nitrogen", "Unthinned", 
                          "P added (Unthinned,Nitrogen)", "Nitrogen x Unthinned"), 
          dv.labels = c("Standing tree C", "Soil C (t ha-1)", "Removed tree C", 
                        "Total net C uptake-excl removed",
                        "Total net C uptake-incl removed"), collapse.ci = TRUE,
          show.re.var = TRUE, show.icc = FALSE, show.r2 = FALSE,
          file="models_fig1.html")
# numbers
fig1.num <- rbind(NewDat.treeStand, newDat.soil, newDat.remTree, NewDat.totstand, newDat.ctot)
fig1.num$model <- c(rep(c("tree_stock","soil_stock","tree_removed",
                          "Tot C seq-excl removed", "Tot C seq-incl removed"), times=c(5,5,3,5,5)))
write.csv(fig1.num, "Figure1_table.csv")

#____Stats predictive values####
# Get treatment estimates
# stand tree carbon
mod.est.stand.tree <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                                 P = c("no_P","no_P","no_P","no_P","Pxtra"))
mod.est.stand.tree$thin <- relevel(mod.est.stand.tree$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  mod.est.stand.tree)
mod.est.stand.tree$estimate <- c(fixef(tree.stand) %*% t(X))
cis <- cbind(mod.est.stand.tree$estimate - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.stand)  %*% x)*1.96)),
             mod.est.stand.tree$estimate + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.stand)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
mod.est.stand.tree <- cbind(mod.est.stand.tree, cis)
mod.est.stand.tree$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                             "Thinning + N + P")
# soil C
mod.est.soil <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                           P = c("no_P","no_P","no_P","no_P","Pxtra"))
mod.est.soil$thin <- relevel(mod.est.soil$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  mod.est.soil)
mod.est.soil$estimate <- c(fixef(soilC) %*% t(X))
cis <- cbind(mod.est.soil$estimate - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilC)  %*% x)*1.96)),
             mod.est.soil$estimate + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilC)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
mod.est.soil <- cbind(mod.est.soil, cis)
mod.est.soil$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                       "Thinning + N + P")

# tree C removed
mod.est.rem.tree <- data.frame(N_mean_g_m2 = c(0,1,1), thin = c("thin", "thin", "thin"),
                           P = c("no_P","no_P","Pxtra"))
mod.est.rem.tree$thin <- relevel(mod.est.rem.tree$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 + P,
                  mod.est.rem.tree)
mod.est.rem.tree$estimate <- c(fixef(tree.removed) %*% t(X))
cis <- cbind(mod.est.rem.tree$estimate - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.removed)  %*% x)*1.96)),
             mod.est.rem.tree$estimate + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(tree.removed)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
mod.est.rem.tree <- cbind(mod.est.rem.tree, cis)
mod.est.rem.tree$treatment <- c("Thinning", "Thinning + N", "Thinning + N + P")

# Total net C uptake - excl removed
mod.est.net.stand.tot <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                               P = c("no_P","no_P","no_P","no_P","Pxtra"))
mod.est.net.stand.tot$thin <- relevel(mod.est.net.stand.tot$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  mod.est.net.stand.tot)
mod.est.net.stand.tot$estimate <- c(fixef(C.stand.stock) %*% t(X))
cis <- cbind(mod.est.net.stand.tot$estimate - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(C.stand.stock)  %*% x)*1.96)),
             mod.est.net.stand.tot$estimate + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(C.stand.stock)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
mod.est.net.stand.tot <- cbind(mod.est.net.stand.tot, cis)
mod.est.net.stand.tot$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                           "Thinning + N + P")

# Total net C uptake - incl removed
mod.est.tot.net <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                              P = c("no_P","no_P","no_P","no_P","Pxtra"))
mod.est.tot.net$thin <- relevel(mod.est.tot.net$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  mod.est.tot.net)
mod.est.tot.net$estimate <- c(fixef(c.tot) %*% t(X))
cis <- cbind(mod.est.tot.net$estimate - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(c.tot)  %*% x)*1.96)),
             mod.est.tot.net$estimate + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(c.tot)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
mod.est.tot.net <- cbind(mod.est.tot.net, cis)
mod.est.tot.net$treatment <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                          "Thinning + N + P")

# make csv of all estimates
fig1.est <- rbind(mod.est.stand.tree, mod.est.soil, mod.est.rem.tree, 
                  mod.est.net.stand.tot, mod.est.tot.net)
fig1.est$model <- c(rep(c("tree_stock","soil_stock","tree_removed",
                             "Tot C seq-excl removed", "Tot C seq-incl removed"), 
                           times=c(5,5,3,5,5)))
write.csv(fig1.est, "Figure1_estimate_table.csv")



# Figure 2 - soil C along latitudinal gradient ####
mean.N <- mean(dat[dat$N_g_m2>0,"N_g_m2"])
dat$N_mean_g_m2 <- dat$N_g_m2/mean.N
data.t <- dat
data.t$thin <- relevel(data.t$thin, ref = "thin")

# other models
#soilC <- lmerTest::lmer((tot.prod.bio*0.5/1000) ~ N_mean_g_m2 * thin + P + 
#                          N_mean_g_m2*lat_wgs84 + thin*lat_wgs84 +
#                          (1|no_exp), 
#                        data=data.t)
#summary(soilC)
# thin * latitude not important so removed for simplicity

#__a) - effect on soil C####
soilC <- lmerTest::lmer((C_m2/1000)*10 ~  N_mean_g_m2 * thin + P + 
                          N_mean_g_m2*lat_wgs84 +
                          + (1|no_exp), 
                        data=data.t)
summary(soilC)

# make plot
lat_grad <- range(data.t$lat_wgs84)
lat_grad_pred <- seq(lat_grad[1], lat_grad[2], 0.1)
newDat.soil <- data.frame(N_mean_g_m2 = rep(c(0,1,0,1,1),111), thin = rep(c("thin", "thin","no_thin","no_thin", "thin"),111),
                     P = rep(c("no_P","no_P","no_P","no_P","Pxtra"),111), lat_grad = lat_grad_pred)

newDat.soil$thin <- relevel(newDat.soil$thin, ref="thin")
X <- model.matrix( ~ N_mean_g_m2 * thin + P + N_mean_g_m2*lat_grad,
                   newDat.soil)
X[,5]=0

newDat.soil$eff <- c(fixef(soilC)[-1] %*% t(X[,-c(1)]))
X[,1] <-rep(c(1,0,0,0,0), 111)
cis <- cbind(newDat.soil$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilC)  %*% x)*1.96)),
             newDat.soil$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilC)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
newDat.soil <- cbind(newDat.soil, cis)
newDat.soil$vars <- rep(c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                     "Thinning + N + P"), 111)
# almost exact the same values for  No vthin+N and Thin+N+P. Hence we need to add a small 
# value so it shows on the plot
newDat.soil$eff <- with(newDat.soil, ifelse(vars == "No thinning + N", eff+0.075, eff))
newDat.soil$eff <- with(newDat.soil, ifelse(vars == "Thinning + N + P", eff-0.075, eff))

yl = expression(Effect ~on ~organic ~soil ~C ~sequestration ~ (t ~ha^{-1}))
cseq.soil.lat <- ggplot(newDat.soil, aes(x=lat_grad, y=eff, colour=vars)) +
  geom_line() +
  xlab("Latitude (degree North)") +
  ylab(yl) +
  labs(tag = "A") +
  theme_classic() +
  theme(axis.text = element_text(size=12, colour="black"),
    axis.title = element_text(size=14),
    #axis.title.x  = element_text(size=14),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.justification=c(0,1),
    legend.position=c(0,1),
    legend.box.background = element_rect(colour = "black")) +
  scale_y_continuous(limits = c(-2, 16), breaks=seq(-2,16, 2)) +
  scale_color_manual(name="Treatments",values=c("#E69F00", "#CC79A7", "black", "#009E73", "#0072B2"), 
                     breaks = c(c("No thinning + N", "Thinning + N + P","Thinning + N",
                                  "Thinning", "No thinning")), 
                     labels= c("No thinning + N", "Thinning + N + P","Thinning + N",
                               "Thinning (reference)", "No thinning"))
#cseq.soil.lat

#__b) effect on tree production####
totprod.c.lat <- lmerTest::lmer((tot.prod.bio*0.5/1000) ~ N_mean_g_m2 * thin + P + 
                                  N_mean_g_m2*lat_wgs84  + 
                                  (1|no_exp), 
                                data=data.t)
summary(totprod.c.lat)
#residuals splitted up into N added and not N added
ggplot(data.frame(lat =data.t$lat_wgs84, res= residuals(totprod.c.lat), Nt=data.t$N_mean_g_m2>0), aes(y=res, x=lat, color=Nt)) +
  geom_point()

lat_grad <- range(data.t$lat_wgs84)
lat_grad_pred <- seq(lat_grad[1], lat_grad[2], 0.1)
newDat <- data.frame(N_mean_g_m2 = rep(c(0,1,0,1,1),111), thin = rep(c("thin", "thin","no_thin","no_thin", "thin"),111),
                     P = rep(c("no_P","no_P","no_P","no_P","Pxtra"),111), lat_grad = lat_grad_pred)

newDat$thin <- relevel(newDat$thin, ref="thin")
X <- model.matrix( ~ N_mean_g_m2 * thin + P + N_mean_g_m2*lat_grad,
                   newDat)
X[,5]=0

newDat$eff <- c(fixef(totprod.c.lat)[-1] %*% t(X[,-c(1)]))
X[,1] <-rep(c(1,0,0,0,0), 111)
cis <- cbind(newDat$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(totprod.c.lat)  %*% x)*1.96)),
             newDat$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(totprod.c.lat)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
newDat <- cbind(newDat, cis)
newDat$vars <- rep(c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                     "Thinning + N + P"), 111)
# almost exact the same values for  No vthin+N and Thin+N+P. Hence we need to add a small 
# value so it shows on the plot
newDat$eff <- with(newDat, ifelse(vars == "No thinning + N", eff+0.075, eff))
newDat$eff <- with(newDat, ifelse(vars == "Thinning + N + P", eff-0.075, eff))

yl = expression(Effect ~on ~total ~tree ~C ~sequestration ~ (t ~ha^{-1}))
cseq.totprod.lat <- ggplot(newDat, aes(x=lat_grad, y=eff, colour=vars)) +
  geom_line() +
  xlab("Latitude (degree North)") +
  ylab(yl) +
  labs(tag = "B") +
  theme_classic() +
  theme(axis.text = element_text(size=12, colour="black"),
        axis.title = element_text(size=14),
        #axis.title.x  = element_text(size=14),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        legend.justification=c(0,1),
        legend.position=c(0,1),
        legend.box.background = element_rect(colour = "black")) +
  scale_y_continuous(limits = c(-5, 75), breaks=seq(-5,75, 10)) +
  scale_color_manual(name="Treatments",values=c("#E69F00", "#CC79A7", "black", "#009E73", "#0072B2"), 
                     breaks = c(c("No thinning + N", "Thinning + N + P","Thinning + N",
                                  "Thinning", "No thinning")), 
                     labels= c("No thinning + N", "Thinning + N + P","Thinning + N",
                               "Thinning (reference)", "No thinning"))
#cseq.totprod.lat

pdf("figure2_soil_tottree.pdf", width=6.5, height=10)
grid.arrange(cseq.soil.lat, cseq.totprod.lat, ncol=1, nrow =2)
dev.off()

#__Stat tables####
tab_model(list(soilC, totprod.c.lat), 
          pred.labels = c("Intercept (thinned)", "Nitrogen", "Unthinned", 
                          "P added (Unthinned,Nitrogen)", "Latitude",
                          "Nitrogen x Unthinned", "Nitrogen x Latitude"), 
          dv.labels = c("Soil C (t ha-1)", "Total tree net C"), collapse.ci = TRUE,
          show.re.var = TRUE, show.icc = FALSE, show.r2 = FALSE,
          file="models_fig2.html")

# Figure 3####
#__a) SOC versus productivity####
data.t <- dat
data.t$prod.C.ton.ha_yr <- data.t$tot.prod.bio*(0.5/10000)/data.t$age_last
data.t$C_m2.ton.ha <- (data.t$C_m2/1000)*10 

#data.t$real_treat <- relevel(data.t$real_treat, ref = "not_thinned+N")
SOCvsProd <- lmerTest::lmer(C_m2.ton.ha~ real_treat*prod.C.ton.ha_yr + 
                              (1|no_exp), 
                            data=data.t)
mod <- SOCvsProd
#data used in model
mod.data <- mod@frame

#if(is.null(vars.plot)) {vars.plot <- attr(terms(mod), "term.labels")}
vars <- attr(terms(mod), "term.labels")

# fix if there is an interaction and if so remove main term
#check.int <- unlist(strsplit(vars, "Species:"))
#if(any(duplicated(check.int)))  {#add <- check.int[anyDuplicated(check.int)]
vars = vars[-which(grepl(":", vars))] #}
p.var <- vars
#plots <- list()

#for (i in 3:length(vars)) {
#  if(!(vars[i] %in% vars.plot)) {next}
#  p.var = c("Species", vars[i])
xa <- p.var[2] # variable on x-axis
#if(grepl(":", xa)) {xa <- unlist(strsplit(xa,":"))[2]}
ya <- colnames(mod.data)[1] # response variable

p1 <- interactions::interact_plot(model=mod, pred = !! xa, modx = real_treat, plot.points = T,
                                  partial.residuals = T, centered="none",
                                  int.type = "confidence", interval = T)

p.dat <- ggplot_build(p1)$data[[3]]
l.dat <- ggplot_build(p1)$data[[1]]
l.dat$real_treat <- p1$data$real_treat

p.dat$real_treat <- data.t$real_treat
agg.dat <-  p.dat[, c("y", "x", "real_treat")]
#agg.dat$se <- NA
colnames(agg.dat)[c(1,2)] <- c(ya, xa)

yalab = expression(Soil~organic~carbon~(Mg))
xalab = expression(Net~tree~carbon~uptake ~ (Mg ~yr^{-1}))
fig3a <- ggplot(agg.dat, aes_string(x=xa, y=ya,color = "real_treat", shape = "real_treat")) +
  geom_point(size=2)+ #geom_linerange(show.legend = FALSE, alpha=0.2) +
  geom_line(data=l.dat, aes_string(y="y", x="x", color="real_treat"), 
            alpha =1, size = 1.1,inherit.aes = FALSE) +
  labs(y= yalab, 
       x= xalab) +
  labs(tag = "A") +
  ylim(c(10, 55)) +
  xlim(c(0, 0.35)) +
  scale_fill_manual(values=c("black", "red"), name="fill") +
  scale_color_manual(name = "Treatment", 
                     breaks = c("not_thinned", "not_thinned+N",
                                "thinning", "thinning+N","thinning+NP"), 
                     values = c("not_thinned" = "#E69F00", "not_thinned+N" = "#CC79A7",
                                "thinning" = "black", "thinning+N" = "#009E73","thinning+NP" = "#0072B2"),
                     labels = c("No thinning", "No thinning + N", 
                                "Thinning", "Thinning + N", "Thinning + N + P")) +
  scale_shape_manual(name = "Treatment", 
                     breaks = c("not_thinned", "not_thinned+N",
                                "thinning", "thinning+N","thinning+NP"), 
                     values = c(1,2,3,4,5), #c("not_thinned" = "#49b7fc", "not_thinned+N" = "#ff7b00",
                     #"thinning" = "#17d898", "thinning+N" = "#ff0083","thinning+NP" = "#0015ff"),
                     labels = c("No thinning", "No thinning + N", 
                                "Thinning", "Thinning + N", "Thinning + N + P")) +
    theme_classic() +
  theme(#legend.justification = c(0, 0), 
    legend.position = c(0.85, 0.81),
    #legend.key.size = unit(2, "line"),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12),
    #legend.position="none",
    legend.box.background = element_rect(colour = "black"),
    axis.text = element_text(size=12, color="black"),
    axis.title = element_text(size=14)) #+

fig3a

#____Stat tables####
# stats
tab_model(SOCvsProd, 
          pred.labels = c("Intercept (Unthinned)", "Nitrogen", "Thinned", 
                          "Nitrogen x Thinned","P (Thinned,Nitrogen)",
                          "Net tree C uptake",  
                          "Nitrogen x Net tree C uptake", "Thinned x Net tree C uptake", 
                          "Nitrogen x Thinned  x Net tree C uptake",
                          "P added (Thinned,Nitrogen) x Net tree C uptake"),
          dv.labels = c("Soil organic carbon (Mg)"), collapse.ci = TRUE,
          show.re.var = TRUE, show.icc = FALSE, show.r2 = FALSE,
          file="models_fig4.html")

#__b) soil respiration####
# save time series data at site 933
ts.dat <- resp[resp$no_exp == "933",]
ts.dat <- droplevels(ts.dat)

# Remove June and August flux measurements at site 933
resp.sub <- resp[!(resp$site_month == "933_june" | resp$site_month == "933_august"), ]

# centralise flux covariate data so they have mean zero
resp.sub$temp.ce <- scale(resp.sub$soil_temp_C, scale=FALSE)
resp.sub$wc.ce <- scale(resp.sub$wc_percent, scale=FALSE)
resp.sub <- droplevels(resp.sub) # remove if there are odd levels left

# soil water content is missing at two sites
# we use imputation to deal with this
# (1) we create many data sets with imputations
# (2) we run individual models
# (3) we average these models and account for extra variation due to imputation
mean.N <- mean(aggregate(N_g_m2 ~ site + plot, resp.sub[resp.sub$N_g_m2>0,], mean, na.rm=T)$N_g_m2)
resp.sub$N_mean_g_m2 <- resp.sub$N_g_m2/mean.N
resp.sub$thin <- relevel(resp.sub$thin, ref = "thin")

varIndex <- c("flux.quad", "wc.ce", "temp.ce", "site", "plot","N_mean_g_m2","thin")
IDS <- varIndex[-c(2,3)]
set.seed(1)
impute.out <- amelia(resp.sub[, varIndex], idvars = IDS, 
                     m = 1000)
modList.perarea <- lmerModList(log(flux.quad) ~ N_mean_g_m2*thin + wc.ce + temp.ce + (1|site/plot), 
                               data = impute.out$imputations)

# pool model output following the amelia package
# first coefs and their SEs
b.out<-NULL
se.out<-NULL
ran.out <- NULL
for(i in 1:length(modList.perarea)) {
  b.out <- rbind(b.out, fixef(modList.perarea[[i]]))
  se.out <- rbind(se.out, coef(summary(modList.perarea[[i]]))[,2])
  ran.out <- rbind(ran.out, c(unlist(VarCorr(modList.perarea[[i]])), attr(VarCorr(modList.perarea[[i]]), "sc")^2))
}
combined.results <- mi.meld(q = b.out, se = se.out)
area.resp.mod <- list(data.frame(coef=combined.results$q.mi[1,],
                            se=combined.results$se.mi[1,],
                            up=combined.results$q.mi[1,]+combined.results$se.mi[1,]*1.96,
                            lo=combined.results$q.mi[1,]-combined.results$se.mi[1,]*1.96,
                            tvalue=combined.results$q.mi[1,]/combined.results$se.mi[1,],
                            Pvalue=2*pt(-abs(combined.results$q.mi[1,]/combined.results$se.mi[1,]),df=Inf)),
                            data.frame(t(apply(ran.out,2,mean))),
                      data.frame(t(summary(modList.perarea[[1]])$ngrps), N=length(resid(modList.perarea[[1]]))))
colnames(area.resp.mod[[2]]) <- c("witin-site", "between-site","within-plot")
write.csv(data.frame(area.resp.mod[[1]],
                     area.resp.mod[[2]][,c(3,1,2)],area.resp.mod[[3]]), "models_fig3.csv")

# # quantify % decrease in unthinned stands
((exp(sum(area.resp.mod[[1]][,1][c(1,2,3,6)])) / 
    exp(area.resp.mod[[1]][,1][1]+ area.resp.mod[[1]][,1][3]))-1)*100

# Plot the mean effect and CI
# pool over all imputation models
# to get correct CIs
b.out<-NULL
se.out<-NULL
for (i in 1:length(modList.perarea)){
newDat <- data.frame(N_mean_g_m2 = c(0,1,0,1), thin = c("thin", "thin","no_thin","no_thin"))
newDat$thin <- relevel(newDat$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin ,
                  newDat)
newDat$eff <- c(area.resp.mod[[1]][-c(4,5),1] %*% t(X)) # remove temp and wc fixef(modList.perarea[[i]])
newDat$eff <- c(fixef(modList.perarea[[i]])[-c(4,5)] %*% t(X)) # remove temp and wc from varcov-matrix
X[,1] <-c(1,1,1,1)

vc <- vcov(modList.perarea[[i]])
ses <- apply(X, 1, function (x) as.vector(sqrt(x %*% vc[c(1,2,3,6),c(1,2,3,6)] %*% x)))
newDat <- cbind(newDat, ses)
newDat$vars <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N")
b.out <- rbind(b.out, newDat$eff)
se.out <- rbind(se.out, newDat$ses)
}
combined.results.plot <- mi.meld(q = b.out, se = se.out)
newDat$eff <- t(combined.results.plot$q.mi)
cis <- cbind(newDat$eff - newDat$ses*1.96,
              newDat$eff + newDat$ses*1.96)
colnames(cis) <- c("lo_95", "up_95")
newDat <- cbind(newDat, cis)
newDat[,c(3,6:7)] <- ((exp(newDat[,c(3,6:7)])/c(exp(newDat[1,3])))-1) *100

fig3b <- ggplot(newDat, aes(x=vars, y=eff)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, colour="black", width=0, position=position_nudge(x=0)) +
  geom_point(size=2, position=position_nudge(x=0)) +
  labs(tag = "B") +
  xlab("") +
  ylab("Effect on soil respiration (% change)") +
  theme_classic() +
  theme(axis.text  = element_text(size=12, colour="black"),
        axis.title.x  = element_text(size=14),
        legend.title=element_text(size=13)) +
  coord_flip() +
  scale_x_discrete(limits= c("Thinning", "Thinning + N", 
                             "No thinning", "No thinning + N"),
                   labels= c("Thinning", "Thinning + N", 
                             "No thinning", "No thinning + N")) +
  scale_y_continuous(limits = c(-35, 35), breaks=seq(-30, 35, 5))

pdf("figure3_prod_resp.pdf", width=6.5, height=10)
grid.arrange(fig3a, fig3b, ncol=1, nrow =2)
dev.off()

#____Stat tables####
write.csv(newDat, file="Figure3b_table.csv", row.names = FALSE)

# Figure 4 - NH4/NO3-N conc####
mean.N <- mean(dat[dat$N_g_m2>0,"N_g_m2"])
dat$N_mean_g_m2 <- dat$N_g_m2/mean.N
data.t <- dat
data.t$thin <- relevel(data.t$thin, ref = "thin")

# Test both N compounds

soilNO3 <- lmerTest::lmer(NO3.N_mg_g_soil ~  N_mean_g_m2 * thin + P +
                            (1|no_exp), 
                          data=data.t)
soilNH4 <- lmerTest::lmer(NH4.N_mg_g_soil ~  N_mean_g_m2 * thin + P +
                            (1|no_exp), 
                          data=data.t)

NewDat.no3 <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                         P = c("no_P","no_P","no_P","no_P","Pxtra"))
NewDat.no3$thin <- relevel(NewDat.no3$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  NewDat.no3)
NewDat.no3$eff <- c(fixef(soilNO3)[-1] %*% t(X[,-c(1)]))
X[,1] <-c(1,0,0,0,0)
cis <- cbind(NewDat.no3$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilNO3)  %*% x)*1.96)),
             NewDat.no3$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilNO3)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
NewDat.no3 <- cbind(NewDat.no3, cis)
NewDat.no3$vars <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                     "Thinning + N + P")

NewDat.nh4 <- data.frame(N_mean_g_m2 = c(0,1,0,1,1), thin = c("thin", "thin","no_thin","no_thin", "thin"),
                         P = c("no_P","no_P","no_P","no_P","Pxtra"))
NewDat.nh4$thin <- relevel(NewDat.nh4$thin, ref="thin")
X <- model.matrix(~N_mean_g_m2 * thin + P,
                  NewDat.nh4)
NewDat.nh4$eff <- c(fixef(soilNH4)[-1] %*% t(X[,-c(1)]))
X[,1] <-c(1,0,0,0,0)
cis <- cbind(NewDat.nh4$eff - apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilNH4)  %*% x)*1.96)),
             NewDat.nh4$eff + apply(X, 1, function (x) as.vector(sqrt(x %*% vcov(soilNH4)  %*% x)*1.96)))
colnames(cis) <- c("lo_95", "up_95")
NewDat.nh4 <- cbind(NewDat.nh4, cis)
NewDat.nh4$vars <- c("Thinning", "Thinning + N", "No thinning", "No thinning + N", 
                     "Thinning + N + P")

soilNeff <- ggplot(NewDat.no3, aes(x=vars, y=eff)) + 
  geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, colour="black", width=0, position=position_nudge(x=0.1)) +
  geom_point(size=2, pch=21, aes(fill="white"), position=position_nudge(x=0.1)) +
  geom_point(data=NewDat.nh4, aes(x=vars, y=eff, fill="black"), position=position_nudge(x=-0.1), 
             size=2, pch=21) +
  geom_errorbar(data=NewDat.nh4, aes(ymin=lo_95, ymax=up_95), 
                lwd=0.7, colour="black", width=0, position=position_nudge(x=-0.1)) +
  #  ylim(c(-10, 100)) +
  xlab("") +
   ylab(expression(Effect ~on ~soil ~NH[4]^{"+"} ~and ~NO[3]^{"-"} ~ (mg ~g^{-1}))) +
  theme_classic() +
  theme(axis.text  = element_text(size=12, colour="black"),
        axis.title.x  = element_text(size=14),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        #legend.key.size = 3,
        legend.justification=c(0.90,.3), 
        legend.position=c(0.9,.3)) +
  coord_flip() +
  scale_x_discrete(limits= c("Thinning", "Thinning + N", "Thinning + N + P",
                             "No thinning", "No thinning + N"),
                   labels= c("Thinning", "Thinning + N", "Thinning + N + P",
                             "No thinning", "No thinning + N")) +
  scale_y_continuous(limits = c(-0.02, 0.05), breaks=seq(-0.02, 0.05, 0.01)) +
  scale_fill_manual(name="Model",values=c("black", "white"), breaks = c("white", "black"), 
                    labels= c(expression(NO[3]^{"-"}), expression(NH[4]^{"+"})))

ggsave("Figure4_NHNO_N.png", soilNeff, height=10, width=18, units="cm") 
pdf("Figure4_NHNO_N.pdf", width=6.5, height=5)
grid.arrange(soilNeff, ncol=1, nrow =1)
dev.off()

# Supp mtrl####
#__Stat tables####
# stats
tab_model(list(soilNO3, soilNH4), 
          pred.labels = c("Intercept (thinned)", "Nitrogen", "Unthinned", 
                          "P added (Unthinned,Nitrogen)", "Nitrogen x Unthinned"), 
          dv.labels = c("Soil NH4+-N (mg g-1)", "NO3-N"), collapse.ci = TRUE,
          show.re.var = TRUE, show.icc = FALSE, show.r2 = FALSE,digits = 4,
          file="models_fig5.html")

# numbers in figure
# numbers
SMfig1 <- rbind(NewDat.no3, NewDat.nh4)
SMfig1$model <- c(rep(c("NO3","NH4"), each=5))
write.csv(SMfig1, "SMfig1_table.csv", row.names = FALSE)

# Area managed boreal forest w/o permafrost####
# not permission to put these files online at the moment
# For now, this will only show the approach. 
# see references for data
library(rgdal)
library(raster)
setwd("~/Downloads")
perm.shp <-readOGR(dsn=".",layer="Permafrost_circumboreal_russe_diss")
str(perm.shp)
perm.shp <- spTransform(perm.shp, CRS("+init=epsg:4326"))
boreal = subset(perm.shp, Extent_DC=="Non")
boreal@data = droplevels(boreal@data)
plot(boreal)
area(boreal)/10000/1000000 # total area boreal biome

rr = raster("Managed_Boreal_Forest_mask.tif")
plot(rr)
rr[rr==0] <- NA
m.bor <- mask(rr, boreal)
plot(m.bor)
m.bor.area <- area(m.bor, na.rm=TRUE)

#stack <- stack(m.bor, m.bor.area)
ar = mask(m.bor, m.bor.area, maskvalue=0)

ar.km2 <- cellStats(m.bor.area, sum)
ar.km2*100/1000000 # million ha


# Map over sites ####
# You can use this link for an overview map of the experimental sites
# https://drive.google.com/open?id=1YI1rB8bTh6yQxq5RTvWSdyZasrk&usp=sharing

# Check residuals ####
# here with soilC model as an example
res.test <- soilC@frame # example with flux1 model
res.test$resid <- residuals(flux1 ,type="pearson")
# numeric predictors
res.test %>%
  select_if(is.numeric) %>%
  gather(preds, value, 2:(ncol(.)-1)) %>% 
  ggplot(aes(x=value, y=resid)) +
  geom_point() +
  facet_wrap(~ preds, scales = "free")
# categorical predictors also
res.test %>%
  select(-site) %>%
  gather(preds, value, 2:(ncol(.)-1)) %>% 
  ggplot(aes(x=value, y=resid)) +
  geom_boxplot() +
  facet_wrap(~ preds, scales = "free")

# check influential points
ggplot(data.frame(lev=hatvalues(flux1),pearson=residuals(flux1,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()

#identify points with high leverage
# set a cutoff, like 0.20 for example, to identify which points that have a large impact on the model 
res.test[which(hatvalues(flux1) >= .20),]

