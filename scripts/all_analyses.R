# Richness and composition of phyllosphere Methylobacterium communities 
# cause variation in Arabidopsis thaliana growth
# Authors of script : Jocelyn Lauzon and Jérémie Pelletier
# This script covers all the statistical analyses and figures of the experiment.

###############################################################################-


# SETTING UP ----

## load packages ----
library(tidyverse)
library(ggplot2)
library(lmerTest)
library(MuMIn)
library(emmeans)

## import data and rearrange data frames ----
biom_data_raw <- read.csv2("data/biom_data.csv", sep = ";")
biom_data_raw <- biom_data_raw %>% remove_rownames %>% column_to_rownames(var="id")
synthcomms.raw <- read.csv2("data/arabi_synthcomms.csv", sep = ";") 
synthcomms.raw <- synthcomms.raw %>% remove_rownames %>% column_to_rownames(var="id")

## remove plants with molds in biomass dataframe ----
which(biom_data_raw$molds)
length(which(biom_data_raw$molds))
biom_data <- subset(biom_data_raw, biom_data_raw$molds != TRUE)
biom_data <- subset(biom_data, select = -4) # remove "molds" column

## biom_data dataframe adjustments ----
biom_data$diversity <- as.numeric(biom_data$diversity)
biom_data$comm[biom_data$comm == "Z"] <- "control" # change letter Z to "control" for control plants
biom_data$comm <- as.factor(biom_data$comm)
biom_data$comm <- fct_inorder(biom_data$comm, ordered = NA)
biom_data$comm <- relevel(biom_data$comm, 'control') # make 'control' as first level of factor
biom_data$shelf_x <- as.factor(biom_data$shelf_x)
biom_data$shelf_y <- as.factor(biom_data$shelf_y)
class(biom_data)
summary(biom_data)
str(biom_data)

## biom_data dataframe preparation for leaf biomass analysis (without diversity 0) ----
data_lb.d <- subset(biom_data, biom_data$diversity != 0) # remove control plants
data_lb.d <- droplevels(data_lb.d)
data_lb.d <- na.omit(data_lb.d) # remove rows with NA
class(data_lb.d)
summary(data_lb.d)
str(data_lb.d)

## create data.lb.d.c dataframe (containing the 12 strains presence/absence data) ----
synthcomms <- synthcomms.raw[rownames(data_lb.d),c(1:12)] # match rows of 'synthcomms' to 'data_lb.d' and remove clades columns
synthcomms_PA <- ifelse(synthcomms!=0, 1, 0) # create object with presence/absence data instead of relative abundance
data_lb.d.c <- cbind(data_lb.d, synthcomms_PA)
data_lb.d.c  <- data_lb.d.c %>% mutate(across(c(A, B, C, D, E, F, G, H, I, J, K, L), as.factor)) # code P/A (1/0) as factor
class(data_lb.d.c)
summary(data_lb.d.c)
str(data_lb.d.c)

  
###############################################################################-


# NULL MODEL ---- 

## Model 0.1 ---- 
mod0.1 <- lm(leaf_biom ~ 1, data = data_lb.d.c)
summary(mod0.1)

## Model 0.2 ----
mod0.2 <- lmer(leaf_biom ~ 1 + (1 | shelf_x) + (1 | shelf_y), data = data_lb.d.c, REML = F)
summary(mod0.2) 

## comparing models ----
model.sel(mod0.1, mod0.2)

## statistical conclusion ----
# Adding shelf position as random effect do not make better predictions than the 
# null model. Therefore, shelf position will not be controlled for.


###############################################################################-


# H1 - GROWTH ~ DIVERSITY ---- 
# Hypothesis 1

## Model H1.1 (linear) ----
modH1.1 <- lm(leaf_biom ~ diversity, data = data_lb.d.c)
summary(modH1.1)

## Model H1.2 (polynomial) ----
modH1.2 <- lm(leaf_biom ~ poly(diversity,2,raw=TRUE), data = data_lb.d.c)
summary(modH1.2)

## comparing models ----
model.sel(mod0.1, modH1.1, modH1.2)

## postulates verification ----

### Appendix 1: Figure S1; residuals ~ predicted ----
tiff(filename = "./figures/modH1-2_res-pred_and_QQplot.tiff", width = 6, height = 3, units = "in", pointsize = 7.5, res=300)
par(mfrow=c(1,2), mar=c(5, 4, 1, 2))
plot(scale(resid(modH1.2), center = T, scale = T) ~ fitted(modH1.2),
     pch = 1, cex = 1.2, col = "black", cex.lab = 1.2,
     xlab = "Predicted values",
     ylab = "Model residuals")
abline(h=0, lty = "dashed" )
### qqplot ----
qqnorm(resid(modH1.2), main=NULL,
       pch = 1, cex = 1.2, col = "black", cex.lab = 1.2)
qqline(resid(modH1.2))
dev.off()

## statistical conclusion ----
# Model 1.2 with polynomial regression is better than Model 1.1 (non polynomial).
# Polynomial regression will then be used for the other hypothesis analysis.

## Figure 1 ----
seq.div <- with(data_lb.d.c, seq(from = 0, to = 10, by = 0.01))

prediction1.2.CI <- data.frame(predict(modH1.2, newdata=data.frame(diversity=seq.div[101:1001]),
                         interval = "confidence", level = 0.95))
prediction1.2.PI <- data.frame(predict(modH1.2, newdata=data.frame(diversity=seq.div[101:1001]), # calculated but not used
                                      interval = "prediction", level = 0.95))

tiff(filename = "./figures/biomass_by_richness.tiff", width = 6, height = 5, units = "in", pointsize = 9, res=300)
par(mfrow=c(1,1), mar=c(4.1, 4.1, 1.1, 1.1), font.lab = 2, cex = 1.3, cex.axis = 1.3) # margins order : bottom, left, top and right
plot(leaf_biom ~ diversity, data = data_lb.d.c,
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction1.2.CI$lwr, rev(prediction1.2.CI$upr)), col = rgb(0.2,0.2,0.2,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=0.8)
points(seq.div[101:1001], prediction1.2.CI$fit, type = 'l', lty = "solid", lwd = 2, col = "black")
dev.off()


###############################################################################-


# H2 - GROWTH ~ STRAINS P/A ----
# Hypothesis 2

## Model H2 exploration ----

global.model.H2 <- lm(leaf_biom ~ A + B + C + D + E + F + G + H + I + J + K + L, data = data_lb.d.c) # saturated model
options(na.action = "na.pass")
models.H2 <-  dredge(global.model.H2, beta = "none", evaluate = TRUE, rank = AICc, extra = "adjR^2")
options(na.action = "na.omit")

## comparing models ----
model.sel(models.H2)

## all models table ----
models.H2.table <- as.data.frame(model.sel(models.H2))
models.H2.table$ID <- rownames(models.H2.table)
models.H2.table <- models.H2.table[c(ncol(models.H2.table), 1:(ncol(models.H2.table) - 1))]
colnames(models.H2.table) <- 
  c("ID","Intercept","E-046","J-078","J-088","J-059","J-043","J-067",
    "J-048","J-076","E-045","J-092","E-005","J-068","adj.R2","df","LL","AICc","ΔAICc","weight")

# table output for Appendix 2
write.csv2(models.H2.table, file = "data/models.H2.csv", row.names = F)
# Convert to pdf afterwards on local computer

## averaging ----
# not used in article, but done just for data exploration
modH2.avg <- model.avg(models.H2, subset = delta <= 2, fit = T,  data = data_lb.d.c)
summary(modH2.avg)

# extract coefficients of models with delta AICc < 2 
# to build table for article
# first model from exploration
modH2.D1 <- lm(leaf_biom ~ F + I + J, data = data_lb.d.c)
summary(modH2.D1)
# second to sixth models
summary(lm(leaf_biom ~ E + F + I + J, data = data_lb.d.c))
summary(lm(leaf_biom ~ F + G + I + J, data = data_lb.d.c))
summary(lm(leaf_biom ~ F + H + I + J, data = data_lb.d.c))
summary(lm(leaf_biom ~ C + F + I + J, data = data_lb.d.c))
summary(lm(leaf_biom ~ B + F + I + J, data = data_lb.d.c))

## postulates verification ----

# prep : make dataframe and calculate residuals for modH2.D1
df.prH2.D1 <- expand.grid(F = c(0,1),I = c(0,1),J = c(0,1))
df.prH2.D1 <- as.data.frame(lapply(df.prH2.D1, as.factor))
pred.prH2.D1 <- data.frame(predicted=predict(modH2.D1, newdata = df.prH2.D1))
pred.prH2.D1 <- cbind(df.prH2.D1,pred.prH2.D1)
data_lb.d.c_resH2.D1 <- subset(data_lb.d.c, select = -c(1:2,4:10,12:13,16:20))
data_lb.d.c_resH2.D1$rownames <- rownames(data_lb.d.c_resH2.D1)
data_resH2.D1 <- merge(data_lb.d.c_resH2.D1, pred.prH2.D1, by = c('F','I','J'))
data_resH2.D1 <- data_resH2.D1 %>% arrange(rownames)
rownames(data_resH2.D1) <- data_resH2.D1$rownames
data_resH2.D1 <- subset(data_resH2.D1, select = -5)
data_resH2.D1$residuals <- data_resH2.D1$leaf_biom-data_resH2.D1$predicted
data_resH2.D1$residuals.st <- scale(data_resH2.D1$residuals, center = T, scale = T)

### Appendix 1: Figure S2; residuals ~ predicted ----
tiff(filename = "./figures/modH2.D1_res-pred_and_QQplot.tiff", width = 6, height = 3, units = "in", pointsize = 7.5, res=300)
par(mfrow=c(1,2), mar=c(5, 4, 1, 2))
plot(residuals.st ~ predicted, data = data_resH2.D1,
     pch = 1, cex = 1.2, col = "black", cex.lab = 1.2,
     xlab = "Predicted values",
     ylab = "Model residuals")
abline(h=0, lty = "dashed")
### qqplot ----
qqnorm(data_resH2.D1$residuals, main=NULL,
       pch = 1, cex = 1.2, col = "black", cex.lab = 1.2)
qqline(data_resH2.D1$residuals)
dev.off()

## statistical conclusion ----
# F and J have a positive effect on leaf biomass, while I has a negative effect

## Figure 2 ----

# Create a new dataframe with a column to represent different combinations of F, I, and J
biom_data.fig2 <- biom_data
biom_data.fig2 <- subset(biom_data.fig2, select = -c(4:5))
biom_data.fig2 <- na.omit(biom_data.fig2)
biom_data.fig2$combination <- ifelse(
  grepl("F", biom_data.fig2$comm) & grepl("I", biom_data.fig2$comm) & grepl("J", biom_data.fig2$comm), "J-067 & E-045 & J-092",
  ifelse(grepl("F", biom_data.fig2$comm) & grepl("I", biom_data.fig2$comm), "J-067 & E-045",
         ifelse(grepl("F", biom_data.fig2$comm) & grepl("J", biom_data.fig2$comm), "J-067 & J-092",
                ifelse(grepl("I", biom_data.fig2$comm) & grepl("J", biom_data.fig2$comm), "E-045 & J-092",
                       ifelse(grepl("F", biom_data.fig2$comm), "J-067",
                              ifelse(grepl("I", biom_data.fig2$comm), "E-045",
                                     ifelse(grepl("J", biom_data.fig2$comm), "J-092",
                                            ifelse(biom_data.fig2$comm == "control", "Control", "Other strains only"))))))))

# Define colors for each combination
combination_colors <- c(
  "J-067 & E-045 & J-092" = "gray27",
  "J-067 & E-045" = "purple",
  "J-067 & J-092" = "green3",
  "E-045 & J-092" = "darkorange",
  "J-067" = "blue",
  "E-045" = "red",
  "J-092" = "yellow1",
  "Other strains only" = "white",
  "Control" = "gray"
)

# Order of the legend labels
legend_order <- c(
  "J-067 & E-045 & J-092",
  "J-067 & E-045",
  "J-067 & J-092",
  "E-045 & J-092",
  "J-067",
  "E-045",
  "J-092",
  "Control",
  "Other strains only"
)

# Plot with different colors for each combination
ggplot(biom_data.fig2, aes(x = comm, y = leaf_biom, fill = combination)) +
  geom_line() +
  geom_point(shape=21, size=3) +
  geom_vline(xintercept = 13.5, color = "gray60", linetype = "dotted", size = 0.5) +
  geom_vline(xintercept = 25.5, color = "gray60", linetype = "dotted", size = 0.5) +
  geom_vline(xintercept = 31.5, color = "gray60", linetype = "dotted", size = 0.5) +
  geom_vline(xintercept = 37.5, color = "gray60", linetype = "dotted", size = 0.5) +
  geom_vline(xintercept = 43.5, color = "gray60", linetype = "dotted", size = 0.5) +
  scale_fill_manual(values = combination_colors, breaks = legend_order) +
  theme_classic() +
  xlab("Community composition") +
  ylab("Dry leaf biomass (mg)") +
  theme(text = element_text(size = 11, face = "bold")) +
  theme(axis.title.y = element_text(vjust = 3.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(fill='Strains combinations') +
  guides(fill = guide_legend(reverse = TRUE)) +
  annotate("text", x = 7, y = 60, label = "monocultures", colour = "black", size = 3.5, vjust = 0) +
  annotate("text", x = 20, y = 60, label = "2 strains", colour = "black", size = 3.5, vjust = 0) +
  annotate("text", x = 28.5, y = 60, label = "4 strains", colour = "black", size = 3.5, vjust = 0) +
  annotate("text", x = 34.5, y = 60, label = "6 strains", colour = "black", size = 3.5, vjust = 0) +
  annotate("text", x = 40.5, y = 60, label = "8 strains", colour = "black", size = 3.5, vjust = 0) +
  annotate("text", x = 46.5, y = 60, label = "10 strains", colour = "black", size = 3.5, vjust = 0)
ggsave("./figures/biomass_by_composition.tiff", path=".", width = 9, height = 5)


###############################################################################-


# H3 - GROWTH ~ DIVERSITY*STRAINS ----
# Hypothesis 3

## Model H3 exploration ----

global.model.H3 <- lm(leaf_biom ~ poly(diversity,2,raw=TRUE) * F * I * J, data = data_lb.d.c)
options(na.action = "na.pass")
models.H3 <-  dredge(global.model.H3, beta = "none", evaluate = TRUE, rank = AICc, extra = "adjR^2")
options(na.action = "na.omit")

## comparing models ----
model.sel(models.H3)

## all models table ----
models.H3.table <- as.data.frame(model.sel(models.H3))
models.H3.table$ID <- rownames(models.H3.table)
models.H3.table <- models.H3.table[c(ncol(models.H3.table), 1:(ncol(models.H3.table) - 1))]
colnames(models.H3.table) <- 
  c("ID","Intercept","J-067","E-045","J-092","D2","J-067:E-045","J-067:J-092","J-067:D2","E-045:J-092",
"E-045:D2","J-092:D2","J-067:E-045:J-092","J-067:E-045:D2","J-067:J-092:D2",
"E-045:J-092:D2","J-067:E-045:J-092:D2","adj.R2","df","LL","AICc","ΔAICc","weight")

# Output for Appendix 3
write.csv2(models.H3.table, file = "data/models.H3.csv", row.names = F)
# Convert to pdf afterwards on local computer

## averaging ----
modH3 <- model.avg(models.H3, subset = delta <= 2, fit = T,  data = data_lb.d.c)
summary(modH3)

## postulates ----

# prep : make dataframe and calculate residuals for modH3
df.pr4 <- expand.grid(diversity = seq.div, F = c(0,1), I = c(0,1), J = c(0,1))
df.pr4 <- data.frame(df.pr4)
df.pr4$F <- as.factor(df.pr4$F)
df.pr4$I <- as.factor(df.pr4$I)
df.pr4$J <- as.factor(df.pr4$J)
df.prH3 <- filter(df.pr4, diversity == 1.00 | diversity == 2.00 | diversity == 4.00 | diversity == 6.00
                 | diversity == 8.00 | diversity == 10.00)
predH3 <- data.frame(predicted=predict(modH3, newdata = df.prH3))
predH3 <- cbind(df.prH3,predH3)
predH3$diversity <- as.numeric(predH3$diversity)
predH3 <- predH3 %>% arrange(diversity)
data_lb.d.c_resH3 <- subset(data_lb.d.c, select = -c(2,4:10,12:13,16:20))
data_lb.d.c_resH3$rownames <- rownames(data_lb.d.c_resH3)
data_resH3 <- merge(data_lb.d.c_resH3, predH3, by = c("diversity", "F", "I", "J"))
data_resH3 <- data_resH3 %>% arrange(rownames)
rownames(data_resH3) <- data_resH3$rownames
data_resH3 <- subset(data_resH3, select = -6)
data_resH3$residuals <- data_resH3$leaf_biom-data_resH3$predicted
data_resH3$residuals.st <- scale(data_resH3$residuals, center = T, scale = T)

### Appendix 1: Figure S3; residuals ~ predicted ----
tiff(filename = "./figures/modH3_res-pred_and_QQplot.tiff", width = 6, height = 3, units = "in", pointsize = 7.5, res=300)
par(mfrow=c(1,2), mar=c(5, 4, 1, 2))
plot(residuals.st ~ predicted, data = data_resH3,
     pch = 1, cex = 1.2, col = "black", cex.lab = 1.2,
     xlab = "Predicted values",
     ylab = "Model residuals")
abline(h=0, lty = "dashed")
### qqplot ----
qqnorm(data_resH3$residuals, main=NULL,
       pch = 1, cex = 1.2, col = "black", cex.lab = 1.2)
qqline(data_resH3$residuals)
dev.off()


## EMM  ----

# Data exploration 
par(mfrow=c(1,2))
seq.div <- with(data_lb.d.c, seq(from = 0, to = max(diversity), by = 0.01))

em4.5 <- emmeans(modH3, data = data_lb.d.c, specs = ~ poly(diversity,2,raw=TRUE) * F * I * J, at=list(diversity=4.5))
plot(em4.5, y, type = "linear.predictor", at = list(diversity = c(1,2,4,6,8,10)))
summary(em4.5)

em6 <- emmeans(modH3, data = data_lb.d.c, specs = ~ poly(diversity,2,raw=TRUE) * F * I * J, at=list(diversity=6))
plot(em6, y, type = "linear.predictor", at = list(diversity = c(1,2,4,6,8,10)))
summary(em6)

# reference figure (for attesting validity of next hand-made figure with 8 panels)
emmip(modH3,  F * I * J ~ poly(diversity,2,raw=TRUE), data = data_lb.d.c, type = "response",  at = list(diversity = seq.div))

## predicted values of modH3 ----

# extract values of fit and 95% confidence intervals
emmip.info <- emmip(modH3,  F * I * J ~ poly(diversity,2,raw=TRUE), data = data_lb.d.c, type = "response", CIs = TRUE,  at = list(diversity = seq.div))

# Fit + CI by conditions using emmip output
modH3.FCI <- data.frame(emmip.info[["data"]])
# subset by strains presence or absence (possible combinations)
LCI_F0.I0.J0 <- filter(modH3.FCI, F == 0 & I == 0 & J == 0) 
LCI_F1.I0.J0 <- filter(modH3.FCI, F == 1 & I == 0 & J == 0) 
LCI_F0.I1.J0 <- filter(modH3.FCI, F == 0 & I == 1 & J == 0) 
LCI_F0.I0.J1 <- filter(modH3.FCI, F == 0 & I == 0 & J == 1) 
LCI_F1.I1.J0 <- filter(modH3.FCI, F == 1 & I == 1 & J == 0) 
LCI_F1.I0.J1 <- filter(modH3.FCI, F == 1 & I == 0 & J == 1)
LCI_F0.I1.J1 <- filter(modH3.FCI, F == 0 & I == 1 & J == 1) 
LCI_F1.I1.J1 <- filter(modH3.FCI, F == 1 & I == 1 & J == 1)

## Figure 3, with 8 panels ----
tiff(filename = "./figures/biomass_by_strainXrichness-panels.tiff", width = 9, height = 9, units = "in", pointsize = 5.75, res=300)
par(mfrow=c(3,3), mar=c(4.1, 4.1, 0.1, 0.1), xpd=TRUE, font.lab = 2, cex = 1.3, cex.axis = 1.3) # margins order : bottom, left, top and right

# F1-IO-J0
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:801], rev(seq.div[101:801])), c(LCI_F1.I0.J0$LCL[101:801], rev(LCI_F1.I0.J0$UCL[101:801])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[101:801], LCI_F1.I0.J0$yvar[101:801], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "J-067", bty = "n", cex = 1.2)

# F0-I1-J0
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "",
     ylab = "",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:601], rev(seq.div[101:601])), c(LCI_F0.I1.J0$LCL[101:601], rev(LCI_F0.I1.J0$UCL[101:601])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[101:601], LCI_F0.I1.J0$yvar[101:601], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "E-045", bty = "n", cex = 1.2)

# F0-I0-J1
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "",
     ylab = "",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:601], rev(seq.div[101:601])), c(LCI_F0.I0.J1$LCL[101:601], rev(LCI_F0.I0.J1$UCL[101:601])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[101:601], LCI_F0.I0.J1$yvar[101:601], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "J-092", bty = "n", cex = 1.2)

# F1-I1-J0
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[201:1001], rev(seq.div[201:1001])), c(LCI_F1.I1.J0$LCL[201:1001], rev(LCI_F1.I1.J0$UCL[201:1001])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[201:1001], LCI_F1.I1.J0$yvar[201:1001], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "J-067 & E-045", bty = "n", cex = 1.2)

# F1-I0-J1
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "",
     ylab = "",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[801:1001], rev(seq.div[801:1001])), c(LCI_F1.I0.J1$LCL[801:1001], rev(LCI_F1.I0.J1$UCL[801:1001])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[801:1001], LCI_F1.I0.J1$yvar[801:1001], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "J-067 & J-092", bty = "n", cex = 1.2)

# F0-I1-J1
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[401:1001], rev(seq.div[401:1001])), c(LCI_F0.I1.J1$LCL[401:1001], rev(LCI_F0.I1.J1$UCL[401:1001])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[401:1001], LCI_F0.I1.J1$yvar[401:1001], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "E-045 & J-092", bty = "n", cex = 1.2)

# F0-IO-J0
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:401], rev(seq.div[101:401])), c(LCI_F0.I0.J0$LCL[101:401], rev(LCI_F0.I0.J0$UCL[101:401])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[101:401], LCI_F0.I0.J0$yvar[101:401], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "Other strains only", bty = "n", cex = 1.2)

# F0-I1-J1
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[601:1001], rev(seq.div[601:1001])), c(LCI_F1.I1.J1$LCL[601:1001], rev(LCI_F1.I1.J1$UCL[601:1001])), col = rgb(0.1,1,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=1)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 1, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "forestgreen", pch = 16, cex=1.5)
points(seq.div[601:1001], LCI_F1.I1.J1$yvar[601:1001], type = 'l', lty = "solid", lwd = 2, col = "forestgreen")
legend(x = "topleft", inset=c(0,0), legend = "J-067, E-045 & J-092", bty = "n", cex = 1.2)

dev.off()
### figure 3 end


## Extra figure (not in article); all strains together to help interpretation ----
tiff(filename = "./figures/biomass_by_strainXrichness-all.tiff", width = 8.5, height = 5, units = "in", pointsize = 6.5, res=300)
par(mfrow=c(1,1), mar=c(4.1, 4.1, 1.1, 16.1), xpd=TRUE, font.lab = 2, cex = 1.3, cex.axis = 1.3) # margins order : bottom, left, top and right
plot(leaf_biom ~ diversity, data = data_lb.d.c, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "hotpink", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 0)$leaf_biom, col = "darkgrey", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 0)$leaf_biom, col = "blue", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 0)$leaf_biom, col = "red", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 0 & J == 1)$leaf_biom, col = "yellow1", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 0)$leaf_biom, col = "purple", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 0 & J == 1)$leaf_biom, col = "green3", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 0 & I == 1 & J == 1)$leaf_biom, col = "orange", pch = 16, cex=1.25)
points(filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$diversity, filter(data_lb.d.c, F == 1 & I == 1 & J == 1)$leaf_biom, col = "black", pch = 16, cex=1.25)

points(seq.div[101:401], LCI_F0.I0.J0$LCL[101:401], type = 'l', lty = "dotted", lwd = 1.25, col = "darkgrey")
points(seq.div[101:401], LCI_F0.I0.J0$UCL[101:401], type = 'l', lty = "dotted", lwd = 1.25, col = "darkgrey")

points(seq.div[101:801], LCI_F1.I0.J0$LCL[101:801], type = 'l', lty = "dotted",  lwd = 1.25, col = "blue")
points(seq.div[101:801], LCI_F1.I0.J0$UCL[101:801], type = 'l', lty = "dotted",  lwd = 1.25, col = "blue")

points(seq.div[101:601], LCI_F0.I1.J0$LCL[101:601], type = 'l', lty = "dotted", lwd = 1.25, col = "red")
points(seq.div[101:601], LCI_F0.I1.J0$UCL[101:601], type = 'l', lty = "dotted", lwd = 1.25, col = "red")

points(seq.div[101:601], LCI_F0.I0.J1$LCL[101:601], type = 'l', lty = "dotted", lwd = 1.25, col = "yellow1")
points(seq.div[101:601], LCI_F0.I0.J1$UCL[101:601], type = 'l', lty = "dotted", lwd = 1.25, col = "yellow1")

points(seq.div[201:1001], LCI_F1.I1.J0$LCL[201:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "purple")
points(seq.div[201:1001], LCI_F1.I1.J0$UCL[201:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "purple")

points(seq.div[801:1001], LCI_F1.I0.J1$LCL[801:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "green3")
points(seq.div[801:1001], LCI_F1.I0.J1$UCL[801:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "green3")

points(seq.div[401:1001], LCI_F0.I1.J1$LCL[401:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "darkorange")
points(seq.div[401:1001], LCI_F0.I1.J1$UCL[401:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "darkorange")

points(seq.div[601:1001], LCI_F1.I1.J1$LCL[601:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "black")
points(seq.div[601:1001], LCI_F1.I1.J1$UCL[601:1001], type = 'l', lty = "dotted", lwd = 1.25, col = "black")

points(seq.div[101:401], LCI_F0.I0.J0$yvar[101:401], type = 'l', lty = "solid", lwd = 2, col = "darkgrey")
points(seq.div[101:801], LCI_F1.I0.J0$yvar[101:801], type = 'l', lty = "solid", lwd = 2, col = "blue")
points(seq.div[101:601], LCI_F0.I1.J0$yvar[101:601], type = 'l', lty = "solid", lwd = 2, col = "red")
points(seq.div[101:601], LCI_F0.I0.J1$yvar[101:601], type = 'l', lty = "solid", lwd = 2, col = "yellow1")
points(seq.div[201:1001], LCI_F1.I1.J0$yvar[201:1001], type = 'l', lty = "solid", lwd = 2, col = "purple")
points(seq.div[801:1001], LCI_F1.I0.J1$yvar[801:1001], type = 'l', lty = "solid", lwd = 2, col = "green3")
points(seq.div[401:1001], LCI_F0.I1.J1$yvar[401:1001], type = 'l', lty = "solid", lwd = 2, col = "darkorange")
points(seq.div[601:1001], LCI_F1.I1.J1$yvar[601:1001], type = 'l', lty = "solid", lwd = 2, col = "black")

legend(x = "topright",
       inset=c(-0.38,0.3),
       legend = c("control plants","other strains only","J-067","E-045","J-092","J-067 & E-045","J-067 & J-092","E-045 & J-092","J-067, E-045 & J-092"), # Vector with the name of each group
       col = c("hotpink","darkgrey","blue","red","yellow1","purple","green3","darkorange","black"),
       lty = "solid", lwd = 2,
       pch = 16,
       bty = "n",
       bg = par("bg"),
       cex = 1.2,
       horiz = FALSE, 
       title = NULL)
dev.off()
# extra figure end


###############################################################################-

# Complimentarity vs selection ----

## create and rearrange new dataframes ----
data.CS <- data_lb.d.c
data.CS <- data.CS[,-c(4:20)]
synthcomms.CS <- synthcomms_PA

# Calculate the observed means of communities
means <- data.CS %>%
  group_by(comm) %>%
  summarise(obs_avg = mean(leaf_biom))

data.CS <- data.CS %>%
  left_join(means, by = "comm")

# Calculate the observed standard deviation of communities
sds <- data.CS %>%
  group_by(comm) %>%
  summarise(obs_sd = sd(leaf_biom))

data.CS <- data.CS %>%
  left_join(sds, by = "comm")

# Values of observed means for 12 strains in monocultures
monocult_means <- unique(data.CS[c(1:34),c(2,4)])
monocult_sd <- unique(data.CS[c(1:34),c(2,5)])

# In synthcomms community matrix, replace 1's with strains monoculture mean
synthcomms.CS2 <- data.frame(A = ifelse(synthcomms.CS[,"A"], monocult_means[1,2],NA),
                          B = ifelse(synthcomms.CS[,"B"], monocult_means[2,2],NA),
                          C = ifelse(synthcomms.CS[,"C"], monocult_means[3,2],NA),
                          D = ifelse(synthcomms.CS[,"D"], monocult_means[4,2],NA),
                          E = ifelse(synthcomms.CS[,"E"], monocult_means[5,2],NA),
                          F = ifelse(synthcomms.CS[,"F"], monocult_means[6,2],NA),
                          G = ifelse(synthcomms.CS[,"G"], monocult_means[7,2],NA),
                          H = ifelse(synthcomms.CS[,"H"], monocult_means[8,2],NA),
                          I = ifelse(synthcomms.CS[,"I"], monocult_means[9,2],NA),
                          J = ifelse(synthcomms.CS[,"J"], monocult_means[10,2],NA),
                          K = ifelse(synthcomms.CS[,"K"], monocult_means[11,2],NA),
                          L = ifelse(synthcomms.CS[,"L"], monocult_means[12,2],NA))
synthcomms.CS3 <- synthcomms.CS2
synthcomms.CS4 <- synthcomms.CS2

###############################################################################-

## Appendix 1: Figure S4; - framework curves ----

x <- seq(0, 10, length.out = 1000)
y <- log10(3*x + 1) / log10(25)
y0 <- seq(0, 0, length.out = 1000)

tiff(filename = "./figures/comp_vs_selection_framework.tiff", width = 6, height = 5, units = "in", pointsize = 7, res=300)
par(mfrow=c(1,1), mar=c(4.1, 4.1, 1.1, 1.1), font.lab = 2, cex = 1.3, cex.axis = 1.3)
plot(1, type="n", data = data.CS,
     pch = 5, cex = 1, col = "blue", cex.lab = 1.3,
     xlab = "Diversity",
     ylab = "Productivity",
     ylim = c(-1.5,1.5),
     xlim = c(0,10))
points(x, y, type = 'l', lty = "solid", lwd = 1, col = "red")
points(x, -y, type = 'l', lty = "solid", lwd = 1, col = "orange")
points(x, y0, type = 'l', lty = "solid", lwd = 1, col = "blue")
text(2.5, 1.23, "positive \n complementarity", cex = 1.3)
text(2.5, -1.1, "negative \n complementarity", cex = 1.3)
text(9, 0.1, "no NBE", cex = 1.3, col = "blue")
text(9, 0.8, "positive \n selection", cex = 1.3, col = "red")
text(9, -0.8, "negative\n selection", cex = 1.3, col = "orange")
text(5.5, 0.5, "indistinct \n positive NBE", cex = 1.3)
text(5.5, -0.5, "indistinct \n negative NBE", cex = 1.3)
dev.off()


###############################################################################-


## Alternative models ----

### Predicted values ----

# Calculate the predicted means of communities, based on the mean of their strains in monocultures
synthcomms.CS2$pred_avg <- rowMeans(synthcomms.CS2[,c(1:12)], na.rm = T)

# Determine the predicted max value of communities, based on the highest mean value of one of its strain in monoculture
row_max <- apply(synthcomms.CS3, 1, max, na.rm = T)
synthcomms.CS2$pred_max <- row_max
rm(synthcomms.CS3)

# Determine the predicted min value of communities, based on the lowest mean value of one of its strain in monoculture
row_min <- apply(synthcomms.CS4, 1, min, na.rm = T)
synthcomms.CS2$pred_min <- row_min
rm(synthcomms.CS4)

### New dataframe for analyses ---- 
data.CS2 <- cbind(data.CS, pred_avg=synthcomms.CS2$pred_avg, pred_max=synthcomms.CS2$pred_max, pred_min=synthcomms.CS2$pred_min)

### Models ----

### Positive selection model based on maximum values predicted ----
lm1.0 <- lm(pred_max ~ 1, data = data.CS2)
lm1.1 <- lm(pred_max ~ diversity, data = data.CS2)
lm1.2 <- lm(pred_max ~ poly(diversity,2,raw=TRUE), data = data.CS2)
model.sel(lm1.0, lm1.1, lm1.2)
# Polynomial regression (lm1.2) is the best model
mod_pos.sel <- lm1.2
summary(mod_pos.sel)
# model predictions
seq.div <- with(data.CS2, seq(from = 0, to = 10, by = 0.01))
prediction.mod_pos.sel <- data.frame(predict(mod_pos.sel, newdata=data.frame(diversity=seq.div[101:1001]),
                                         interval = "confidence", level = 0.95))

### No diversity effect model based on average values predicted ----
lm2.0 <- lm(pred_avg ~ 1, data = data.CS2)
lm2.1 <- lm(pred_avg ~ diversity, data = data.CS2)
lm2.2 <- lm(pred_avg ~ poly(diversity,2,raw=TRUE), data = data.CS2)
model.sel(lm2.0, lm2.1, lm2.2)
# Null model (lm2.0) is marginal a better model than lm2.1, but we'll use lm2.1 
# since it is more representative of the alternative model we want to use
mod_no.div.fx <- lm2.1
summary(mod_no.div.fx)
# model predictions
prediction.mod_no.div.fx <- data.frame(predict(mod_no.div.fx, newdata=data.frame(diversity=seq.div[101:1001]),
                                         interval = "confidence", level = 0.95))

### Negative selection model based on minimum values predicted ----
lm3.0 <- lm(pred_min ~ 1, data = data.CS2)
lm3.1 <- lm(pred_min ~ diversity, data = data.CS2)
lm3.2 <- lm(pred_min ~ poly(diversity,2,raw=TRUE), data = data.CS2)
model.sel(lm3.0, lm3.1, lm3.2)
# Polynomial regression (lm3.2) is the best model
mod_neg.sel <- lm3.2
summary(mod_neg.sel)
# model predictions
prediction.mod_neg.sel <- data.frame(predict(mod_neg.sel, newdata=data.frame(diversity=seq.div[101:1001]),
                                         interval = "confidence", level = 0.95))

### observed data model : Model H1.2 (polynomial) ----
mod_obs <- lm(leaf_biom ~ poly(diversity,2,raw=TRUE), data = data.CS2)
# model predictions
prediction.mod_obs <- data.frame(predict(mod_obs, newdata=data.frame(diversity=seq.div[101:1001]),
                                       interval = "confidence", level = 0.95))


###############################################################################-

## Figures ----

### Appendix 1: Figure S5; panels of the four models ----

tiff(filename = "./figures/alternative_and_obs_models.tiff", width = 6, height = 6, units = "in", pointsize = 4, res=300)
par(mfrow=c(2,2), mar=c(4.1, 4.1, 1.1, 1.1), font.lab = 2, cex = 1.3, cex.axis = 1.3) # margins order : bottom, left, top and right

plot(pred_avg ~ diversity, data = data.CS2,
     pch = 5, cex = 1, col = "blue", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(0, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_no.div.fx$upr, rev(prediction.mod_no.div.fx$lwr)), col = rgb(0,0,1,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=0.8)
points(seq.div[101:1001], prediction.mod_no.div.fx$fit, type = 'l', lty = "solid", lwd = 1, col = "blue")
text(0, 64, "a", cex=1.8)

plot(pred_max ~ diversity, data = data.CS2,
     pch = 2, cex = 1, col = "red", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(0, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_pos.sel$upr, rev(prediction.mod_pos.sel$lwr)), col = rgb(1,0,0,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=0.8)
points(seq.div[101:1001], prediction.mod_pos.sel$fit, type = 'l', lty = "solid", lwd = 1, col = "red")
text(0, 64, "b", cex=1.8)

plot(pred_min ~ diversity, data = data.CS2,
     pch = 6, cex = 1, col = "orange", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(0, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_neg.sel$upr, rev(prediction.mod_neg.sel$lwr)), col = rgb(1,0.8,0,0.3), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=0.8)
points(seq.div[101:1001], prediction.mod_neg.sel$fit, type = 'l', lty = "solid", lwd = 1, col = "orange")
text(0, 64, "c", cex=1.8)

plot(leaf_biom ~ diversity, data = data.CS2,
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(0, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_obs$upr, rev(prediction.mod_obs$lwr)), col = rgb(0.2,0.2,0.2,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=0.8)
points(seq.div[101:1001], prediction.mod_obs$fit, type = 'l', lty = "solid", lwd = 1, col = "black")
text(0, 64, "d", cex=1.8)

dev.off()
# figure S5 end


### Figure 4 - Observed data with all 4 models predictions ----

tiff(filename = "./figures/comp_vs_select_predictions.tiff", width = 6, height = 5, units = "in", pointsize = 9, res=300)
par(mfrow=c(1,1), mar=c(4.1, 4.1, 1.1, 1.1), font.lab = 2, cex = 1.3, cex.axis = 1.3) # margins order : bottom, left, top and right
plot(pred_avg ~ diversity, data = data.CS2, type = "n",
     pch = 1, cex = 1, col = "black", cex.lab = 1.3,
     xlab = "Strain richness",
     ylab = "Dry leaf biomass (mg)",
     ylim = c(5, 65),
     xlim = c(-0.25,10.5))
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_no.div.fx$upr, rev(prediction.mod_no.div.fx$lwr)), col = rgb(0,0,1,0.2), border = NA)
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_pos.sel$upr, rev(prediction.mod_pos.sel$lwr)), col = rgb(1,0,0,0.2), border = NA)
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_neg.sel$upr, rev(prediction.mod_neg.sel$lwr)), col = rgb(1,0.8,0,0.3), border = NA)
polygon(c(seq.div[101:1001], rev(seq.div[101:1001])), c(prediction.mod_obs$upr, rev(prediction.mod_obs$lwr)), col = rgb(0.2,0.2,0.2,0.2), border = NA)
points(filter(biom_data, diversity == 0)$diversity, filter(biom_data, diversity == 0)$leaf_biom, col = "black", pch = 4, cex=0.8)
points(data.CS2$diversity, data.CS2$leaf_biom, col = "black", pch = 1, cex=1)
points(seq.div[101:1001], prediction.mod_no.div.fx$fit, type = 'l', lty = "solid", lwd = 2, col = "blue")
points(seq.div[101:1001], prediction.mod_pos.sel$fit, type = 'l', lty = "solid", lwd = 2, col = "red")
points(seq.div[101:1001], prediction.mod_neg.sel$fit, type = 'l', lty = "solid", lwd = 2, col = "orange")
points(seq.div[101:1001], prediction.mod_obs$fit, type = 'l', lty = "solid", lwd = 2, col = "black")
dev.off()


### SCRIPT END ###