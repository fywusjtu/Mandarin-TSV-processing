## JML 2020, setwd("~/Fuyun_Documents/Topic_processing/github/Github_reanalysis_201911/2020_reanalysis/JML_2020_report")
## topic Exp. 1, 52id

library(readr)
library(ggplot2)
library(dplyr)
library(lme4)
library(MASS)
library(lattice)

topic <- read.csv("topic_e1_52id.csv") ## 8528 rows of datapoints

##target data without Qs
topic.noQ <- subset(topic, topic$region != "?") ##7696 rows of datapoints
xtabs(~condition + item, data = topic.noQ)
xtabs(~condition + subj, data = topic.noQ)
topic = topic.noQ

topic$subj <- as.factor(topic$subj)
topic$item <- as.factor(topic$item)
str(topic)
nrow(topic)

boxcox(topic$RT~ topic$condition * topic$subj) ##lamda peaks near -1, thus use rrt
topic$rrt <- -1000/topic$RT
head(topic, n=5)

topic.trim <- subset(topic, topic$RT > 100 & topic$ RT < 3000) ## eliminated 13 (0.17%) datapoinjts from the original 7696 datapoints nrow(topic.trim)
qqmath(~rrt|subj, data=topic.trim)
qqmath(~rrt|item, data=topic.trim)

write.csv(topic.trim, file = "topic1_52id_trim.csv")

## define sum contrast coding
contrasts(topic.trim$order) <- contr.sum(2)

### at W5: region =4 (MV vs. agent-NP)
W5 <- subset(topic.trim, region == "4")
## add RTs of previous region
W4 <- subset(topic.trim, region =="3")
df_prevW5 <- dplyr::select(W4, subj, item, rrt)
df_prevW5 <- dplyr::rename(df_prevW5, prev_rrt = rrt)
W5 <- left_join(W5, df_prevW5)
head(W5, n=5)
str(W5)

m.W5 <-lmer(rrt~order + (1+order|subj) + (1+order|item) + prev_rrt+ logW + Wlength + Nstrokes, 
             W5) ## fail to converge
m1.W5 <-lmer(rrt~order + (1+order|subj) + (1|item) + prev_rrt+ logW + Wlength + Nstrokes, 
            W5) 
summary(m1.W5)  ## SVO faster than  TSV (t =-2.48, p=0.016); ME of previous rrt (t=16.03)  

## at W6: region = 5 (all vs. this)
W6 <- subset(topic.trim, region == "5")
## add RTs of previous region
df_prevW6 <- dplyr::select(W5, subj, item, rrt)
df_prevW6 <- dplyr::rename(df_prevW6, prev_rrt = rrt)
W6 <- left_join(W6, df_prevW6)
head(W6, n=5)
str(W6)

m.W6 <-lmer(rrt~order + (1+order|subj) + (1+order|item) + prev_rrt+ logW + Wlength + Nstrokes, 
           W6)  ## drop Wlength
m.W6 <-lmer(rrt~order + (1+order|subj) + (1+order|item) + prev_rrt+ logW + Nstrokes, 
            W6)
summary(m.W6) ## no diff. betw. TSV and SVO; only ME of previous rrt (t=11.37)

## at W7: region =6 (MV vs. NP2)
W7 <- subset(topic.trim, region == "6")
## add RTs of previous region
df_prevW7 <- dplyr::select(W6, subj, item, rrt)
df_prevW7 <- dplyr::rename(df_prevW7, prev_rrt = rrt)
W7 <- left_join(W7, df_prevW7)
head(W7, n=5)
str(W7)

m.W7 <-lmer(rrt~order + (1+order|subj) + (1+order|item) + prev_rrt+ logW + Wlength + Nstrokes, 
            W7)
summary(m.W7) ## TSV faster than SVO (t=2.19, p=0.033; SVO slower than TSV at MV); ME of previous rrt (t=17.977)

## at W8: region =7 (conjution 'and')
W8 <- subset(topic.trim, region == "7")
## add RTs of previous region
df_prevW8 <- dplyr::select(W7, subj, item, rrt)
df_prevW8 <- dplyr::rename(df_prevW8, prev_rrt = rrt)
W8 <- left_join(W8, df_prevW8)
head(W8, n=5)
str(W8)

m.W8 <-lmer(rrt~order + (1+order|subj) + (1+order|item) + prev_rrt+ logW + Wlength + Nstrokes, 
            W8)
summary(m.W8) ## no diff.

## exp. 2, 71id
topic2 <- read.csv("topic2_71id_reanalysis.csv")  ##20448 rows of datapoints

### target sentences
topic2.noQ <- subset(topic2, topic2$position != "?") ##18744 rows of datapoints
xtabs(~condition + item, data = topic2.noQ)
xtabs(~condition + subj, data = topic2.noQ)
topic2 = topic2.noQ

topic2$subj <- as.factor(topic2$subj)
topic2$item <- as.factor(topic2$item)
topic2$constituent <- as.factor(topic2$constituent)

boxcox <- boxcox(topic2$RT~ topic2$condition * topic2$subj) ##lamda peaks around -1, thus use rrt to normalize the distribution of RT
topic2$rrt <- -1000/topic2$RT

topic2.trim <-subset(topic2, topic2$RT > 100 & topic2$ RT < 3000)  ## remove 29 (0.15%) outliers
qqmath(~rrt|subj, data=topic2.trim)  ## better use trimmed data, rrts for both subj and item are linear 
qqmath(~rrt|item, data=topic2.trim)


## define sum contrast coding
contrasts(topic2.trim$order) <- contr.sum(2)
contrasts(topic2.trim$clauseType) <- contr.sum(2)

str(topic2.trim)

write.csv (topic2.trim, file = "topic2_71id_trim_for_plot.csv")

## first, regress on extraneous predictors, including log frequency, number of strokes, position
l <- lmer(rrt ~ logW + Nstrokes + Wlength  + position +
            (1 | subj) + (1 | item), topic2.trim)
topic2.trim$res <- residuals(l) 

write.csv(topic2.trim, file = "topic2_71id_trim_residualRRT_forPlot.csv") ## for residual plot

## second, for each word region (only look at W5 'RC-V', W7 'head noun', W11 'MV'), build linear mixed effects model, using fixed effects and previous rrt

## W5, pos = 4
W5 <- subset(topic2.trim, position =="4")
# 前一个region的rrt
W4 <- subset(topic2.trim, position =="3")
df_prevreg <- dplyr::select(W4, subj, item, rrt)
df_prevreg <- dplyr::rename(df_prevreg, prev_rrt = rrt)
W5 <- left_join(W5, df_prevreg)
head(W5, n=5)

# 再用residual做回归
model.W5 <- lmer(res ~ order*clauseType +
                   (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W5) ## fail to converge

model.W5.1 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order+clauseType | item)+ prev_rrt, W5) ## fail to converge

model.W5.2 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W5)

model.W5.3 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order+ clauseType | item) + prev_rrt, W5)

model.W5.4 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 | item) + prev_rrt, W5)

model.W5.5 <- lmer(res ~ order*clauseType +
                     (1 + order+clauseType | subj) + (1 | item) + prev_rrt, W5)
summary(model.W5.5) ## ME of word order (t= -5.08), ME of prev rrt (t=20.14); no other effects

## W6 (DE), pos = 5
W6 <- subset(topic2.trim, position =="5")
# 前一个region的rrt
df_prevreg <- dplyr::select(W5, subj, item, rrt)
df_prevreg <- dplyr::rename(df_prevreg, prev_rrt = rrt)
W6 <- left_join(W6, df_prevreg)
head(W6, n=5)
# 再用residual做回归
model.W6 <- lmer(res ~ order*clauseType +
                   (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W6) 

summary(model.W6) ## ME of word order (t=-5.758), ME of prev_rrt (t=14.046)


## W7 (HN **critical region of NP2), pos = 6
W7 <- subset(topic2.trim, position =="6")
# 前一个region的rrt
df_prevreg7 <- dplyr::select(W6, subj, item, rrt)
df_prevreg7 <- dplyr::rename(df_prevreg7, prev_rrt = rrt)
W7 <- left_join(W7, df_prevreg7)
head(W7, n=5)

# 再用residual做回归
model.W7 <- lmer(res ~ order*clauseType +
                   (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W7) ##fail to converge

model.W7.1 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order+clauseType | item) + prev_rrt, W7)

model.W7.2 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order | item) + prev_rrt, W7)

model.W7.3 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 | item) + prev_rrt, W7)
model.W7.4 <- lmer(res ~ order*clauseType +
                     (1 + order+clauseType | subj) + (1 | item) + prev_rrt, W7)
model.W7.5 <- lmer(res ~ order*clauseType +
                     (1 + order+clauseType | subj) + prev_rrt, W7)
summary(model.W7.5)  ## ME of order (t=-7.11), ME of previous rrt (t = 20.59)

## W8 (ususally, HN + 1; pos = 7)
W8 <- subset(topic2.trim, position =="7")
# 前一个region的rrt
W7 <- subset(topic2.trim, position =="6")
df_prevreg8 <- dplyr::select(W7, subj, item, rrt)
df_prevreg8 <- dplyr::rename(df_prevreg8, prev_rrt = rrt)
W8 <- left_join(W8, df_prevreg8)
head(W8, n=5)

# 再用residual做回归
model.W8 <- lmer(res ~ order*clauseType +
                   (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W8) ##fail to converge

model.W8.1 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order+clauseType | item) + prev_rrt, W8)

model.W8.2 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order | item) + prev_rrt, W8)


summary(model.W8.2) ##ME of word order, t=-2.845; ME of clausetype, t=-2.265; prev_rrt, t=19.895)

## W9 (all, pos=8)
W9 <- subset(topic2.trim, position =="8")
# 前一个region的rrt
W8 <- subset(topic2.trim, position =="7")
df_prevreg9 <- dplyr::select(W8, subj, item, rrt)
df_prevreg9 <- dplyr::rename(df_prevreg9, prev_rrt = rrt)
W9 <- left_join(W9, df_prevreg9)
head(W9, n=5)

# 再用residual做回归
model.W9 <- lmer(res ~ order*clauseType +
                   (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W9) ##fail to converge

model.W9.1 <- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order+clauseType | item) + prev_rrt, W9)

summary(model.W9.1) ##no effects; only prev_rrt: t=19.544

## W10 (MV1)
W10 <- subset(topic2.trim, position =="9")
# 前一个region的rrt
W9 <- subset(topic2.trim, position =="8")
df_prevreg10 <- dplyr::select(W9, subj, item, rrt)
df_prevreg10 <- dplyr::rename(df_prevreg10, prev_rrt = rrt)
W10 <- left_join(W10, df_prevreg10)
head(W10, n=5)

# 再用residual做回归
model.W10 <- lmer(res ~ order*clauseType +
                    (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W10) 

summary(model.W10)  ## no effects of fixed factors; only ME of previous rrt (t= 17.49)

##W11 (MV2)
W11 <- subset(topic2.trim, position =="10")
# 前一个region的rrt
df_prevreg11 <- dplyr::select(W10, subj, item, rrt)
df_prevreg11 <- dplyr::rename(df_prevreg11, prev_rrt = rrt)
W11 <- left_join(W11, df_prevreg11)
head(W11, n=5)

# 再用residual做回归
model.W11 <- lmer(res ~ order*clauseType +
                    (1 + order*clauseType | subj) + (1 + order*clauseType | item) + prev_rrt, W11)

model.W11.1<- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order+clauseType | item) + prev_rrt, W11)

model.W11.2<- lmer(res ~ order*clauseType +
                     (1 + order*clauseType | subj) + (1 + order | item) + prev_rrt, W11)

summary(model.W11.2) ## ME of order (t=-2.43), ME of complexity (t=-2.85), sig. interaction (t= -2.55); ME of previous rrt (t=16.002)

## unpack this interaction
W11$WO = as.numeric(ifelse((W11$condition=="a" | W11$condition == "c"),-0.5,0.5))
W11$complexity = as.numeric(ifelse((W11$condition=="a" | W11$condition == "b"),-0.5,0.5))

str(W11)

W11.SVO = W11[which(W11$WO == 0.5),]
W11.TSV = W11[which(W11$WO == -0.5),]
model.W11.interaction = lmer(res ~ complexity + (1+complexity | subj) + (1+complexity | item) + prev_rrt, data = W11.SVO)
summary(model.W11.interaction) ## ME of complexity for the two SVO conditions (a vs. c, t= -3.135); ME of prev_rrt (t=8.641)

model.W11.interTSV = lmer(res ~ complexity + (1+complexity | subj) + (1+complexity | item) + prev_rrt, data = W11.TSV)
summary(model.W11.interTSV)  ## no effect of complexity for the two TSV conditions (b vs. d, t=-0.156); ME of pre_rrt (t = 12.73)

W11.RC = W11[which(W11$complexity == -0.5),]
W11.MC = W11[which(W11$complexity == 0.5),]
model.W11.interRC = lmer(res ~ WO + (1+WO | subj) + (1+WO | item) + prev_rrt, data = W11.RC)
summary(model.W11.interRC)  ## no effect of word order for the two RC conditions (a vs. b, t=-0.064)

model.W11.interMC = lmer(res ~ WO + (1+WO | subj) + (1+WO | item) + prev_rrt, data = W11.MC)
summary(model.W11.interMC)  ## ME of word order for the two MC conditions (c vs. d: t=-3.331)

### exp. 3, 53id
topic4 <- read.csv("topic_e3_3conditions_53id.csv")  ##19080 datapoints

## target sentences
topic4.noQ <- subset(topic4, topic4$position != "?") ##17808 rows of datapoints
xtabs(~condition + item, data = topic4.noQ)
xtabs(~condition + subj, data = topic4.noQ)
topic4 = topic4.noQ

topic4$subj <- as.factor(topic4$subj)
topic4$item <- as.factor(topic4$item)


boxcox <- boxcox(topic4$RT~ topic4$condition * topic4$subj) ##lamda peaks around -1, thus use rrt to normalize the distribution of RT
topic4$rrt <- -1000/topic4$RT

topic4.trim <-subset(topic4, topic4$RT > 100 & topic4$ RT < 3000)  ## remove 38 (0.21%) outliers; 11770 datapoints
qqmath(~rrt|subj, data=topic4.trim)  
qqmath(~rrt|item, data=topic4.trim)

topic4= topic4.trim 

write.csv(topic4.trim, file = "topic_new_e4_53id_trim.csv")

## first, regress on extraneous predictors, including log frequency, number of strokes, position
l.topic4 <- lmer(rrt ~ logW + Nstrokes + Wlength  + position +
                   (1 | subj) + (1 | item), topic4)
topic4$res <- residuals(l.topic4)  

write.csv(topic4, file = "topic_new_e4_53id_trim_residuals_forPlot.csv")

## second, for each word region (only look at W5 'RC-V', W7 'head noun', W11 'MV'), build linear mixed effects model, using fixed effects and previous rrt

## W4 (V/adj), pos = 3
W4 <- subset(topic4, position =="3")
# 前一个region的rrt
W3 <- subset(topic4, position =="2")
df_prevreg <- dplyr::select(W3, subj, item, rrt)
df_prevreg <- dplyr::rename(df_prevreg, prev_rrt = rrt)
W4 <- left_join(W4, df_prevreg)
head(W4, n=5)

# 再用residual做回归
model.W4 <- lmer(res ~ condition +
                   (1 + condition| subj) + (1 + condition| item) + prev_rrt, W4) 
summary(model.W4)  ## no differences


## W5 (RC-V/N), pos = 4
W5 <- subset(topic4, position =="4")
# 前一个region的rrt
W4 <- subset(topic4, position =="3")
df_prevreg <- dplyr::select(W4, subj, item, rrt)
df_prevreg <- dplyr::rename(df_prevreg, prev_rrt = rrt)
W5 <- left_join(W5, df_prevreg)
head(W5, n=5)

# 再用residual做回归
model.W5 <- lmer(res ~ condition +
                   (1 + condition| subj) + (1 + condition| item) + prev_rrt, W5) 

summary(model.W5) ## no differences betw. a and b or a and c, ME of prev rrt (t=19.441)

## W6 (DE), pos = 5
W6 <- subset(topic4, position =="5")
# 前一个region的rrt
W5 <- subset(topic4, position =="4")
df_prevreg6 <- dplyr::select(W5, subj, item, rrt)
df_prevreg6 <- dplyr::rename(df_prevreg6, prev_rrt = rrt)
W6 <- left_join(W6, df_prevreg6)
head(W6, n=5)

# 再用residual做回归
model.W6 <- lmer(res ~ condition +
                   (1 + condition| subj) + (1 + condition| item) + prev_rrt, W6) 

summary(model.W6)  ## b faster than a (t = -2.47, p=0.021); no diff. betw. a and c. (t=1.16, p=0.26)


## W7 (HN), pos = 6
W7 <- subset(topic4, position =="6")
# 前一个region的rrt
W6 <- subset(topic4, position =="5")
df_prevreg7 <- dplyr::select(W6, subj, item, rrt)
df_prevreg7 <- dplyr::rename(df_prevreg7, prev_rrt = rrt)
W7 <- left_join(W7, df_prevreg7)
head(W7, n=5)

# 再用residual做回归
model.W7 <- lmer(res ~ condition+
                   (1 + condition| subj) + (1 + condition| item) + prev_rrt, W7) 


summary(model.W7)  ## b faster than a, t = -2.40 (p=0.03), ME of previous rrt (t = 21.58)

## W8 (ADV, spillover region 1), pos = 7
W8 <- subset(topic4, position =="7")
# 前一个region的rrt
W7 <- subset(topic4, position =="6")
df_prevreg8 <- dplyr::select(W7, subj, item, rrt)
df_prevreg8 <- dplyr::rename(df_prevreg8, prev_rrt = rrt)
W8 <- left_join(W8, df_prevreg8)
head(W8, n=5)
# 再用residual做回归
model.W8 <- lmer(res ~ condition+
                   (1 + condition| subj) + (1 + condition| item) + prev_rrt, W8) 
model.W8.1 <- lmer(res ~ condition+
                     (1 + condition| subj) + (1 | item) + prev_rrt, W8) 
summary(model.W8.1) ## in the spillover region of ADV, b faster than a (t=-2.80, p=0.006); no diff. betw. a and c (t=1.69, p=0.09)

model.W8.2 <-lmer(res ~ condition+
                    (1 + condition| subj)  + prev_rrt, W8) 
summary(model.W8.2)

## W9 (all, spillover region 2), pos = 8
W9 <- subset(topic4, position =="8")
# 前一个region的rrt
W8 <- subset(topic4, position =="7")
df_prevreg9 <- dplyr::select(W8, subj, item, rrt)
df_prevreg9 <- dplyr::rename(df_prevreg9, prev_rrt = rrt)
W9 <- left_join(W9, df_prevreg9)
head(W9, n=5)
# 再用residual做回归
model.W9 <- lmer(res ~ condition+
                   (1 + condition| subj) + (1 + condition| item) + prev_rrt, W9) 
model.W9.1 <- lmer(res ~ condition+
                     (1 + condition| subj) + (1 | item) + prev_rrt, W9) 
summary(model.W9.1)  ## no differences (a vs. b: t=-0.26, t=0.795; a vs. c: t = -1.12, p=0.27)
model.W9.2 <- lmer(res ~ condition+
                     (1 + condition| subj) + prev_rrt, W9)
summary(model.W9.2)

## W10 (MV1)
W10 <- subset(topic4, position =="9")
# 前一个region的rrt
W9 <- subset(topic4, position =="8")
df_prevreg10 <- dplyr::select(W9, subj, item, rrt)
df_prevreg10 <- dplyr::rename(df_prevreg10, prev_rrt = rrt)
W10 <- left_join(W10, df_prevreg10)
head(W10, n=5)

# 再用residual做回归
model.W10 <- lmer(res ~ condition+
                    (1 + condition| subj) + (1 + condition| item) + prev_rrt, W10) 

summary(model.W10)  ## no differences; only ME of previous rrt (t= 20.82)

##W11 (MV2)
W11 <- subset(topic4, position =="10")
# 前一个region的rrt
df_prevreg11 <- dplyr::select(W10, subj, item, rrt)
df_prevreg11 <- dplyr::rename(df_prevreg11, prev_rrt = rrt)
W11 <- left_join(W11, df_prevreg11)
head(W11, n=5)

# 再用residual做回归
model.W11 <- lmer(res ~ condition+
                    (1 + condition| subj) + (1 + condition| item) + prev_rrt, W11)

summary(model.W11) ## no differences; ME of previous rrt (t=16.35)

##W12 (and)
W12 <- subset(topic4, position =="11")
# 前一个region的rrt
df_prevreg12 <- dplyr::select(W11, subj, item, rrt)
df_prevreg12 <- dplyr::rename(df_prevreg12, prev_rrt = rrt)
W12 <- left_join(W12, df_prevreg12)
head(W11, n=5)

# 再用residual做回归
model.W12 <- lmer(res ~ condition+
                    (1 + condition| subj) + (1 + condition| item) + prev_rrt, W12)

summary(model.W12)
