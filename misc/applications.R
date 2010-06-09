# Tests on real data looking for benefits of smooth terms
#
# Note: using the AIC values returned in anova.mer is conservative, 
#		since lme4 computes the marginal AIC 
# (see
# "On the Behaviour of Marginal and Conditional Akaike Information Criteria in Linear Mixed Models",
# Sonja Greven, Thomas Kneib
# http://www.bepress.com/jhubiostat/paper202/)
#
# Author: fabians
###############################################################################
library(amer)
if(require(mlmRev)){
	
	summary(Contraception)
	(fm1 <- lmer(use ~ urban + age + livch + (1|district), Contraception, family = binomial))
	(fm2 <- lmer(use ~ urban + age + livch + (urban|district), Contraception, family = binomial))	
	(fm3 <- amer(use ~ urban + bsp(age) + livch + (urban|district), Contraception, family = binomial))
	(fm4 <- amer(use ~ urban + bsp(age, by=urban) + livch + (urban|district), Contraception, family = binomial))
	anova(fm1 ,fm2, fm3, fm4)
	plotF(fm3)
	plotF(fm4)
	
	summary(Chem97)	
	(fm1 <- lmer(score ~ gcsecnt + (1 | school) + (1 | lea), Chem97))
	(fm2 <- amer(score ~ tpU(gcsecnt) + (1 | school) + (1 | lea), Chem97))
	anova(fm1, fm2)
	plotF(fm2)
	
	xyplot(height ~ age|Subject, type="b", data=Oxboys)
	(fm1 <- lmer(height ~ age + I(age^2) + I(age^3) + I(age^4) + (age + I(age^2) | Subject), data=Oxboys))
	(fm2 <- amer(height ~ tp(age) + tp(age, k=5, by=Subject, allPen=T), data=Oxboys))
	#compare AIC
	anova(fm1, fm2)
	
	summary(ScotsSec)
	(fm1 <- lmer(attain ~ sex  + (1 | primary) + (1 |second), ScotsSec))
	(fm2 <- lmer(attain ~ sex  + verbal + (1 | primary) + (1 |second), ScotsSec))
	(fm3 <- amer(attain ~ sex  + tpU(verbal) + (1 | primary) + (1 |second), ScotsSec))
	 anova(fm1, fm2, fm3)
	plotF(fm3)
	
}

if(require(SASmixed)){
 
 summary(SIMS)
 (fm1 <- lmer(Gain ~ Pretot + (Pretot | Class), SIMS))
 (fm2 <- amer(Gain ~ tp(Pretot) + (Pretot | Class), SIMS))
 plotF(fm2)
 anova(fm1, fm2)
 
 summary(Demand)
 (fm1 <- lmer(log(d) ~ log(y) + log(rd) + log(rt) + log(rs) + (1 | State) + (1 | Year), Demand))
 Demand$logD <- log(Demand$d)
 (fm2 <- amer(logD ~ tp(y) + tp(rd) + tp(rt) + tp(rs) + (1 | State) + (1 | Year), Demand))
 anova(fm1, fm2)
 plotF(fm2)
 
 summary(HR)
 (fm1 <- lmer(HR ~ Time * Drug + baseHR + (Time | Patient), HR))
 (fm2 <- amer(HR ~ Drug + tp(Time, by=Drug) + baseHR + (Time | Patient), HR))
 anova(fm1, fm2)
 
}
