


#---------------------------------
#Chargement objets et packages
#---------------------------------

#install.packages("glmnet")
library(glmnet)
library(boot)
library(MASS)

final <- read.csv2("data/aorte.csv")

#---------------------
#evenements tardifs
tardif <- final[final$CJPhosp != 1, ]
table(tardif$CJPtotal)

#evenements précoces
prec <- final
prec$CJPtotal <- prec$CJPhosp
table(prec$CJPtotal)


#---------------------------------
#selection des variables d'interet
#---------------------------------
sub <- final[, c("tnt", "creatpre","age","imc", "met", "sex", "asa", "lee", "hta", "coronaroP","icc","aomi","diab", "acfa", "avc", "CJPtotal")]

sub <- tardif[, c("tnt", "creatpre","age","imc", "met", "sex", "asa", "lee", "hta", "coronaroP","icc","aomi","diab", "acfa", "avc", "CJPtotal")]
sub <- prec[, c("tnt", "creatpre","age","imc", "met", "sex", "asa", "lee", "hta", "coronaroP","icc","aomi","diab", "acfa", "avc", "CJPtotal")]


#-----------------
#Données manquantes
#-----------------
#nombre de patients avec au moins une donnée manquante : 15 patients soit <5%
prop.table(table(apply(sub,1, function(x)sum(is.na(x)))>0))
#nombre de NA par colonnes
apply(sub,2, function(x)sum(is.na(x)))

#Je supprime les NA
work <-  na.omit(sub)


#----------------
#Datamanagement
#----------------
work$sex <- as.numeric(work$sex)-1
work$CJPtotal <- as.factor(work$CJPtotal)
var <- colnames(work)[colnames(work) != "CJPtotal"]
#var <- c("tnt", "creatpre","age","imc", "met", "sex", "asa", "lee", "hta", "coronaroP","icc","aomi","diab", "acfa", "avc")

#----------------
#virer les extrêmes
#----------------

formula.vec <- paste0("CJPtotal ~ ", paste(var, collapse="+"))
mod.lg <- glm(as.formula(formula.vec), family = "binomial", data = work)
a.obj <- stepAIC(mod.lg)
plot(glm(formula = CJPtotal ~ tnt + creatpre + age + sex + hta + aomi + acfa + avc, family = "binomial", data = work))

summary(get(a.obj$call))
#a supprimer pour evenements tardifs
work <- work[! rownames(work) %in% c(74, 110), ] 
#a supprimer pour evenements précoces
work <- work[! rownames(work) %in% c(110, 168, 322), ]

#CJPtotal ~ tnt + age + hta + aomi + acfa + avc

#2e tour de validation
mod.lg <- glm(as.formula(formula.vec), family = "binomial", data = work)
a.obj <- stepAIC(mod.lg)
plot( glm(as.formula(formula.vec), family = "binomial", data = work))
#a supprimer pour evenements précoces
work <- work[! rownames(work) %in% 10, ]

#3e tour de validation
mod.lg <- glm(as.formula(formula.vec), family = "binomial", data = work)
a.obj <- stepAIC(mod.lg)
#Il n'y a plus tnt...
plot( glm(as.formula(formula.vec), family = "binomial", data = work))


#Je reprend modèle avec tnt
plot(glm(formula = CJPtotal ~ tnt + creatpre + age + sex + hta + aomi + acfa + avc, family = "binomial", data = work))
work <- work[! rownames(work) %in% c(262, 165, 179), ]
mod.lg <- glm(formula = CJPtotal ~ creatpre + age + sex + hta + aomi + acfa + avc, family = "binomial", data = work)
a.obj <- stepAIC(mod.lg)

#-----------------
#Lasso
#-----------------
x <- as.matrix(work[ , var])
y <- work$CJPtotal


fit <- glmnet(x, y, alpha = 1, family = "binomial") #1 is lasso (it's the default anyway)
plot(fit,label = TRUE)
set.seed(12345)
set.seed(123)
cvfit = cv.glmnet(x, y, family = "binomial")
cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "auc")
plot(cvfit)
#recuperer le lambda
lambda1 <- cvfit$lambda.1se
lambda1 <- cvfit$lambda.min
a <- coef(cvfit, s = lambda1)
a


# #Pour récupérer le tableau des beta obtenus par lasso
vec <- rep(NA, a@Dim[1])
vec[a@i] <- a@x[-1]
names.a <- a@Dimnames[[1]]
data.frame(variable = as.character(names.a[(a@i)+1])[-1], beta = as.numeric(vec[(a@i)+1])[-length(a@i)], stringsAsFactors = FALSE)


#paramètres utilisés
cvfit = cv.glmnet(x, y, family = "binomial")
lambda1 <- cvfit$lambda.1se
a <- coef(cvfit, s = lambda1)
a

vec <- rep(NA, a@Dim[1])
vec[a@i] <- a@x[-1]
names.a <- a@Dimnames[[1]]
var.df <- data.frame(variable = as.character(names.a[(a@i)+1])[-1], beta = as.numeric(vec[(a@i)+1])[-length(a@i)], stringsAsFactors = FALSE)
boot.var <- var.df[!is.na(var.df$beta), "variable"]



#répéter le processus en gardant le même lambda
cvfit2 <- cv.glmnet(x, y, family = "binomial")
a2 <- coef(cvfit2, s = lambda1)
a2

#------------------------------------
# intervalle de confiance : bootstrap (version Benjamin)
coefboot <- matrix(rep(NA, 20000), ncol=length(boot.var), nrow=10000)
for (i in 1:10000) {
  x <- sample(1:nrow(work), nrow(work), replace=TRUE)
  #predicteur <- matrix[x, ]
  predicteur <- as.matrix(work[x, boot.var])
  predit <- as.character(work$CJPtotal[x])
  if (all(predit == 1) | all(predit == 0)) {
    coefboot.df[i, ] <- rep(NA, length(boot.Var))
  } else {
    fit <- glmnet(predicteur, predit, family="binomial", 
                  alpha=1, lambda=0.04655163)
    coefboot.df[i, ] <- coef(fit)[-1] 
    cat(i, "\n")
  }
}
res <- c(paste(round(exp(coef(lasso)[2]),3), " [", round(exp(quantile(coefboot[, 1], probs=c(0.025))),3),";",round(exp(quantile(coefboot[, 1], probs=c(0.975))),3) , "]", sep=""),
         paste(round(exp(coef(lasso)[3]),3), " [", round(exp(quantile(coefboot[, 2], probs=c(0.025))),3),";",round(exp(quantile(coefboot[, 2], probs=c(0.975))),3), "]", sep=""),
         paste(round(exp(coef(lasso)[4]),3), " [", round(exp(quantile(coefboot[, 3], probs=c(0.025))),3),";",round(exp(quantile(coefboot[, 3], probs=c(0.975))),3), "]", sep=""),
         paste(round(exp(coef(lasso)[9]),3), " [", round(exp(quantile(coefboot[, 4], probs=c(0.025))),3),";",round(exp(quantile(coefboot[, 4], probs=c(0.975))),3), "]", sep=""))
write.csv2(res, "res_OR_IC95.csv")
#------------------------------------

data <- work

boot.lasso <- function(data, indices){
  data <- data[indices, ]
  data <- na.omit(data)
  x <- as.matrix(data[ , boot.var])
  y <- data$CJPtotal
  cvfit <- cv.glmnet(x, y, family = "binomial")
  a <- coef(cvfit, s = lambda1)
  res <- as.numeric(a)[-1]
  # browser()
  # vec <- rep(NA, a@Dim[1])
  # vec[a@i] <- a@x[-1]
  # res <- as.numeric(vec) ; names(res) <- as.character(a@Dimnames[[1]])
  
  #res <- as.list(res)
  # res <- data.frame(variable = as.character(a@Dimnames[[1]]), beta = as.numeric(vec), stringsAsFactors = FALSE)
  # res <- res[3,2]
  return(res)
}
res <- boot(data = work, statistic = boot.lasso, R = 1000)


#intervalle de confiance
n <- length (res$t0)
.list.ci <- lapply(1:n, function(x) boot.ci(res,index=x,type="perc"))
.list.ci <- lapply(1:n, function(x) boot.ci(res,index=x,type="bca"))




