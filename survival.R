##### Overall survival ######### 
library(survival)
library(rms)


data_uchl1$OS_M <- as.numeric(as.character(data_uchl1$CURATED_DAYS_TO_DEATH_OR_LAST_FU))/30.4
fit = npsurv(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
               UCHL1_G, data = data_uchl1)
fit


strata = levels(data_uchl1$UCHL1_G)

library(Cairo)
CairoSVG(file = "OS.svg",  width = 6, height = 6, 
         onefile = TRUE, bg = "transparent",
         pointsize = 12)
par(mar=c(6,4,2,8), mgp = c(2, 1, 0))
survplot(fit,
         time.inc = 12,
         xlab = 'Months',
         lty = c(1:2),
         conf="none", add=FALSE, 
         label.curves=FALSE, abbrev.label=FALSE,
         levels.only=TRUE, lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE,n.risk=TRUE,logt=FALSE,
         dots=FALSE,
         grid=FALSE,
         srt.n.risk=0, sep.n.risk=0.04, adj.n.risk=0.5, 
         y.n.risk=-0.25, cex.n.risk=0.6, pr=FALSE       
)

legend(60, 1.0, strata, lty = c(1:2), cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

legend(5.5*12, 0.85, 'P-value = 0.445', cex = 0.8,
       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1,
       trace = TRUE,
       bty = 'n')

dev.off()

diff = survdiff(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
                  UCHL1_G, data = data_uchl1)
diff


data_uchl1$stage <- 
  factor(data_uchl1$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE %in% 
           c('Stage III' ,'Stage IV'),
         levels = c(FALSE, TRUE),
         labels = c('Stage 0-II', 'Stage III-IV'))
data_uchl1$stage[is.na(data_uchl1$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE)] <- NA
summary (data_uchl1$CURATED_PATHOLOGIC_STAGE_AJCC7_AT_DIAGNOSIS_SIMPLE)

cox <- coxph(Surv(OS_M, CURATED_VITAL_STATUS == 'Dead')~
               stage + 
               CURATED_AGE_AT_TCGA_SPECIMEN+
               UCHL1_G, data = data_uchl1)

summary(cox)

sink('cox_analysis_output.txt')