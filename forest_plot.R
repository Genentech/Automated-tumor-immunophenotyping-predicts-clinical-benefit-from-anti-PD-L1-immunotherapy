library(gClinBiomarker)
library(survival)

imp130=read.csv("IMP130_combined.csv")
oak=read.csv("OAK_combined.csv")

outn=c("PFS","OS")
varn=c("MOCHA_BITE","Manual_IP")

#dat=oak
dat=imp130

dat$MOCHA_BITE1=as.character(dat$MOCHA_BITE)
dat$MOCHA_BITE1[dat$MOCHA_BITE!="Inflamed"]="Non-Inflamed"

dat$Manual_IP1=as.character(dat$Manual_IP)
dat$Manual_IP1[dat$Manual_IP!="Inflamed"]="Non-Inflamed"

for(i in 1:2){
for(j in 1:2){
#pdf(paste("Oak_",outn[i],"_",varn[j],".pdf",sep=""),width=12)
pdf(paste("Impassion_",outn[i],"_",varn[j],".pdf",sep=""),width=12)

pltf_p <- PlotTabForestBiomarker(data=dat,show.itt = FALSE, show.bep = FALSE,
                                    outcome.class="survival",
                                    outcome.var=c(paste(outn[i],"_MONTHS",sep=""),paste(outn[i],"_EVENT",sep="")),
                                    trt="ARM",ties="exact",#across.and.within=TRUE,
                                    var=varn[j],var.class="categorical", less=TRUE)

pltf_p <- PlotTabForestBiomarker(data=dat,show.itt = FALSE, show.bep = FALSE,
                                    outcome.class="survival",
                                    outcome.var=c(paste(outn[i],"_MONTHS",sep=""),paste(outn[i],"_EVENT",sep="")),
                                    trt="ARM",ties="exact",#across.and.within=TRUE,
                                    var=paste(varn[j],"1",sep=""),var.class="categorical", less=TRUE)

PlotKM(data=dat, tte=paste(outn[i],"_MONTHS",sep=""),cen=paste(outn[i],"_EVENT",sep=""), 
         main=paste(outn[i],"by",varn[j]), 
         trt="ARM", var=paste(varn[j],"1",sep=""),legend.loc = "topright")
#fit=coxph(Surv(oak$OS_MONTH,oak$OS_EVENT)~oak$MOCHA1*oak$ARM,ties="exact")
#text(0.9,0.8,paste("Interaction p-value =",round(summary(fit)$coefficients[3,5],3)))
dev.off()
}
}

