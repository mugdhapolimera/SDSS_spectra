require("astro", quietly=FALSE)
require("magicaxis", quietly=FALSE)
require("stats",quietly=FALSE)

rdat<-read.csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpt1_filter.csv")
# a simple Te method calib. is the NII/Halpha PP04, for ???2.5 < N2 < ???0.3
# Pettini & Pagel 2004

N2<-log10(rdat$nii_6584_flux_ext/rdat$h_alpha_flux_ext)

#selgd<-which(N2>(-2.5) & N2<(-0.3))
selbd<-which(N2<(-2.5) | N2>(-0.3))
PP04_12logOH<-9.37+2.03*N2+1.26*N2^2+0.32*N2^3
rep<-which(is.na(PP04_12logOH))
PP04_12logOH[rep]<-(-99.)
PP04_12logOH[selbd]<-(-99.)

# alternative PP04 calib using O3, hbeta, and N, Halpha
O3N2=log10((rdat$oiii_5007_flux_ext/rdat$h_beta_flux_ext)/(rdat$nii_6584_flux_ext/rdat$h_alpha_flux_ext))

PP04_12logOH_v2<-8.73-0.32*O3N2

rep<-which(is.na(PP04_12logOH_v2))
PP04_12logOH_v2[rep]<-(-99.)
selbd<-which(O3N2 > 1.9) # supposedly this calibration is bad in this regime, don't use those values
PP04_12logOH_v2[selbd]<-(-99.)

aplot(PP04_12logOH, PP04_12logOH_v2, pch=16, xlab="[NII]/Halpha Z", ylab="[OIII]/Hbeta &[NII]/Halpha Z", xlim=c(7.5,10.5), ylim=c(7.5,10.5))

lines(c(7.5,11),c(7.5,11), col="red")




out<-cbind(as.vector(rdat[,"NAME"]), PP04_12logOH, PP04_12logOH_v2)

#write.csv(out, "RES_PP04Zest.csv", row.names=FALSE, quote=FALSE)



iziout<-read.table("RESOLVE_izioutv1.txt", strip.white=TRUE, skip=1)

matchi<-match(iziout$V1, rdat$NAME)
#PP04_12logOH

png(file="ResZcomp_PP04.png", width=6, height=6, units="in", res=200)

gdpts=which(!is.na(matchi))

aplot(iziout[gdpts,"V2"], PP04_12logOH[matchi[gdpts]], pch=16, xlim=c(7.5,9),ylim=c(7.5,9), xlab="12 + log(O/H)  [IZI]", ylab="12 + log(O/H)  [PP04]", col=acol("black", alpha=0.25))

lines(c(7.5,9), c(7.5,9), col="red", lwd=2)

cat("\n\n")
graphics.off()