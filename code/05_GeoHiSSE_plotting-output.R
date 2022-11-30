# This is an R script to plot output from log files of GeoHiSSE/MuSSE/MuHiSSE analyses
# presented in McCullough et al. 2022 systematic biology. 

# written by Jenna McCullough and Rosana Zenil-Ferguson 
# last edited on 15 June 2022 by Jenna McCullough 

# Modified by Putter Tiatragul
# Last modified November 2022

library(devtools)
install_github("revbayes/RevGadgets")
library(RevGadgets)
library(ggplot2)

#multiplot.R needs to be in your working directory for this to work. 
source("code/utility/multiplot.R")

#Remember, colours need to be in alphabetical order of your parameters!

###########################################
###########################################
# Model 1 ANALYSIS: GEOHISSE non_arid
###########################################
###########################################
non_arid_1.sse<-read.table("data/intermediate_data/geohisse/00_GeoHiSSE-blindsnake/output/comb_geohisse_blindsnake.log", header=TRUE)
head(non_arid_1.sse)
str(non_arid_1.sse)

#anagenetic dispersal for model 1 
anagenetic.disp_1<-data.frame(dens=c(non_arid_1.sse$anagenetic_dispersal_21, non_arid_1.sse$anagenetic_dispersal_31, non_arid_1.sse$anagenetic_dispersal_54, 
                                     non_arid_1.sse$anagenetic_dispersal_64),Type=rep(c("arid to widespread A","Island to widespread A","arid to widespread B","Island to widespread B"), 
                                                                                     each=length(non_arid_1.sse$anagenetic_dispersal_21)))
anagenetic.colors=c("#336600","#330033","#66CC66","#996699") 

#arids A = #336600 
#non_arid A = #66CC66 
#arids B = #330033
#non_arid B = #996699 

p1_1<-ggplot(anagenetic.disp_1, aes(x=dens, fill=Type))+labs(title="Anagenetic dispersal for non_arid v. arids",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =anagenetic.colors ) +xlim(0,2)
p1_1

#Extinction for model 1 
extinction_1<-data.frame(dens=c(non_arid_1.sse$extinction_rates.2.,non_arid_1.sse$extinction_rates.3., non_arid_1.sse$extinction_rates.5., non_arid_1.sse$extinction_rates.6.),Type=rep(c("arids A", "non_arid A", "arids B", "non_arid B"), each=length(non_arid_1.sse$extinction_rates.5.)))
extinction.colors=c("#336600","#330033","#66CC66","#996699") 
p2_1<-ggplot(extinction_1, aes(x=dens, fill=Type))+labs(title="Extinction for non_arid v. contients",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = extinction.colors )+xlim(0,1.3)
p2_1

#Diversification for model 1 

# non_arid A : speciation_rates.2. + speciation_rates.5. + speciation_rates.6. + speciation_rates.7. + speciation_rates.8.  - extinction_rates.3.
# non_arid B : speciation_rates.10. + speciation_rates.13. + speciation_rates.14. + speciation_rates.15. + speciation_rates.16. - extinction_rates.6.
# arids A : speciation_rates.1. + speciation_rates.3. + speciation_rates.4. + speciation_rates.7. + speciation_rates.8. - extinction_rates.2.
# arids B : speciation_rates.9. + speciation_rates.11. + speciation_rates.12. + speciation_rates.15. + speciation_rates.16. - extinction_rates.5.
diversification_1<-data.frame(dens=c(non_arid_1.sse$speciation_rates.2. + non_arid_1.sse$speciation_rates.5. + non_arid_1.sse$speciation_rates.6. + non_arid_1.sse$speciation_rates.7. + non_arid_1.sse$speciation_rates.8.  - non_arid_1.sse$extinction_rates.3.,
                                              non_arid_1.sse$speciation_rates.1. + non_arid_1.sse$speciation_rates.3. + non_arid_1.sse$speciation_rates.4. + non_arid_1.sse$speciation_rates.7. + non_arid_1.sse$speciation_rates.8. - non_arid_1.sse$extinction_rates.2.,
                                              non_arid_1.sse$speciation_rates.9. + non_arid_1.sse$speciation_rates.11. + non_arid_1.sse$speciation_rates.12. + non_arid_1.sse$speciation_rates.15. + non_arid_1.sse$speciation_rates.16. - non_arid_1.sse$extinction_rates.5.,
                                              non_arid_1.sse$speciation_rates.10. + non_arid_1.sse$speciation_rates.13. + non_arid_1.sse$speciation_rates.14. + non_arid_1.sse$speciation_rates.15. + non_arid_1.sse$speciation_rates.16.  - non_arid_1.sse$extinction_rates.6.),
                                       Type=rep(c("non_arid A", "arids A", "arids B", "non_arid B" ), each=length(non_arid_1.sse$speciation_rates.1.)))

diversification.colors=c("#336600","#330033","#66CC66","#996699")
p3_1<-ggplot(diversification_1, aes(x=dens, fill=Type))+labs(title="Diversification of non_arid v. contients",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = diversification.colors )+xlim(0,4)
p3_1

# uncomment when it's time to make PDFs 
#pdf("final-output/1_ana_ext_div-v3.pdf", width = 5, height = 10)
#multiplot(p1_1, p2_1, p3_1)
#dev.off()
 
 
#hidden rates of model 1
hidden_rates_1<-data.frame(dens=c(non_arid_1.sse$hidden_rate_1, non_arid_1.sse$hidden_rate_2),Type=rep(c("Hidden rate 1","Hidden rate 2"), each=length(non_arid_1.sse$anagenetic_dispersal_21)))
hidden.colors=c("black", "gray") 
p4_1<-ggplot(hidden_rates_1, aes(x=dens, fill=Type))+labs(title="hidden rates of model 1",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = hidden.colors )+xlim(0,2)
p4_1

#pdf("final-output/1_hidden-rates.pdf", width = 5, height = 3.3)
#dev.off()

###########################################
###########################################
# model 2 ANALYSIS: GEOHISSE non_arid EXTENDED MODEL 
###########################################
###########################################


non_arid_2_extended.sse<-read.table("02_non_arid-extended/geohisse-non_arid-680K-25burnin-resamp45.log", header=TRUE)
head(non_arid_2_extended.sse)# just to check what I read
str(non_arid_2_extended.sse)

#anagenetic dispersal for model 2 
anagenetic.disp_2<-data.frame(dens=c(non_arid_2_extended.sse$anagenetic_dispersal_21, non_arid_2_extended.sse$anagenetic_dispersal_31, non_arid_2_extended.sse$anagenetic_dispersal_54, non_arid_2_extended.sse$anagenetic_dispersal_64),Type=rep(c("arid to widespread A","Island to widespread A","arid to widespread B","Island to widespread B"), each=length(non_arid_2_extended.sse$anagenetic_dispersal_21)))
anagenetic.colors=c("#336600","#330033","#66CC66","#996699") 

p1_2<-ggplot(anagenetic.disp_2, aes(x=dens, fill=Type))+labs(title="Anagenetic dispersal for non_arid v. arids extended",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values =anagenetic.colors )+xlim(0,1.3)+ylim(0,75)
p1_2

#Extinction for model 2 
extinction_2<-data.frame(dens=c(non_arid_2_extended.sse$extinction_rates.2.,non_arid_2_extended.sse$extinction_rates.3., non_arid_2_extended.sse$extinction_rates.5., non_arid_2_extended.sse$extinction_rates.6.),Type=rep(c("arids A", "non_arid A", "arids B", "non_arid B"), each=length(non_arid_2_extended.sse$extinction_rates.5.)))
extinction.colors=c("#336600","#330033","#66CC66","#996699") 
p2_2<-ggplot(extinction_2, aes(x=dens, fill=Type))+labs(title="Extinction for non_arid v. arids Extended",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = extinction.colors )+xlim(0,1.3)
p2_2


# Diversification for  model 2 
# non_arid A : speciation_2 + speciation_5 + speciation_7 - extinction_rates.3.
# non_arid B : speciation_10 +speciation_13 + speciation_15 - extinction_rates.6.
# arids A : speciation_1 + speciation_3 + speciation_7 - extinction_rates.2.
# arids B : speciation_9 + speciation_11 + speciation_15 - extinction_rates.5.

diversification_2 <- data.frame(dens=c(non_arid_2_extended.sse$speciation_2 + non_arid_2_extended.sse$speciation_5 + non_arid_2_extended.sse$speciation_7 - non_arid_2_extended.sse$extinction_rates.3.,
                                       non_arid_2_extended.sse$speciation_10 + non_arid_2_extended.sse$speciation_13 + non_arid_2_extended.sse$speciation_15 - non_arid_2_extended.sse$extinction_rates.6.,
                                       non_arid_2_extended.sse$speciation_1 + non_arid_2_extended.sse$speciation_3 + non_arid_2_extended.sse$speciation_7 - non_arid_2_extended.sse$extinction_rates.2.,
                                       non_arid_2_extended.sse$speciation_9 + non_arid_2_extended.sse$speciation_11 + non_arid_2_extended.sse$speciation_15 - non_arid_2_extended.sse$extinction_rates.5.),
                                Type=rep(c("non_arid A", "non_arid B", "arids A", "arids B"), each=length(non_arid_2_extended.sse$extinction_rates.3.)))

diversification.colors=c("#336600","#330033","#66CC66","#996699")
p3_2<-ggplot(diversification_2, aes(x=dens, fill=Type))+labs(title="diversification of non_arid v. arids extended model",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = diversification.colors )+xlim(0,5)
p3_2

# comparing types of speciation in non_arid v. arids  for model 2 in extended speciations 
# arids A					    	# non_arid A
# widespread sympatry = 1 		# widespread sympatry = 2
# subset sympatry = 3				  # subset sympatry = 5
# allopatric = 7				    	# allopatric = 7 

# arids B					  	# non_arid B
# widespread sympatry = 9 	# widespread sympatry = 10
# subset sympatry = 11			# subset sympatry = 13
# allopatric = 15					  # allopatric = 15 

speciation_2 <- data.frame(dens=c(non_arid_2_extended.sse$speciation_1, non_arid_2_extended.sse$speciation_3, non_arid_2_extended.sse$speciation_7, 
                                  non_arid_2_extended.sse$speciation_2 , non_arid_2_extended.sse$speciation_5, non_arid_2_extended.sse$speciation_7, 
                                  non_arid_2_extended.sse$speciation_9, non_arid_2_extended.sse$speciation_11, non_arid_2_extended.sse$speciation_15,
                                  non_arid_2_extended.sse$speciation_10, non_arid_2_extended.sse$speciation_13, non_arid_2_extended.sse$speciation_15), 
                           Type=rep(c("arids A WS", "arids A SS", "arids A Allo", "non_arid A WS", "non_arid A SS", "non_arid A Allo", "arids B WS", "arids B SS", "arids B Allo",  "non_arid B WS", "non_arid B SS", "non_arid B Allo"), 
                                    each=length(non_arid_2_extended.sse$speciation_15)))    

p4_2<-ggplot(speciation_2, aes(x=dens, fill=Type))+labs(title="Types of speciation in non_arid V. arids Extended",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlim(-0.1,2.25)
p4_2

# make a violin plot 
dp1 <- ggplot(speciation_2, aes(x=Type, y=dens,fill=Type)) + geom_violin(trim=FALSE)+
  labs(title="Types of speciation in non_arid V. arids Ext",x="State", y = "Posterior Density")
dp1<-dp1 + theme_classic() +ylim(0,2)
dp1

#pdf("final-output/2_speciation-types.pdf", width = 4, height = 3)
#dev.off()

#pdf("final-output/2_ana_ext_div-v2.pdf", width = 5, height = 10)
#multiplot(p1_2, p2_2, p3_2)
#dev.off()
# 

# hidden rates of model 2 
hidden_rates_2a<-data.frame(dens=c(non_arid_2_extended.sse$hidden_rate_1, non_arid_2_extended.sse$hidden_rate_2),Type=rep(c("Hidden rate 1","Hidden rate 2"), each=length(non_arid_2_extended.sse$anagenetic_dispersal_21)))

p4_2_hidden<-ggplot(hidden_rates_2a, aes(x=dens, fill=Type))+labs(title="hidden rates of model 2",x="Rate", y="Posterior Density")+geom_density(alpha=0.7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_fill_manual( values = hidden.colors )+scale_x_continuous(limits=c(0,2.5))
p4_2_hidden


#pdf("final-output/2_hidden-rates.pdf", width = 5, height = 3.3)
#dev.off()


