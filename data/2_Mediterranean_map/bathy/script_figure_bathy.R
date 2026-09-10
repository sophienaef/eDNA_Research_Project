### Loading of library
library(rgdal)
library(gpclib)
library(maptools)
gpclibPermit()
library(plotrix)

setwd("/Volumes/Data_perso/These/Gap_analyses/New-gap")


### Loading of Library gpc function: getGPC;rotate.poly;centroids.poly ;avCentro.poly ;translate.poly;inverse.poly
source("gpcLibrary_20052011.R")

### Loading of Mediterranean countries shape file.
essai=readShape2Gpc('Med_cont_pays.shp')

### Loading of Mediterranean shelf
Gpc_plato <- readShape2Gpc('med_et2_groupe.shp')

### Loadinf of ecoregion
#ecoreg <-  readShape2Gpc('/Volumes/Data_model/MED/Matos/Ecoregion/Eco_region_med.shp')

ecoregion <-  readShape2Gpc('/Volumes/Data_perso/These/Beta_temporel/essai_inter_med_ecoreg.shp')

setwd("/Volumes/Data_perso/These/Redaction_these/Figures/Figure_bathy_correction/")


Coris_poly  <- readShape2Gpc('Coris_julis.shp')
Bathy  <- readShapeSpatial('med_etopo2_proj.shp')

coords= read.csv2("coord_plato.csv",sep="\t",dec=".")

data=data.frame(coordX=coords[,1],coordY = coords[,2])
coord_m <- get_proj(data=data[,1:2],ref1 ="+proj=longlat +datum=WGS84",  ref2="+init=epsg:3035")

#Recupération des centroides
res=lapply(Bathy@polygons,function(x){
data= x@Polygons[[1]]@coords
data =data[-dim(data)[1],]
res = apply(data,2,mean)
return(res)
})

essai1=do.call(rbind,res)

#Recupération des données de bathy

prof=Bathy[[1]]

essai2 = data.frame(coordx=essai1[,1],coordy= essai1[,2],prof=prof)
essai2 = essai2[essai2$prof<0,]


tiff(filename = "Figure_Bathy.tiff",width = 7.5, height = 4, units = "in", res=600,compression = "lzw")

################# Plot coris julis de base #################
layout(mat= rbind(c(1,2),c(3,3)))

par(mar=c(0.75,1.5,1.5,1))

lala(coordX=coord_m[,1] ,coordY=coord_m[,2],vect=0,cex_pt = 0.21,lwd_bb = 0.7,tcl=-0.25,cex_axis=0.1,cex_legend=0.65,names_fig= "a)",cex_text_fig =1.2,ecoreg=FALSE,plato=TRUE,col.plat="lightskyblue",poly.sup=FALSE,text.legend="Aire de répartition")


############## Plot_bathy##################
 par(mar=c(0.75,1.5,1.5,1))

## plot bathy
lala_bathy(x=essai2[,3])

############## Plot coris julis modifiér

par(mar=c(1,12,1.5,12))
lala(coordX=coord_m[,1] ,coordY=coord_m[,2],vect=0,cex_pt = 0.21,lwd_bb = 0.7,tcl=-0.25,cex_axis=0.1,cex_legend=0.65,names_fig= "c)",cex_text_fig =1.2,ecoreg=FALSE,plato=TRUE,col.plat="lightskyblue",col.poly="royalblue3",poly.sup=TRUE, names.poly.sup=Coris_poly,text.legend=c("Aire de répartition corrigée","Ancienne aire de réaprtition"),col.legend=c("royalblue3","lightskyblue"))


dev.off()

###############################################################
#Fonction graphique

lala=function(coordX=coord_m[,1] ,coordY=coord_m[,2],vect=seq(160,260,10),x= dat$Richesse,include.lowest=TRUE,colour=c("#FFFF80FF","#FFFF00FF","#FFCC00FF","#FF9900FF","#FF6600FF","#FF3300FF","#FF0000FF"),cex_pt = 0.21,lwd_bb = 0.7,tcl=-0.25,cex_axis=0.1,cex_legend=0.5,names_fig= "a)",cex_text_fig =1.2,ecoreg=TRUE,plato=FALSE,col.plat="lightskyblue",col.legend="lightskyblue",poly.sup=TRUE, names.poly.sup=Coris_poly,col.poly="lightskyblue",text.legend=""){
	
   
   if (sum(vect!=0)){
    vect= vect
	richesse = cut(x, breaks=vect,include.lowest=include.lowest)

	colour= colour
    
	plot(rnorm(100),col="white",ylim = c(840000,2524929),xlim = c( 2874611,6683749),type="n",axes=F,xlab="",ylab="")
	
    if(ecoreg==TRUE){
	lapply(ecoregion,function(x){plot(x,poly.args=list(col="white",lwd=0.2),add=TRUE,axis=F)})
    }
    
	if(plato==TRUE){
   	plot(Gpc_plato,poly.args=list(col=col.plat,lwd=0.1),add=T)
    }
    
    lapply ( essai,function (x) { plot ( x,poly.args = list (col = "#EFEFEF",lwd = 0.3),add = TRUE,axis = F) })	
    
        
    points(coordX,coordY,col = colour[richesse], pch = 15,cex = 0.19) 
    
    
    mtext(c("0°","10°E","20°E","30°E"),side=1,line=0,at=c(3e+06,4e+06,5e+06,6e+06),cex=0.5,bg="white")
    mtext(c("30°N","35°N","40°N","45°N"),side=2,line=0.4,at=c(1e+06,1.5e+06,2e+06,2.5e+06),cex=0.5,bg="white")
    
    axis(side=1,line=0,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)
    axis(side=2,line=0,las=2,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)
    
    mtext(names_fig,side=2,line=0.3,at= 2.75e+06,cex=0.6,bg="white",las=2)    
    
    box(lwd=lwd_bb)

	
	legend("topright",leg= levels(richesse),pch=15,col=colour,pt.cex=1.2,bg="white",cex=cex_legend,lwd=0.2,box.lwd=lwd_bb,lty=0 )
}
else {
	
	    
	plot(rnorm(100),col="white",ylim = c(840000,2524929),xlim = c( 2874611,6683749),type="n",axes=F,xlab="",ylab="")
	
    if(ecoreg==TRUE){
	lapply(ecoregion,function(x){plot(x,poly.args=list(col="white",lwd=0.2),add=TRUE,axis=F)})
    }
    
	if(plato==TRUE){
   	plot(Gpc_plato,poly.args=list(col=col.plat,lwd=0.1),add=T)
    }
    
    if(poly.sup==TRUE){
   	plot(names.poly.sup,poly.args=list(col=col.poly,lwd=0.1),add=T)
    }

    
    lapply ( essai,function (x) { plot ( x,poly.args = list (col = "#EFEFEF",lwd = 0.3),add = TRUE,axis = F) })	
        
    
    mtext(c("0°","10°E","20°E","30°E"),side=1,line=0,at=c(3e+06,4e+06,5e+06,6e+06),cex=0.5,bg="white")
    mtext(c("30°N","35°N","40°N","45°N"),side=2,line=0.4,at=c(1e+06,1.5e+06,2e+06,2.5e+06),cex=0.5,bg="white")
    
    axis(side=1,line=0,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)
    axis(side=2,line=0,las=2,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)
    
    mtext(names_fig,side=2,line=0.3,at= 2.75e+06,cex=0.6,bg="white",las=2)    
    
    legend("topright",leg= text.legend,pch=15,col=col.legend,pt.cex=1.2,bg="white",cex=cex_legend,lwd=0.2,box.lwd=lwd_bb,lty=0 )

    
    box(lwd=lwd_bb)	
	
	}

	
	
} # end of function lala 


get_proj = function(data,ref1 ="+proj=longlat +datum=WGS84",  ref2="+init=epsg:3035"){ 
    
	coordinates(data)<- ~coordX+coordY
	proj4string(data) <- CRS(ref1) 
	DataUTM <- spTransform(data, CRS(ref2))
    
    return(data.frame(coordX=coordinates(DataUTM)[,1],coordY=coordinates(DataUTM)[,2]))
    
} 
############################################


res=lapply(Bathy@polygons,function(x){
data= x@Polygons[[1]]@coords
data =data[-dim(data)[1],]
res = apply(data,2,mean)
return(res)
})

essai1=do.call(rbind,res)


prof=Bathy[[1]]

essai2 = data.frame(coordx=essai1[,1],coordy= essai1[,2],prof=prof)

essai2 = essai2[essai2$prof<0,]


lala_bathy =function (x=essai2[,3]){

breaks = cut(x,breaks=c(-5000,-4000,-3000,-2000,-1000,-500,-250,-100,-50,0),dig.lab = 4)

colour= c("gray28", "#757575", "#848484", "#939393", "#A3A3A3","#B2B2B2", "#C1C1C1", "#D1D1D1", "#E0E0E0", "#EFEFEF", "#FFFFFF")

plot(rnorm(100),col="white",ylim = c(840000,2524929),xlim = c( 2874611,6683749),type="n",axes=F,xlab="",ylab="")
	points(essai2[,1],essai2[,2],col=colour[breaks],pch=21,cex=0.19)
   lapply ( essai,function (x) { plot ( x,poly.args = list (col = "#EFEFEF",lwd = 0.3),add = TRUE,axis = F) })
   legend("topright",leg= levels(breaks),pch=15,col=colour,pt.cex=1.2,bg="white",cex=0.4,lwd=0.2,box.lwd=0.7,lty=0,ncol=2 )

 mtext(c("0°","10°E","20°E","30°E"),side=1,line=0,at=c(3e+06,4e+06,5e+06,6e+06),cex=0.5,bg="white")
    mtext(c("30°N","35°N","40°N","45°N"),side=2,line=0.4,at=c(1e+06,1.5e+06,2e+06,2.5e+06),cex=0.5,bg="white")
    
    axis(side=1,line=0,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)
    axis(side=2,line=0,las=2,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)
    
    mtext("b)",side=2,line=0.3,at= 2.75e+06,cex=0.6,bg="white",las=2)    
    
  box()
	
}

################## figure occurence bathy ###############


data=read.csv("/Volumes/Data_perso/These/Aire_repartition/Data_esp_biomod/Fichier_demersal_projection_2080_2099/Dat_2080_2099_Pred_ESp_374")

tiff(filename = "Figure_Bathy_2.tiff",width = 3.5, height = 2, units = "in", res=600,compression = "lzw")


data1=data.frame(coordX=data$long,coordY=data$lat)

coord = get_proj(data1)
par(mar=c(1,1.5,1.25,1))


plot(rnorm(100),col="white",ylim = c(840000,2524929),xlim = c( 2874611,6683749),type="n",axes=F,xlab="",ylab="")
	points(coord[,1],coord[,2],col="yellowgreen",pch=15,cex=0.19)
   lapply ( essai,function (x) { plot ( x,poly.args = list (col = "#EFEFEF",lwd = 0.3),add = TRUE,axis = F) })

 mtext(c("0°","10°E","20°E","30°E"),side=1,line=0,at=c(3e+06,4e+06,5e+06,6e+06),cex=0.35,bg="white")
    mtext(c("30°N","35°N","40°N","45°N"),side=2,line=0.4,at=c(1e+06,1.5e+06,2e+06,2.5e+06),cex=0.35,bg="white")
    
    axis(side=1,line=0,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)
    axis(side=2,line=0,las=2,cex.axis=0.1,lwd=0.35,tcl=-0.25,bg="white",labels=F)

box()
dev.off()