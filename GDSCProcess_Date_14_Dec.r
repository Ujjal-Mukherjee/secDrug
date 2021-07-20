rm(list=ls())
setwd("C:/Users/ukm/Box/GDSC release-8.2 11-24-2020")

library(sqldf)
library(readxl)

d1 = read_xlsx("GDSC1_Fitted_Dose_Response_25Feb20.xlsx")

d2 = read_xlsx("GDSC2_Fitted_Dose_Response_25Feb20.xlsx")

x=unique(d1$CELL_LINE_NAME)
y=unique(d2$CELL_LINE_NAME)

unx = as.data.frame(intersect(x,y))
colnames(unx) = c("CELL_LINE_NAME")

dux = sqldf("SELECT * FROM d1
			WHERE d1.CELL_LINE_NAME IN (SELECT CELL_LINE_NAME FROM unx)")
duy = sqldf("SELECT * FROM d2
			WHERE d2.CELL_LINE_NAME IN (SELECT CELL_LINE_NAME FROM unx)")
duz = rbind(dux, duy)
save(duz, file="Combined_GDSC_Jan_28_21.RData")

cline = read.csv("BCellLineNames.csv")
colnames(cline)=c("CELL_ID", "CELL_LINE_NAME", "DATABASE")

duw = sqldf("SELECT * FROM duz
			WHERE duz.CELL_LINE_NAME IN (SELECT CELL_LINE_NAME FROM cline)");
save(duw, file = "B_CELL_LINE_ONLY_GDSC_JAN_28_21.RData")

dtemp = duw[,c("COSMIC_ID", "CELL_LINE_NAME", "TCGA_DESC", "DRUG_NAME", "LN_IC50")]

library(dplyr)

d3 = reshape(dtemp, timevar = "DRUG_NAME", idvar = "CELL_LINE_NAME", v.names="LN_IC50", direction = "wide")

save(d3, file = "B_CELL_WIDE_FORMAT_ALL_DRUGS_JAN_28_21.RData")

fracNA = apply(d3, 2, FUN = function(x)return(1-length(na.omit(x))/length(x)))

ind = which(fracNA >= 0.1)

d4 = d3[,-ind]

save(d4, file = "B_CELL_WIDE_FORMAT_ALL_DRUGS_REMOVED_MISSING_GE_10_PERC_JAN_28_21.RData")


########################################
#IMPUTE d4
########################################

load("B_CELL_WIDE_FORMAT_ALL_DRUGS_REMOVED_MISSING_GE_10_PERC_JAN_28_21.RData")

library(Amelia)

miss_row = apply(d4, 1, FUN=function(x)return(1 - length(na.omit(x))/length(x)))
miss_col = apply(d4, 2, FUN=function(x)return(1 - length(na.omit(x))/length(x)))

barplot(miss_row)
grid(nx=0, ny=8, col=2)

barplot(miss_col)
grid(nx=0, ny=8, col=2)

d5 = d4[,c(1:3,4:20)]
aout = amelia(d5, m = 25, idvars = c(1:3))

x = (aout$imputations$imp1[,4:20]+
	aout$imputations$imp2[,4:20]+
	aout$imputations$imp3[,4:20]+
	aout$imputations$imp4[,4:20]+
	aout$imputations$imp5[,4:20]+
	aout$imputations$imp6[,4:20]+
	aout$imputations$imp7[,4:20]+
	aout$imputations$imp8[,4:20]+
	aout$imputations$imp9[,4:20]+
	aout$imputations$imp10[,4:20]+
	aout$imputations$imp11[,4:20]+
	aout$imputations$imp12[,4:20]+
	aout$imputations$imp13[,4:20]+
	aout$imputations$imp14[,4:20]+
	aout$imputations$imp15[,4:20]+
	aout$imputations$imp16[,4:20]+
	aout$imputations$imp17[,4:20]+
	aout$imputations$imp18[,4:20]+
	aout$imputations$imp19[,4:20]+
	aout$imputations$imp20[,4:20]+
	aout$imputations$imp21[,4:20]+
	aout$imputations$imp22[,4:20]+
	aout$imputations$imp23[,4:20]+
	aout$imputations$imp24[,4:20]+
	aout$imputations$imp25[,4:20])/25

d4[,c(4:20)]=x


index = seq(5,435,5)

for(i in index){

d5 = d4[,c(1:3,i:(min(440,i+15)))]

if(dim(na.omit(d5))[1] < dim(d5)[1]){

aout = amelia(d5, m = 25, idvars = c(1:3))

x = (aout$imputations$imp1[,4:19]+
	aout$imputations$imp2[,4:19]+
	aout$imputations$imp3[,4:19]+
	aout$imputations$imp4[,4:19]+
	aout$imputations$imp5[,4:19]+
	aout$imputations$imp6[,4:19]+
	aout$imputations$imp7[,4:19]+
	aout$imputations$imp8[,4:19]+
	aout$imputations$imp9[,4:19]+
	aout$imputations$imp10[,4:19]+
	aout$imputations$imp11[,4:19]+
	aout$imputations$imp12[,4:19]+
	aout$imputations$imp13[,4:19]+
	aout$imputations$imp14[,4:19]+
	aout$imputations$imp15[,4:19]+
	aout$imputations$imp16[,4:19]+
	aout$imputations$imp17[,4:19]+
	aout$imputations$imp18[,4:19]+
	aout$imputations$imp19[,4:19]+
	aout$imputations$imp20[,4:19]+
	aout$imputations$imp21[,4:19]+
	aout$imputations$imp22[,4:19]+
	aout$imputations$imp23[,4:19]+
	aout$imputations$imp24[,4:19]+
	aout$imputations$imp25[,4:19])/25

d4[,c(i:(min(440,i+15)))]=x
}
}

for(i in 4:440){

	d4[,i]=exp(d4[,i])
	
}


save(d4, file = "B_CELL_WIDE_FORMAT_ALL_IMPUTED_FEB_05_21.RData")


########################################
#GROUP AS S, N, and R
######################################## 


AssignGroup=function(x,l){

	ina=which(is.na(x));
	x[ina]=mean(x[-ina]);

	n=length(x);
	n1=ceiling(l*n);
	n2=floor((1-l)*n);

	class=rep("N",n);
	ind=order(x);
	class[ind[1:n1]]=rep("S",n1);
	class[ind[n2:n]]=rep("R",((n-n2)+1));

	return(class);

}





load("Data_Processed_GDSC1_Selected_Columns_Below_60_Per_Missing.RData")

missing = apply(d6, 2, FUN=function(x){return(1-length(na.omit(x))/length(x))})

library(sqldf)

grby = function(x, r){
	ur = unique(r)
	missing = c()
	for(i in 1:length(ur)){
		x1 = x[which(r==ur[i])]
		miss = length(na.omit(x1))/length(x1)
		missing = c(missing, miss)
	}
	return(missing)
}

missbygr = apply(d6[,-c(1:4)], 2, FUN = grby, r = d6$TCGA_DESC)
missbygr = cbind(Tissue = unique(d6$TCGA_DESC), missbygr)

View(missbygr)

d6 = d6[,order(missing)]

rowmiss = apply(d6[,c(5:445)], 1, FUN=function(x){return(1-length(na.omit(x))/length(x))})

ind = which(rowmiss > 0.6)

d6 = d6[-ind,]

library(Amelia)

for(i in seq(95,445,1)){
	d7 = amelia(d6[,c(1:4,5:i)], m = 25, idvars = c(1:4))
	Impute = (d7$imputations$imp1[,c(5:i)] + 
				d7$imputations$imp2[,c(5:i)] +
				d7$imputations$imp3[,c(5:i)] + 
				d7$imputations$imp4[,c(5:i)] +
				d7$imputations$imp5[,c(5:i)] + 
				d7$imputations$imp6[,c(5:i)] +
				d7$imputations$imp7[,c(5:i)] + 
				d7$imputations$imp8[,c(5:i)] +
				d7$imputations$imp9[,c(5:i)] + 
				d7$imputations$imp10[,c(5:i)] +
				d7$imputations$imp11[,c(5:i)] + 
				d7$imputations$imp12[,c(5:i)] +
				d7$imputations$imp13[,c(5:i)] + 
				d7$imputations$imp14[,c(5:i)] +
				d7$imputations$imp15[,c(5:i)] + 
				d7$imputations$imp16[,c(5:i)] +
				d7$imputations$imp17[,c(5:i)] + 
				d7$imputations$imp18[,c(5:i)] +
				d7$imputations$imp19[,c(5:i)] + 
				d7$imputations$imp20[,c(5:i)] +
				d7$imputations$imp21[,c(5:i)] + 
				d7$imputations$imp22[,c(5:i)] +
				d7$imputations$imp23[,c(5:i)] + 
				d7$imputations$imp24[,c(5:i)] +
				d7$imputations$imp25[,c(5:i)])/25
		d6[,c(5:i)] = Impute
		save(Impute, file=paste("Interim_Impute_Columns_5_",i,"_.RData",sep=""))
}














