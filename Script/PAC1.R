## ----echo=FALSE,message=FALSE,error=FALSE,warning=FALSE-------------------------------------------------------------
library(stringr)
library(dplyr)
library(tidyverse)
library(SummarizedExperiment)
# LLlegim l'arxiu de dades que volem analitzar
data <- readLines("C:/Users/41583789L/Desktop/UOC_2nQ/Anàlisis dades òmiques/PAC1/Dades/ST000688_AN001061.txt")

## Les dades de metabolòmica comencen a la fila 124 i acaben a la 288
# Creem matriu de dades de metabolòmica (filas 124-288)
data_matrix <- read.table(text = data[124:288], sep = "\t", header = TRUE, row.names = 1)

## Copiem la informació de la primera fila i la eliminem de la martiu perquè no són dades, sino informació de quin grup pertany cada subjecte
primera_fila <-data_matrix[1,]
data_matrix <- data_matrix[-1,]

## Canviem el nom de les columnes, passem a tenir només el nom dels individus, a tenir a quin grup pertanyen i assignem un codi
colnames(data_matrix) <-  paste(c(rep("Benign MS",13),rep("SPMS",20)), c("01", "02", "03","04","05","06","07","08","09","10","11","12","13","01", "02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20"), sep="")

## Modifiquem una mica els noms de les files perquè no quedin noms tan llargs
noms <- str_extract(rownames(data_matrix), "^[^;]+")
rownames(data_matrix) <- noms
data_matrix <- as.data.frame(apply(data_matrix,2,as.numeric))
rownames(data_matrix) <- noms


# A partir de la fila 292 trobem les metadates dels metabòlits

metabolite_metadata <- read.table(text = data[292:455], sep = "\t", header = TRUE, row.names = 1,quote = "",fill=T)
metabolite_metadata <- cbind(rownames(metabolite_metadata),metabolite_metadata)
rownames(metabolite_metadata) <- NULL


# De la fila 67-100 troben la informació de les covariables dels individus
#Extraer información de individuos y factores 
sample_metadata <- read.table(text = data[67:100], sep = "\t", header = F, row.names = 2)
## Eliminem la primera columna perquè no ens interesa
sample_metadata <- sample_metadata[,-1]

## Separem la columna de totes els factors en funció de ;
sample_metadata <- sample_metadata %>%
  separate(col = 3, into = c("Age", "BMI", "Gender", "Race"), sep = "; ", extra = "merge") %>%
  mutate(Age = gsub("Age=", "", Age),
         BMI = gsub("BMI kg/m2=", "", BMI),
         Gender = gsub("Gender=", "", Gender),
         Race = gsub("Race=", "", Race))


## Canviem UNKNOWN per NA
sample_metadata <- sample_metadata %>% mutate(Age = if_else(Age=="UNKNOWN",NA,Age),
                                              BMI = if_else(BMI=="UNKNOWN",NA,BMI),
                                              Gender = if_else(Gender=="UNKNOWN",NA,Gender),
                                              Race = if_else(Race=="UNKNOWN",NA,Race))

colnames(sample_metadata) <- c("Subject","Tipus MS","Age","BMI","Gender","Race")

rownames(sample_metadata) <- NULL
sample_metadata$`Tipus MS` <- as.factor(str_extract(sample_metadata$`Tipus MS`, "(?<=:).*"))


# Desde la fila 1 a la fila 122, hi pot haver informació de l'experiment, així que la seleccione,
general_metadata <-  read.table(text = data[c(5:15,52:66,106:122)], sep = "\t", header = F)
general_metadata$V2 <- paste0(general_metadata$V1,": ",general_metadata$V2)
colnames(general_metadata)[2] <- "Experiment"

## Podem veure alguna de les matrius que hem obtingut

head(data_matrix)[1:4,c(1,3,5,14,25,33)]

head(metabolite_metadata)[1:4,]

#Guardem les metadades.

#metadades <- bind_rows(data.frame(metadata(se)), data.frame(colData(se)), data.frame(rowData(se)))
#metadades[39:nrow(metadades),1] <- rownames(metadades)[39:nrow(metadades)]
#write_csv(metadades, "metadades_data.csv")





## ----echo=FALSE,message=FALSE,error=FALSE,warning=FALSE-------------------------------------------------------------

## Creem l'objecte SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(data_matrix)),
  rowData = metabolite_metadata,
  colData = sample_metadata,
  metadata = list(general_metadata = general_metadata$Experiment)
)


## Donat a que ens poden ser útils per després, podem crear dos dades més a partir
## de les dades obtingudes, en primer lloc les podem 'escalar' aplicant logaritmes
## i també calcularem la abundància relativa 'across groups' per a poder 
## replicar un dels boxplots que realitzen

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(data_matrix),escalados = t(log(t(assays(se)$counts)))),
  rowData = metabolite_metadata,
  colData = sample_metadata,
  metadata = list(general_metadata = general_metadata$Experiment)
)

abund_real <- function(data){
  suma<- sum(data)
  abundance_rel <- data/suma*100
  return(c(abundance_rel))
}
med_vals <- apply(t(assays(se)$escalados),2,median)
data_rel <- as.data.frame(sweep(t(assays(se)$escalados),2,med_vals,"-"))
data_rel <- as.data.frame(t(data_rel))

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(data_matrix),escalados = t(log(t(assays(se)$counts)))
                ,log_abund=data_rel),
  rowData = metabolite_metadata,
  colData = sample_metadata,
  metadata = list(general_metadata = general_metadata$Experiment)
)


## -------------------------------------------------------------------------------------------------------------------
se



## ----fig.align='center'---------------------------------------------------------------------------------------------
boxplot(assays(se)$counts,
        las = 2,
        names = colData(se)$`Tipus MS`, 
        cex.axis=0.8,
        main="Distribució dels valors dels metabolits"

)


## ----fig.align='center', echo=FALSE,fig.height=4--------------------------------------------------------------------
boxplot(assays(se)$escalados,
        las = 2,
        names = colData(se)$`Tipus MS`,  
        cex.axis=0.8,
        main="Distribució dels valors dels metabolits escalats Log"

)


## ----fig.align='center',fig.height=3.8------------------------------------------------------------------------------
pcX<-prcomp(t(assays(se)$escalados), scale=FALSE) # Ja s'han escalat les dades
loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)

summary(pcX)$importance[,1:4]
colores <- c("red", "forestgreen")  

xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))
plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colores, 
     main ="Principal components (PCA)")
names2plot<-colnames(se)

text(pcX$x[,1],pcX$x[,2],names2plot, pos=3, cex=.6)


## ----fig.align='center',fig.height=3.5------------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(assays(se)$escalados)),method="average")
plot(clust.euclid.average, hang=-1,,xlab="Dendograma")


## ----fig.align='center',message=FALSE,error=FALSE,warning=FALSE-----------------------------------------------------
library(ggplot2)
library(reshape2) 

boxplot(assays(se)$log_abund,
        las = 2,
        col = colores[as.numeric(colData(se)$`Tipus MS`)], 
        ylab = "Log Abundance",
        ylim = c(-4, 4)      

)

# Añade la línea horizontal en y = 0
abline(h = 0, col = "black", lty = 1, lwd = 1)  


## ----fig.align='center',message=FALSE,error=FALSE,warning=FALSE,echo=FALSE------------------------------------------
ttest=function(x){tt=t.test(x[1:13],x[14:33])
return(c(tt$statistic,
         tt$p.value,
         tt$estimate[2]/tt$estimate[1]))
}
fc=function(x){tt=t.test(x[1:13],x[14:33])
return(c(tt$statistic,
         tt$p.value,
         tt$estimate[1]/tt$estimate[2]))
}

ans <- apply(assays(se)$escalados,1,ttest)
ts <- ans[1,]

pvals<-ans[2,]
pvals_adj <- p.adjust(pvals,method = "BH")
fc <-  apply(assays(se)$counts,1,ttest)

fc2 <- log2(fc[3,])

data_grafic <- data.frame(meta=rownames(se),log2FC=fc2,log10p=(pvals))
library(EnhancedVolcano)

p1 <- EnhancedVolcano(data_grafic,
    lab = data_grafic$meta,  
    x = 'log2FC',  
    y = 'log10p',  
    pCutoff = 0.05,  
    FCcutoff = 0.585, 
    title = 'Volcano Plot',
    selectLab = data_grafic$meta[data_grafic$log10p < 0.05], 
    subtitle = 'Diferencial de Expresión',
    pointSize = 3.0,  
    labSize = 2.5,   
    colAlpha = 0.8,   
    col = c("black", "blue", "red", "orange"),
    ylim = c(0,4),
    xlim=c(-2,2.5),   
    drawConnectors = TRUE, 
    widthConnectors = 0.5,  
    colConnectors = "grey50", 
    boxedLabels = TRUE, 
    max.overlaps = 20  
)
p1


## ----fig.align='center',message=FALSE,error=FALSE,warning=FALSE,echo=FALSE------------------------------------------
data_grafic <- data.frame(meta=rownames(se),log2FC=fc2,log10p=(pvals_adj))

p2 <- EnhancedVolcano(data_grafic,
    lab = data_grafic$meta,  
    x = 'log2FC',  
    y = 'log10p',  
    pCutoff = 0.05,  
    FCcutoff = 0.585, 
    title = 'Volcano Plot',
    selectLab = data_grafic$meta[data_grafic$log10p < 0.05], 
    subtitle = 'Diferencial de Expresión, p adjusted',
    pointSize = 3.0,  
    labSize = 2.5,   
    colAlpha = 0.8,   
    col = c("black", "blue", "red", "orange"),
    ylim = c(0,4),
    xlim=c(-2,2.5),   
    drawConnectors = TRUE, 
    widthConnectors = 0.5,  
    colConnectors = "grey50", 
    boxedLabels = TRUE, 
    max.overlaps = 20  
)
p2

