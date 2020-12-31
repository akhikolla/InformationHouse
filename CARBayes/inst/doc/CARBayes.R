### R code from vignette source 'CARBayes.Rnw'

###################################################
### code chunk number 1: CARBayes.Rnw:63-64
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: CARBayes.Rnw:403-409
###################################################
library(CARBayesdata)
library(shapefiles)
library(sp)
data(lipdata)
data(lipdbf)
data(lipshp)


###################################################
### code chunk number 3: CARBayes.Rnw:414-417
###################################################
library(CARBayes)
lipdbf$dbf <- lipdbf$dbf[ ,c(2,1)]
data.combined <- combine.data.shapefile(data=lipdata, shp=lipshp, dbf=lipdbf)


###################################################
### code chunk number 4: CARBayes.Rnw:431-436
###################################################
library(CARBayesdata)
library(sp)
data(GGHB.IG)
data(pricedata)
head(pricedata)


###################################################
### code chunk number 5: CARBayes.Rnw:460-462
###################################################
library(dplyr)
pricedata <- pricedata %>% mutate(logprice = log(pricedata$price))


###################################################
### code chunk number 6: CARBayes.Rnw:468-470
###################################################
library(GGally)
ggpairs(data = pricedata, columns = c(8, 3:7))


###################################################
### code chunk number 7: CARBayes.Rnw:483-484
###################################################
pricedata.sp  <- merge(x=GGHB.IG, y=pricedata, by="IG", all.x=FALSE)


###################################################
### code chunk number 8: CARBayes.Rnw:492-495
###################################################
library(rgdal)
pricedata.sp  <- spTransform(pricedata.sp , 
                    CRS("+proj=longlat +datum=WGS84 +no_defs"))


###################################################
### code chunk number 9: CARBayes.Rnw:500-510
###################################################
library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = pricedata.sp@data$price)
map1 <- leaflet(data=pricedata.sp) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(price), color="", weight=1, 
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = pricedata.sp@data$price, opacity = 1, 
                title="Price") %>%
    addScaleBar(position="bottomleft")
map1


###################################################
### code chunk number 10: CARBayes.Rnw:526-529
###################################################
form <- logprice~crime+rooms+sales+factor(type) + driveshop
model <- lm(formula=form, data=pricedata.sp@data)
summary(model)


###################################################
### code chunk number 11: CARBayes.Rnw:534-538
###################################################
library(spdep)
W.nb <- poly2nb(pricedata.sp, row.names = rownames(pricedata.sp@data))
W.list <- nb2listw(W.nb, style="B")
moran.mc(x=residuals(model), listw=W.list, nsim=1000)


###################################################
### code chunk number 12: CARBayes.Rnw:548-549
###################################################
W <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 13: CARBayes.Rnw:731-736
###################################################
library(CARBayesdata)
library(sp)
data(GGHB.IG)
data(respiratorydata)
head(respiratorydata)


###################################################
### code chunk number 14: CARBayes.Rnw:742-743
###################################################
respiratorydata.sp <- merge(x=GGHB.IG, y=respiratorydata, by="IG", all.x=FALSE)


###################################################
### code chunk number 15: CARBayes.Rnw:748-750
###################################################
respiratorydata.sp <- spTransform(respiratorydata.sp, 
                          CRS("+proj=longlat +datum=WGS84 +no_defs"))


###################################################
### code chunk number 16: CARBayes.Rnw:757-767
###################################################
library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = respiratorydata.sp@data$SMR)
map2 <- leaflet(data=respiratorydata.sp) %>% 
  addTiles() %>% 
  addPolygons(fillColor = ~colours(SMR), color="", weight=1, 
              fillOpacity = 0.7) %>%
  addLegend(pal = colours, values = respiratorydata.sp@data$SMR, opacity = 1, 
            title="SMR") %>%
  addScaleBar(position="bottomleft")
map2


###################################################
### code chunk number 17: CARBayes.Rnw:782-784
###################################################
W.nb <- poly2nb(respiratorydata.sp, row.names = rownames(respiratorydata.sp@data))
W <- nb2mat(W.nb, style="B")


###################################################
### code chunk number 18: CARBayes.Rnw:799-801
###################################################
income <- respiratorydata.sp@data$incomedep
Z.incomedep <- as.matrix(dist(income, diag=TRUE, upper=TRUE)) 


