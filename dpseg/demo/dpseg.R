library("dpseg")

## NOTE: library loads bacterial growth curve data as data.frame oddata

## calculate linear segments in semi-log bacterial growth data
Sj <- dpseg(x=oddata$Time, y=log(oddata$A1), minl=5, P=0.0001, verb=1, store.matrix=TRUE)

## inspect resulting segments
print(Sj)

## plot results
plot(Sj, delog=TRUE, log="y")

## predict method
plot(predict(Sj), type="l") 

## view the algorithm in action, use file.name to save as <file.name>.mpeg
## NOTE: store.matrix in the dpseg call above stores the scoring matrix,
## this is required for nicer educational movie
par(ask=FALSE) # turns off asking to `Hit <Return> to see the next plot`
movie(Sj, delay=0.1)


