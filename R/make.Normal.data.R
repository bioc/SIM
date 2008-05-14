`make.Normal.data` <-
function(data)
{
col <- ncol(data)
splits <- round(col/3)

Norm1 <- 1:splits
Norm2 <- (splits+1):(splits*2)
Norm3 <- (splits*2+1):col

NormalData1 <- apply(data[,Norm1], 1, median)
NormalData2 <- apply(data[,Norm2], 1, median)
NormalData3 <- apply(data[,Norm3], 1, median)
NormalData <- cbind(NormalData1, NormalData2, NormalData3)
return(NormalData)
}

