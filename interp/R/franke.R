franke.fn <- function(x,y,fn=1){
    switch(fn,
           "1"=0.75*exp(-((9*x-2)^2+(9*y-2)^2)/4)+
               0.75*exp(-((9*x+1)^2)/49-(9*y+1)/10)+
               0.5*exp(-((9*x-7)^2+(9*y-3)^2)/4)-
               0.2*exp(-(9*x-4)^2-(9*y-7)^2),
           "2"=(tanh(9*y-9*x)+1)/9,
           "3"=(1.25+cos(5.4*y))/(6*(1+(3*x-1)^2)),
           "4"=exp(-81*((x-0.5)^2+(y-0.5)^2)/16)/3,
           "5"=exp(-81*((x-0.5)^2+(y-0.5)^2)/4)/3,
           "6"=sqrt(64-81*((x-0.5)^2+(y-0.5)^2))/9-0.5)
}

franke.data <- function(fn=1,ds=1,data){
    ret <- cbind(x=data[[ds]]$x,y=data[[ds]]$y,
                 z=franke.fn(data[[ds]]$x,y=data[[ds]]$y,fn))
    list(x=ret[,"x"],y=ret[,"y"],z=ret[,"z"])
}
