locpoly <- function(x, y, z,
                    xo = seq(min(x), max(x), length = nx),
                    yo = seq(min(y), max(y), length = ny),
                    nx = 40, ny = 40,
                    input = "points", output = "grid",
                    h = 0, kernel = "uniform", solver = "QR", degree = 3,
                    pd = ""){

    ## secondary use of the partial derivatives estimate for Akimas splines:
    ## use them directly grid- or pointwise.

    if(!(output %in% c("grid","points"))){
        stop("unknown value for \"output\"!")
    }
    if(!(input %in% c("grid","points"))){
        stop("unknown value for \"output\"!")
    }

    if(input=="grid"){
        nx <- length(x)
        ny <- length(y)
        if(dim(z)[1]!=nx | dim(z)[2]!=ny)
            stop("wrong dimensions of x, y, and z!")
        x <- matrix(rep(x,ny),nx,ny)
        y <- t(matrix(rep(y,nx),ny,nx))
    }
    
    if(pd=="all"){
        if(output=="grid"){
            ans <- locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)
            ans$x=xo
            ans$y=yo
        } else {
            ans <- locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)
            ans$x=xo
            ans$y=yo
        }
    } else if(pd==""){
        if(output=="grid")
            ans <- list(x=xo,y=yo,z=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$z)
        else
            ans <- list(x=xo,y=yo,z=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$z)

    } else if(pd=="x"){
        if(degree>0){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zx=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zx)
            else
                ans <- list(x=xo,y=yo,zx=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zx)
        } else stop("need degree>0 for pd=\"x\"")
    } else if(pd=="y"){
        if(degree>0){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zy=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zy)
            else
                ans <- list(x=xo,y=yo,zy=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zy)
        } else stop("need degree>0 for pd=\"y\"")
    } else if(pd=="xx"){
        if(degree>1){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zxx=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zxx)
            else
                ans <- list(x=xo,y=yo,zxx=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zxx)
        } else stop("need degree>1 for pd=\"xx\"")
    } else if(pd=="yy"){
        if(degree>1){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zyy=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zyy)
            else
                ans <- list(x=xo,y=yo,zyy=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zyy)
        } else stop("need degree>1 for pd=\"yy\"")
    } else if(pd=="xy"){
        if(degree>1){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zxy=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zxy)
            else
                ans <- list(x=xo,y=yo,zxy=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zxy)
        } else  stop("need degree>1 for pd=\"xy\"")
    } else if(pd=="xxx"){
        if(degree>2){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zxxx=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zxxx)
            else
                ans <- list(x=xo,y=yo,zxxx=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zxxx)
        } else  stop("need degree>2 for pd=\"xxx\"")
    } else if(pd=="yyy"){
        if(degree>2){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zyyy=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zyyy)
            else
                ans <- list(x=xo,y=yo,zyyy=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zyyy)
        } else stop("need degree>2 for pd=\"yyy\"")
    } else if(pd=="xxy"){
        if(degree>2){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zxxy=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zxxy)
            else
                ans <- list(x=xo,y=yo,zxxy=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zxxy)
        } else stop("need degree>2 for pd=\"xxy\"")
    } else if(pd=="xyy"){
        if(degree>2){
            if(output=="grid")
                ans <- list(x=xo,y=yo,zxyy=locpoly.partderiv.grid(xo,yo,x,y,z,kernel,h,solver,degree)$zxyy)
            else
                ans <- list(x=xo,y=yo,zxyy=locpoly.partderiv.points(xo,yo,x,y,z,kernel,h,solver,degree)$zxyy)
        } else stop("need degree>2 for pd=\"xyy\"")
    } else
        stop(paste("unsupported value for pd: ", pd,
                   "\nonly partial derivatives of order up to 3 can be estimated (needs degree=3)!"))

    ans
}
