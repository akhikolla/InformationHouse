#' Retrieve data from Budowle and Moretti (1999) from the web
#' 
#' Retreives the Budowle and Moretti (1999) and compiles the allele frequency
#' tables needed for the other parts of this package such as \code{sim}.
#' 
#' The first three populations have data on 20 loci, the second three on 13
#' loci. The missing values (0's in the raw data) have been dropped and are not
#' used in calculating the frequencies. This function will not work if you are
#' not connected to the internet, or access to the internet is blocked.
#' 
#' @param url - the location of the webpage this data is stored on. If
#'   \code{NULL} then hardcoded values in the function are used. The argument
#'   allows the user to get the function to work given this values may change.
#'   However, it is unlikely that it will help.
#' @param id - the id of the HTML element where the data is stored on the
#'   webpage. If \code{NULL} then a hardcoded value in the function is used. The
#'   argument allows the user to get the function to work given this values may
#'   change. However, it is unlikely that it will help. NOTE: this is super
#'   flakey because of how the FBI web authors have created it. This may change
#'   in which case, the function will more likely change than the id.
#'   
#' @return A list consisting of six elements corresponding to the six
#' populations detailed in the data set. Each of the list elements is a list in
#' itself with three further elements named \code{loci}, \code{profiles} and \code{freqs}. 
#' \code{loci} is a vector of the 13-20 STR locus names. \code{freqs} is a list of 13-20 
#' vectors, each vector contains the allele frequencies. \code{profiles} contains the 
#' raw profiles that the allele frequency tables were constructed from.
#' @author James M. Curran
#' @seealso fbiCaucs, USCaucs
#' @references Budowle, B. and Moretti, T.R. (1999), \emph{Genotype Profiles
#' for Six Population Groups at the 13 CODIS Short Tandem Repeat Core Loci and
#' Other PCR Based Loci}, Forensic Science Communications 1(2).
#' @keywords datasets
#' @examples
#' 
#' ## not run
#' \dontrun{
#' db = fetchBMdata()
#' names(db)
#' f = db[["TRINIDADIAN"]]$freqs
#' dbExpect(f, k = "UN", collapse = TRUE)
#' }
#' 
#' @importFrom stats xtabs
#' @importFrom xml2 read_html 
#' @importFrom rvest html_nodes html_text
#' @importFrom stringr str_trim
#' @export fetchBMdata
fetchBMdata = function(url = NULL, id = NULL){
    
    if(is.null(url)){
        url = "http://www.fbi.gov/about-us/lab/forensic-science-communications/fsc/july1999/dnaloci.txt"
    }
    if(is.null(id)){
        id = "parent-fieldname-text-2baa15a4d877814a4a62588114966028"
    }
    
    ## Note this is the flakiest part of this function. If the FBI decide
    ## to move the data again, or reformat the website this will almost certainly
    ## fail. Your tax dollars at work folks!
    fetchLines = function(){
        webpage = read_html(url)
        xpath = paste0('//div[@id="', id, '"]')
        data = html_text(html_nodes(webpage, xpath = xpath))
        Lines = unlist(strsplit(data, "[\r]+"))
        i1 = grep("Table 1.", Lines)
        i2 = grep("^3409", Lines)
        Lines = Lines[i1:i2]
    }
    
    Lines = fetchLines()
    Lines = gsub("[&][#]13;", "", Lines)
    Lines = stringr::str_trim(Lines)

    ## All rares are encoded 108.1
    Lines = gsub("(<|>)[0-9]+","108.1",Lines)

    ## drop the empty lines

    Lines = Lines[nchar(Lines)>0]

    ## drop the first line

    Lines = Lines[-1]

    ## The lines that have length < 50 are the population names

    popLocations = which(nchar(Lines) < 50)
    popNames = Lines[popLocations]

    ## Locus names are the same in all six populations
    ## However, last 3 pops only have CODIS loci
    ## First tag on locus line is ID# so drop that

    Loci = unlist(strsplit(Lines[2], "[\t ]+"))[-1]
    nLoci = length(Loci)

    ## get the alleles possible at each locus

    Alleles = vector(length = nLoci, mode = "list")

    allPops = strsplit(Lines[-c(popLocations,popLocations+1)],"[\t ]+")
    
    allPops = lapply(allPops, function(profile){
                            l = 2 * nLoci + 1 - length(profile)
                            if(l != 0){
                                profile = c(profile, rep(NA, l))
                            }
                            profile[profile == "0"] = NA
                            profile
                         }    
    )
    
    Alleles = vector(mode = "list", length = nLoci)
    apm = do.call("rbind", allPops)
    
    for(loc in 1:nLoci){
        i1 = 2 * loc
        i2 = i1 + 1
        alleles = c(apm[,i1], apm[,i2])
        alleles = alleles[!is.na(alleles)]
        Alleles[[loc]] = sort(unique(alleles))
        
        if(all(grepl("^[0-9.]+$", Alleles[[loc]]))){
            Alleles[[loc]] = as.character(sort(as.numeric(Alleles[[loc]])))
        }
    }
    names(Alleles) = Loci
    
    numPops = length(popNames)
    db = vector(length = numPops, mode = "list")
    names(db) = popNames

    i1 = popLocations + 2
    i2 = c(popLocations[-1] - 1, length(Lines))
    popSizes = i2 - i1 + 1
    start = cumsum(c(1, popSizes[-numPops]))
    end = cumsum(popSizes)


    for(pop in 1:numPops){
        if(pop > 3){
            nLoci = 13 ## the last three pops only have 13 loci
            Loci = Loci[1:nLoci]
        }

        db[[pop]] = list(loci = Loci,
                         profiles = apm[start[pop]:end[pop], 1:(2 * nLoci + 1)],
                         freqs = vector(mode = "list", length = nLoci))

        names(db[[pop]]$freqs) = Loci

        #browser()
        for(loc in 1:nLoci){
            i1 = 2 * loc
            i2 = i1 + 1
            a = c(db[[pop]]$profiles[,i1], db[[pop]]$profiles[,i2])
            a = factor(a, levels = Alleles[[loc]])
            tbl = xtabs(~a)
            db[[pop]]$freqs[[loc]] = tbl / sum(tbl)
        }
    }

        # popLines = strsplit(Lines[i1[pop]:i2[pop]],"[\t ]+")
        # popM = matrix(0, nrow = length(popLines), ncol = 2*nLoci)
        # 
        # 
        # j1 = which(sapply(popLines, length) == 2*nLoci + 1)
        # 
        # for(j in j1){
        #     popM[j,] = popLines[[j]][-1]
        # }
        # 
        # j1 = (1:nrow(popM))[-j1]
        # 
        # for(i in j1){
        #     popM[j,1:26] = popLines[[j]][-1]
        # }

    #     for(loc in 1:nLoci){
    #         j1 = 2*(loc - 1)  + 1
    #         j2 = 2*loc
    # 
    #         tbl = NULL
    # 
    #         ## 1st 14 loci have numeric alleles
    #         if(loc <=14){
    #             tbl = table(as.numeric(popM[,j1:j2]))
    #             A = as.numeric(names(tbl))
    #             tbl = tbl[A!=0]
    #             A = A[A!=0]
    # 
    #             tbl = tbl/sum(tbl)
    # 
    #             pos = match(A, Alleles[[loc]])
    # 
    #             db[[pop]]$freqs[[loc]][pos] = tbl
    #         }else{
    #             tbl = table(popM[,j1:j2])
    #             A = names(tbl)
    #             tbl = tbl[A!="0"]
    #             A = A[A!="0"]
    # 
    #             tbl = tbl/sum(tbl)
    # 
    #             pos = match(A, Alleles[[loc]])
    # 
    #             db[[pop]]$freqs[[loc]][pos] = tbl
    #         }
    #     }
    # }

    
    # allPopsM = matrix(0, ncol = 2*nLoci, nrow = length(allPops))
    # browser()
    # 
    # i1 = which(sapply(allPops, length) == 2*nLoci + 1)
    # 
    # for(i in i1){
    #     allPopsM[i,] = allPops[[i]][-1]
    # }
    # 
    # i1 = (1:nrow(allPopsM))[-i1]
    # browser()
    # 
    # for(i in i1){
    #     allPopsM[i,1:26] = allPops[[i]][-1]
    # }
    # 
    # allPops = allPopsM
    # 
    # for(loc in 1:nLoci){
    #     i1 = 2*(loc - 1)  + 1
    #     i2 = 2*loc
    # 
    #     tbl = table(allPops[,i1:i2])
    #     A = names(tbl)
    #     A = A[A!="0"]
    # 
    #     if(!any(grepl("[^0-9.]", A))){
    #         A = sort(as.numeric(A))
    #     }
    # 
    #     Alleles[[loc]] = A
    # }
    # 
    # names(Alleles) = Loci
    # 
    # db = vector(length = 6, mode = "list")
    # names(db) = popNames
    # 
    # i1 = popLocations + 2
    # i2 = c(popLocations[-1] - 1, length(Lines))
    # 
    # for(pop in 1:6){
    #     if(pop > 3){
    #         nLoci = 13 ## the last three pops only have 13 loci
    #         Loci = Loci[1:nLoci]
    #     }
    # 
    #     db[[pop]] = list(loci = Loci,
    #                      freqs = vector(mode = "list", length = nLoci))
    # 
    #     names(db[[pop]]$freqs) = Loci
    #     #browser()
    #     
    #     for(loc in 1:nLoci){
    #         db[[pop]]$freqs[[loc]] = rep(0, length(Alleles[[loc]]))
    #         names(db[[pop]]$freqs[[loc]]) = Alleles[[loc]]
    #     }
    # 
    #     popLines = strsplit(Lines[i1[pop]:i2[pop]],"[\t ]+")
    #     popM = matrix(0, nrow = length(popLines), ncol = 2*nLoci)
    # 
    # 
    #     j1 = which(sapply(popLines, length) == 2*nLoci + 1)
    # 
    #     for(j in j1){
    #         popM[j,] = popLines[[j]][-1]
    #     }
    # 
    #     j1 = (1:nrow(popM))[-j1]
    # 
    #     for(i in j1){
    #         popM[j,1:26] = popLines[[j]][-1]
    #     }
    # 
    #     for(loc in 1:nLoci){
    #         j1 = 2*(loc - 1)  + 1
    #         j2 = 2*loc
    # 
    #         tbl = NULL
    # 
    #         ## 1st 14 loci have numeric alleles
    #         if(loc <=14){
    #             tbl = table(as.numeric(popM[,j1:j2]))
    #             A = as.numeric(names(tbl))
    #             tbl = tbl[A!=0]
    #             A = A[A!=0]
    # 
    #             tbl = tbl/sum(tbl)
    # 
    #             pos = match(A, Alleles[[loc]])
    # 
    #             db[[pop]]$freqs[[loc]][pos] = tbl
    #         }else{
    #             tbl = table(popM[,j1:j2])
    #             A = names(tbl)
    #             tbl = tbl[A!="0"]
    #             A = A[A!="0"]
    # 
    #             tbl = tbl/sum(tbl)
    # 
    #             pos = match(A, Alleles[[loc]])
    # 
    #             db[[pop]]$freqs[[loc]][pos] = tbl
    #         }
    #     }
    # }

    return(db)
}
