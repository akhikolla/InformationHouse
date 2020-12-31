#' Read in a file of allele frequencies
#' 
#' Reads in a file of alleles in a particular format.
#' 
#' This function reads frequencies in the rectangular allele freqency table
#' format used by FSI Genetics and other journals.  This file format assumes a
#' comma separated value file (CSV) (although the column delimeter can be
#' specified). The first column should be labelled 'Allele' and contain the STR
#' allele designations that are used in the data set. The remaining columns
#' will have the locus name as a header, and frequencies that are either blank,
#' zero, or non-zero.  Blanks or zeros are used to specify that the allele is
#' not observed (and not used) at the locus. The final row of the file should
#' start with 'N' or 'n' in the first column and give the number of individuals
#' typed (or the number of alleles recorded) in assessing the frequency of the
#' alleles.
#' 
#' The second format is a very particular 'Curran' text format. The first line
#' contains the number of loci in the multiplex. The next line will contain the
#' name of the first locus and the number of alleles, nA, the locus separated
#' by a comma. The next nA lines contain the allele number (from 1 to nA), the
#' STR designation of the allele, and the frequency separated by commas. This
#' pattern is repeated for each locus. In the future this function will read
#' the rectangular allele freqency table used by FSI Genetics and other
#' journals.
#' 
#' @param strPath The file from which to read the frequencies
#' @param FSIGenFormat Tells the function whether the file is either in FSI
#' Genetics format (see below) or 'Curran' format
#' @param delim This argument is used when \code{FSIGenFormat} is \code{TRUE},
#' and is the regular expression used to delimit columns of the table. it is
#' set to a single comma by default, and multiple delimiters are considered
#' empty separate fields. There probably should be an additional argument which
#' specifies the missing or empty cell symbol, but I won't programme this
#' unless somebody asks for it
#' @return a list containing two vectors and a list, loci, counts, and freqs.
#' The vector loci is a vector of the locus names in the frequency file. The
#' vector counts is a vector of the number of individuals (or sometimes
#' alleles) typed at each locus. This will null if the 'Curran' format is used.
#' The list freqs, is a list of vectors with each vector containing the
#' frequencies of the alleles at the locus. The names of the elements of the
#' vectors are the STR allele designations.
#' @author James M. Curran
#' @export readFreqs
readFreqs = function(strPath, FSIGenFormat = TRUE, delim = ','){
    Freqs = list(loci = NULL, freqs = NULL, counts = NULL)

    if(FSIGenFormat){
        ## expect the information to start with Alleles, and end with n
        f1 = file(strPath, 'r')

        if(isOpen(f1)){
            Lines = readLines(f1)
            close(f1)

            alleleStart = -1

            ## look for a line starting with alleles

            if(!any(grepl('^[Aa]llele', Lines))){
                stop("The header line should start with 'Allele'")
            }else{
                alleleStart = grep('^[Aa]llele', Lines)
                if(length(alleleStart) > 1){
                    stop("There are two or more lines starting with Allele")
                }
            }


            ## look for line starting with n
            countStart = -1

            if(!any(grepl('^[Nn]', Lines))){
                stop("The frequency file should have locus allele counts")
            }else{
                countStart = grep('^[Nn]', Lines)
                if(length(countStart) > 1){
                    warning("There should only one locus allele count line. Using first ")
                    countStart = countStart[1]
                }
            }

            if(countStart == -1 || alleleStart == -1){
                stop("Header and footer lines not found")
            }else{
                n = length(Lines)
                Lines = Lines[alleleStart:countStart]

                if((n - length(Lines)) > 0){
                    cat(paste("Dropped", n - length(Lines), "lines\n"))
                }

                locusLine = Lines[1]
                Tokens = unlist(strsplit(locusLine, delim))
                Loci = toupper(Tokens[-1])
                nLoci = length(Loci)
                cat(paste("Detected", nLoci, "loci\n"))

                countLine = Lines[length(Lines)]
                Tokens = unlist(strsplit(countLine, delim))
                counts = as.numeric(Tokens[-1])

                Lines = Lines[-c(1,length(Lines))]

                freqs = vector(length = nLoci, mode = "list")
                names(freqs) = Loci

                for(line in Lines){
                    Tokens = unlist(strsplit(line, delim))

                    a = Tokens[1]
                    fx = as.numeric(Tokens[-1])


                    for(loc in 1:nLoci){
                        if(!is.na(fx[loc]) && fx[loc] > 0){
                            f = fx[loc]

                            if(f > 1){
                                stop("Error: frequency greater than 1")
                            }

                            freqs[[loc]] = c(freqs[[loc]], f)
                            names(freqs[[loc]])[length(freqs[[loc]])] = a
                        }
                    }
                }

                cat(paste("Successfully processed", length(Lines), "lines\n"))

                Freqs$freqs = freqs
                Freqs$loci = Loci
                Freqs$counts = counts
            }
        }else{
            stop(paste("Couldn't open", strPath, "for reading"))
        }
    }else{ ## Curran format
        f1 = file(strPath, "r")

        if(isOpen(f1)){
            Lines = readLines(f1)
            nLines = length(Lines)
            ## cat(paste("Read", nLines, "lines from", strPath,"\n"))
            close(f1)
        }else{
            stop(paste("Couldn't open :", strPath, "for reading"))
        }

        Freqs = NULL

        nLoci = as.numeric(Lines[1])
        ## cat(paste(nLoci, "\n"))

        Loci = rep("", nLoci)
        freqs = vector(nLoci, mode = "list")
        currLine = 2

        for(nLoc in 1:nLoci){
            Tokens = unlist(strsplit(Lines[currLine], ","))
            currLine = currLine + 1

            Loci[nLoc] = Tokens[1]
            nAlleles = as.numeric(Tokens[2])
            Alleles = rep("", nAlleles)

            freqs[[nLoc]] = rep(0, nAlleles)

            for(nA in 1:nAlleles){
                Tokens = unlist(strsplit(Lines[currLine], ","))
                currLine = currLine + 1
                Alleles[nA] = Tokens[2]
                freqs[[nLoc]][nA] = as.numeric(Tokens[3])
            }

            names(freqs[[nLoc]]) = Alleles
        }

        names(freqs) = Loci
        Freqs$loci = Loci
        Freqs$freqs = freqs
    }
    return(Freqs)
}
