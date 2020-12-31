randomProfiles = function(Freqs, BlockSize = 1000){
    Profile = vector(mode = "list", length = BlockSize)
    profileVec = .randomProfiles(Freqs$freq, BlockSize)
    nLoci = nLoci = length(Freqs$loci)
    
    for(b in 1:BlockSize){
        i1 = (b - 1)*2*nLoci + 1
        i2 =  b*2*nLoci

        Profile[[b]] = matrix(profileVec[i1:i2], ncol = 2, nrow = nLoci, byrow = T)
        class(Profile[[b]]) = "profile"
    }

    return(Profile)
}
