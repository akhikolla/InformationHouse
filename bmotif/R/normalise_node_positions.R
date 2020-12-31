normalise_node_positions <- function(pc,type,six_node){
  if(!ncol(pc) %in% c(46, 148)){stop("Something has gone very wrong: pc does not have 46 or 148 columns")}
  if(six_node == TRUE & ncol(pc) != 148){stop("Something has gone very wrong: six_node is TRUE but ncol != 148. Contact the package maintainer.")}
  if(six_node == FALSE & ncol(pc) != 46){stop("Something has gone very wrong: six_node is FALSE but ncol != 46. Contact the package maintainer.")}
  if(type == "none"){
    return(pc)
  } else if(type == "sizeclass"){
    pc[,paste0("np",1:2)] <- pc[,paste0("np",1:2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1:2)]))
    pc[,paste0("np",3:6)] <- pc[,paste0("np",3:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3:6)]))
    pc[,paste0("np",7:16)] <- pc[,paste0("np",7:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7:16)]))
    pc[,paste0("np",17:46)] <- pc[,paste0("np",17:46)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17:46)]))
    if(six_node){
      pc[,paste0("np",47:148)] <- pc[,paste0("np",47:148)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",47:148)]))
    }
    pc[do.call(cbind, lapply(pc, is.nan))] <- NA
  } else if(type == "sizeclass_plus1"){
    pc <- pc + 1
    pc[,paste0("np",1:2)] <- pc[,paste0("np",1:2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1:2)]))
    pc[,paste0("np",3:6)] <- pc[,paste0("np",3:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3:6)]))
    pc[,paste0("np",7:16)] <- pc[,paste0("np",7:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7:16)]))
    pc[,paste0("np",17:46)] <- pc[,paste0("np",17:46)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17:46)]))
    if(six_node){
      pc[,paste0("np",47:148)] <- pc[,paste0("np",47:148)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",47:148)]))
    }
  } else if (type == "sizeclass_NAzero"){
    pc[,paste0("np",1:2)] <- pc[,paste0("np",1:2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1:2)]))
    pc[,paste0("np",3:6)] <- pc[,paste0("np",3:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3:6)]))
    pc[,paste0("np",7:16)] <- pc[,paste0("np",7:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7:16)]))
    pc[,paste0("np",17:46)] <- pc[,paste0("np",17:46)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17:46)]))
    if(six_node){
      pc[,paste0("np",47:148)] <- pc[,paste0("np",47:148)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",47:148)]))
    }
    pc[is.na(pc)] <- 0
  } else if(type == "sum"){
    # pc[,paste0("np",1:ncol(pc))] <- pc[,paste0("np",1:ncol(pc))]/sapply(1:nrow(pc), function(x) sum(pc[x,]))
    pc <- pc/apply(pc, 1, sum)
  } else if(type == "position"){
    pc <- apply(pc, MARGIN = 2, FUN = function(x) x/sum(x))
    pc[do.call(cbind, lapply(pc, is.nan))] <- NA
  } else if(type == "levelsize"){
    pc[,paste0("np",1:2)] <- pc[,paste0("np",1:2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1:2)]))
    pc[,paste0("np",3:4)] <- pc[,paste0("np",3:4)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3:4)]))
    pc[,paste0("np",5:6)] <- pc[,paste0("np",5:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",5:6)]))
    pc[,paste0("np",7:8)] <- pc[,paste0("np",7:8)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7:8)]))
    pc[,paste0("np",9:14)] <- pc[,paste0("np",9:14)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",9:14)]))
    pc[,paste0("np",15:16)] <- pc[,paste0("np",15:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",15:16)]))
    pc[,paste0("np",17:18)] <- pc[,paste0("np",17:18)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17:18)]))
    pc[,paste0("np",19:31)] <- pc[,paste0("np",19:31)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",19:31)]))
    pc[,paste0("np",32:44)] <- pc[,paste0("np",32:44)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",32:44)]))
    pc[,paste0("np",45:46)] <- pc[,paste0("np",45:46)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",45:46)]))
    if(six_node){
      pc[,paste0("np",47:48)] <- pc[,paste0("np",47:48)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",47:48)]))
      pc[,paste0("np",49:70)] <- pc[,paste0("np",49:70)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",49:70)]))
      pc[,paste0("np",71:124)] <- pc[,paste0("np",71:124)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",71:124)]))
      pc[,paste0("np",125:146)] <- pc[,paste0("np",125:146)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",125:146)]))
      pc[,paste0("np",147:148)] <- pc[,paste0("np",147:148)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",147:148)]))
    }
    pc[do.call(cbind, lapply(pc, is.nan))] <- NA
  } else if(type == "levelsize_plus1") {
    pc <- pc + 1
    pc[,paste0("np",1:2)] <- pc[,paste0("np",1:2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1:2)]))
    pc[,paste0("np",3:4)] <- pc[,paste0("np",3:4)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3:4)]))
    pc[,paste0("np",5:6)] <- pc[,paste0("np",5:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",5:6)]))
    pc[,paste0("np",7:8)] <- pc[,paste0("np",7:8)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7:8)]))
    pc[,paste0("np",9:14)] <- pc[,paste0("np",9:14)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",9:14)]))
    pc[,paste0("np",15:16)] <- pc[,paste0("np",15:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",15:16)]))
    pc[,paste0("np",17:18)] <- pc[,paste0("np",17:18)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17:18)]))
    pc[,paste0("np",19:31)] <- pc[,paste0("np",19:31)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",19:31)]))
    pc[,paste0("np",32:44)] <- pc[,paste0("np",32:44)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",32:44)]))
    pc[,paste0("np",45:46)] <- pc[,paste0("np",45:46)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",45:46)]))
    if(six_node){
      pc[,paste0("np",47:48)] <- pc[,paste0("np",47:48)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",47:48)]))
      pc[,paste0("np",49:70)] <- pc[,paste0("np",49:70)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",49:70)]))
      pc[,paste0("np",71:124)] <- pc[,paste0("np",71:124)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",71:124)]))
      pc[,paste0("np",125:146)] <- pc[,paste0("np",125:146)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",125:146)]))
      pc[,paste0("np",147:148)] <- pc[,paste0("np",147:148)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",147:148)]))
      }
    } else if(type == "levelsize_NAzero"){
      pc[,paste0("np",1:2)] <- pc[,paste0("np",1:2)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",1:2)]))
      pc[,paste0("np",3:4)] <- pc[,paste0("np",3:4)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",3:4)]))
      pc[,paste0("np",5:6)] <- pc[,paste0("np",5:6)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",5:6)]))
      pc[,paste0("np",7:8)] <- pc[,paste0("np",7:8)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",7:8)]))
      pc[,paste0("np",9:14)] <- pc[,paste0("np",9:14)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",9:14)]))
      pc[,paste0("np",15:16)] <- pc[,paste0("np",15:16)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",15:16)]))
      pc[,paste0("np",17:18)] <- pc[,paste0("np",17:18)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",17:18)]))
      pc[,paste0("np",19:31)] <- pc[,paste0("np",19:31)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",19:31)]))
      pc[,paste0("np",32:44)] <- pc[,paste0("np",32:44)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",32:44)]))
      pc[,paste0("np",45:46)] <- pc[,paste0("np",45:46)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",45:46)]))
      if(six_node){
        pc[,paste0("np",47:48)] <- pc[,paste0("np",47:48)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",47:48)]))
        pc[,paste0("np",49:70)] <- pc[,paste0("np",49:70)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",49:70)]))
        pc[,paste0("np",71:124)] <- pc[,paste0("np",71:124)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",71:124)]))
        pc[,paste0("np",125:146)] <- pc[,paste0("np",125:146)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",125:146)]))
        pc[,paste0("np",147:148)] <- pc[,paste0("np",147:148)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",147:148)]))
      }
      pc[is.na(pc)] <- 0
    } else if(type == "motif"){
      nm <- ifelse(six_node, yes = 44, no = 17) # number of motifs
      mp <- sapply(1:nm, function(x) motif_info(x, link = FALSE)) # motif positions
      for(i in mp){
        pc[,paste0("np",i)] <- pc[,paste0("np",i)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",i)]))
      }
      pc[do.call(cbind, lapply(pc, is.nan))] <- NA
    } else if(type == "motif_plus1"){
      pc <- pc + 1
      nm <- ifelse(six_node, yes = 44, no = 17) # number of motifs
      mp <- sapply(1:nm, function(x) motif_info(x, link = FALSE)) # motif positions
      for(i in mp){
        pc[,paste0("np",i)] <- pc[,paste0("np",i)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",i)]))
      }
    } else if(type == "motif_NAzero"){
      nm <- ifelse(six_node, yes = 44, no = 17) # number of motifs
      mp <- sapply(1:nm, function(x) motif_info(x, link = FALSE)) # motif positions
      for(i in mp){
        pc[,paste0("np",i)] <- pc[,paste0("np",i)]/sapply(1:nrow(pc), function(x) sum(pc[x,paste0("np",i)]))
      }
      pc[is.na(pc)] <- 0
    } else {
    stop("'type' must be a character string equal to 'sum', 'sizeclass', 'position' or 'levelsize'") # if 'type' does not equal 'sum' or 'sizeclass' return an error
  }
  return(pc)
}
