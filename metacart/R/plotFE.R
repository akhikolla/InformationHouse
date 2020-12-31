#' Visualisation of a FE meta-tree
#'
#' Plot function for a \code{FEmrt} object. The plot shows the result of \code{FEmrt}.
#' The plot function uses the plot method from the package \pkg{ggplot2}
#'
#' For categorical variables we recommend to use short names for levels to avoid overlapping
#' labels at split points.
#' @method plot FEmrt
#' @param x A FEmrt object.
#' @param ... additional arguments to pass
#' @import ggplot2
#' @import gridExtra
#' @export
plot.FEmrt <- function(x, ...){
  if (length(x$n) < 2) {
    warning("No tree was detected")
  } else{
  # Extract necessary information from the FEmrt object

  frame <- x$tree$frame
  # re-label the nodes to eliminate the gap of numbers
  frame0 <- frame
  rownames(frame0) <- rank(as.numeric(rownames(frame)))
  # a vector to translate rpart node labels to our own labels
  trans.labels <- rownames(frame)
  names(trans.labels) <- rownames(frame0)

  # subgrouping the studies by the final tree
  term.nodes <- rownames(frame0)[x$tree$where]
  # find the parent leaves and the order of split
  inx.pleaves <- which(frame0$var != "<leaf>")
  pleaves <- rownames(frame0)[inx.pleaves]
  # order the splits by the first child leaf of the parent leaves
  pleaves <- pleaves[order(as.numeric(rownames(frame0)[inx.pleaves+1]))]
  splits.all <- labels(x$tree, minlength = 0L)
  splits <- splits.all[inx.pleaves+1][order(as.numeric(rownames(frame0)[inx.pleaves+1]))]
  ntree <- data.frame(split = as.character(splits), pleaf = as.numeric(pleaves))
  # create an ojbect with all plotting information
  object <- list()
  object$data <- x$data
  object$data$term.node <- term.nodes
  object$tree <- ntree
  object$n <- table(term.nodes)
  object$tree$split <- as.character(object$tree$split)

  # plot
  if (length(object$n) < 2) {stop("no tree was detected")}
  else {
    # transparent theme of ggplot2
    transparent_theme <- ggplot2::theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # legend.position="none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )
    transparent_theme2 <- ggplot2::theme(
      #axis.line = element_blank(),
      axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # legend.position="none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )


    y <- NULL
    term.node <- NULL
    yi <- NULL
    leaf.no <- NULL
    tree <- object$tree
    tree <- tree[!is.na(tree$pleaf), ]
    # first, grow the tree from the root node with order and
    # find out the order of terminal nodes
    count <- 1
    nodes <- data.frame(leaf=1, pleaf=0, x=0, y=0, w=1)

    for(pleaf in tree$pleaf) {
      pleaf_row <- nodes[nodes$leaf == pleaf, ]

      # Add left child.
      count <- count + 1
      # reduce the width in half, in order to calculate the
      # terminal X-axis values
      nodes <- updateNodes(
        nodes,
        data.frame(
          leaf = count,
          pleaf = pleaf,
          x = pleaf_row$x - pleaf_row$w / 2,
          y = pleaf_row$y - 1,
          w = pleaf_row$w / 2
        )
      )
      # Add right child.
      count <- count + 1
      nodes <- updateNodes(
        nodes,
        data.frame(
          leaf = count,
          pleaf = pleaf,
          x = pleaf_row$x + pleaf_row$w / 2,
          y = pleaf_row$y - 1,
          w = pleaf_row$w / 2
        )
      )
    }
    # Add split conditions
    nodes$split = NA
    nodes$split[tree$pleaf] = tree$split

    # Second, find the new x.coordinates
    nodes$x.new <- rep(NA, nrow(nodes))
    # fix the x.coord for terminal nodes
    inx.term <- !(nodes$leaf %in% nodes$pleaf)
    nodes.term <- nodes[inx.term,]
    nodes$x.new[inx.term] <- rank(nodes$x[inx.term])
    nodes$leaf.no[inx.term] <- object$n[as.character(nodes$leaf[inx.term])]
    # fix the x.coord for parent nodes
    for (i in min(nodes$y):-1){
      inx.pleaf <- which(nodes$y == i)
      coords <- sapply(split(nodes$x.new[inx.pleaf], nodes$pleaf[inx.pleaf]), mean)
      leaf.no <- sapply(split(nodes$leaf.no[inx.pleaf], nodes$pleaf[inx.pleaf]), sum)
      inx.replace <- names(coords[!is.na(coords)])
      nodes$x.new[as.numeric(inx.replace)] <- coords[inx.replace]
      nodes$leaf.no[as.numeric(inx.replace)] <- leaf.no[inx.replace]
    }
    # replace the x.coord
    nodes$x <- nodes$x.new


    config.leaf_width_scale <- 0.9

    # find the good size for the ovals representing nodes
    # half of the minimun distance between adjacent node centroids
    x.scale <- config.leaf_width_scale / 2 *
      min(sapply(split(nodes[-1, ]$x,f =nodes[-1, ]$y), function(x) min(diff(sort(x)))))
    y.scale <- x.scale*diff(range(nodes$y))/diff(range(nodes$x))


    # Build the plot
    vis <- ggplot()
    # Add lines first
    for(i in 1:nrow(nodes)){
      node <- nodes[i, ]
      # skip root
      if(node$pleaf == 0){
        next
      }
      parent = nodes[nodes$leaf == node$pleaf, ]
      data_line = data.frame(x = c(node$x, parent$x),
                             y = c(node$y, parent$y))
      vis <- vis + geom_line(data = data_line, aes(x, y), color = "black")
    }

    config.branch_text_left_dx = -0.2
    config.branch_text_right_dx = 0.2
    config.branch_text_left = "Yes"
    config.branch_text_right = "No"
    config.branch_text_size = 3

    config.leaf_oval_ratio = 1.3
    config.leaf_text_size = 5

    config.split_text_dy = -0.33
    config.split_text_size = 3
    config.split_label = T

    # TODO: customize the leaf text
    for (i in 1:nrow(nodes)) {
      node <- nodes[i, ]
      parent = nodes[nodes$leaf == node$pleaf,]
      # Add nodes
      vis <- oval_draw(vis, node$x, node$y, config.leaf_oval_ratio, x.scale, y.scale) +
        geom_text(
          data = data.frame(x = node$x, y = node$y),
          aes(x, y),
          label = paste("K =",node$leaf.no),
          size = config.leaf_text_size
        )

      # The height difference between two levels is used as a baseline
      # TODO: customizable
      h = 1

      # Add split conditions in case it's a splitting node
      if(!is.na(node$split)){
        dy <- h * config.split_text_dy
        data_text = data.frame(x = node$x, y = node$y + dy)
        show_text = ifelse(config.split_label, geom_label, geom_text)
        vis <- vis +
          show_text(
            data = data_text,
            aes(x, y),
            label = encodeHtml(node$split),
            size = config.split_text_size
          )
      }


      # Add `yes/no`` text in branch lines
      # Calculate offset of text, to avoid overlapping lines of branch.
      dx = h * ifelse(node$leaf %% 2 == 0,
                      config.branch_text_left_dx,
                      config.branch_text_right_dx)
      data_text = data.frame(x = (node$x + parent$x) / 2 + dx,
                             y = (node$y + parent$y) / 2)
      vis <- vis +
        geom_text(
          data = data_text,
          aes(x, y),
          label = ifelse(
            node$leaf %% 2 == 0,
            config.branch_text_left,
            config.branch_text_right
          ),
          size = config.branch_text_size
        )
    }

    vis <- vis + transparent_theme

    # plot the confident intervals
    term <- nodes[is.na(nodes$split),]
    # with ascending x
    term <- term[ordered(term$x.new),]
    yi <- model.response(x$data)
    p <- ggplot()
    p <- p + geom_hline(
      yintercept = c(min(yi), max(yi)),
      linetype = "solid"
    )
    p <- p + geom_hline(
      yintercept = 0,
      linetype = "dashed"
    )
    p <- p + scale_x_discrete(limits = as.factor(term$leaf))
    # # can be adjusted in future
    CI.ratio = 2
    for (i in unique(term.nodes)) {
      y.coord2 = frame0[i, ]$yval
      x.coord2 = nodes[i, ]$x
      # node names need to be changed for x$se
      p <- CI_draw(p, x = x.coord2, y = y.coord2, b = 1.96* x$se[trans.labels[i]], a = 1.96* x$se[trans.labels[i]]/CI.ratio)

    }
    p <- p + transparent_theme2
    # Finally put two plots together
    # X:    ggplotGrob + annotation_custom
    # XXX:  gtable_matrix. require fixed size.
    # X:    viewport
    # O:    gridExtra

    grid.arrange(vis, p, nrow = 2, as.table=T, heights = c(3,1))


  }
  }

}

