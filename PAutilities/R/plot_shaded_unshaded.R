#' @rdname plot.paired_equivalence
shaded_equivalence_plot <- function(results, ...) {

  results$x_label <- paste0(
    results$variable,
    ifelse(results$CI_sig == "*", "*", "")
  )

  ggplot(results, aes(x_label)) +

    geom_rect(
      ymin = results[1,"region_low"],
      ymax = results[1,"region_high"],
      xmin = -Inf,
      xmax = Inf,
      fill = "gray",
      alpha = 0.2
    ) +

    geom_errorbar(
      aes(
        ymin = CI_low, ymax = CI_high,
        linetype = factor(
          CI_sig, levels = c("*", "NS")
        )
      ),
      width = 0.25,
      size = 1.1,
      show.legend = FALSE
    ) +

    scale_linetype_manual(
      breaks = c("*", "NS"), guide = FALSE,
      drop = FALSE, values = c("solid", "dotted")
    ) +

    coord_flip() +
    theme_classic() +
    theme(
      ...,
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +

    expand_limits(
      y = c(
        min(0, results$region_low*1.1),
        results$region_high*1.1
      )
    )

}

#' @rdname plot.paired_equivalence
unshaded_equivalence_plot <- function(results, ...) {

  names(results) <- gsub(
    "^variable$", "varname", names(results)
  )
  results$x_label <- paste0(
    results$varname,
    ifelse(results$CI_sig == "*", "*", "")
  )

  results <- reshape2::melt(
    results, measure.vars = c(
      "region_low", "region_high",
      "CI_low", "CI_high"
    )
  )

  results$variable <- as.character(results$variable)
  results$vartype <- ifelse(
    grepl("^region_", results$variable),
    results$y_type,
    "prediction"
  )
  results$variable <- gsub("^.*_", "", results$variable)

  results <- reshape2::dcast(results, ...~variable)

  ggplot(results, aes(x_label)) +
    geom_errorbar(
      aes(
        ymin = low, ymax = high,
        group = vartype, linetype = vartype
      ),
      width = 0.25,
      size = 1.1,
      position = position_dodge(0.5)
    ) +

    coord_flip() +
    theme_classic() +
    theme(
      ...,
      axis.line = element_line(size = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )

}
