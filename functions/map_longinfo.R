#this script uses functions and code from "BA.2 and BA.5 omicron differ immunologically from both BA.1 omicron and pre-omicron variants" by Roessler, A., Netzl, A., et al. 2022
#available on repository: https://github.com/acorg/roessler_netzl_et_al2022/
#DOI: 10.5281/zenodo.7341691..
# Get names with notes including any reactivity adjustments
agNamesAdjusted <- function(map) {
  
  ag_names_adjusted <- agNames(map)
  ag_adjustments <- agReactivityAdjustments(map)
  ag_names_adjusted[ag_adjustments != 0] <- paste(
    ag_names_adjusted[ag_adjustments != 0], "∆", ag_adjustments[ag_adjustments != 0]
  )
  ag_names_adjusted
  
}

# Split apart antigen or serum groups into multiple groups based by ;
split_groups <- function(info, info_group) {

  groups_split <- strsplit(as.character(info[[info_group]]), ";")
  levels_split <- strsplit(levels(info[[info_group]]), ";")

  for (n in seq_along(groups_split[[1]])) {

    group_n  <- vapply(groups_split, function(x) trimws(x[[n]]), character(1))
    levels_n <- unique(vapply(levels_split, function(x) trimws(x[[n]]), character(1)))
    info[[paste0(info_group, n)]] <- factor(group_n, levels_n)

  }

  info

}

# Get ag labels for a map to use with ggplot
ag_labels <- function(map) {
  labels <- agNames(map)
  names(labels) <- seq_len(numAntigens(map))
  labels
}

# Get sr labels for a map to use with ggplot
sr_labels <- function(map) {
  labels <- srNames(map)
  names(labels) <- seq_len(numSera(map))
  labels
}

# Get antigen colors
ag_colors <- function(map, color_names = "ag_num", color = "fill") {
  if (color == "fill") colors <- agFill(map)
  else                 colors <- agOutline(map)
  if (color_names == "ag_num")  names(colors) <- seq_len(numAntigens(map))
  if (color_names == "ag_name") names(colors) <- agNames(map)
  colors
}

# Get sera colors
sr_colors <- function(map, color_names = "sr_num", color = "outline") {
  if (color == "fill") colors <- agFill(map)
  else                 colors <- agOutline(map)
  if (color_names == "sr_num")  names(colors) <- seq_len(numSera(map))
  if (color_names == "sr_name") names(colors) <- srNames(map)
  colors
}

# Get sera group colors
sr_group_colors <- function(map) {
  group_levels <- levels(srGroups(map))
  colors <- srOutline(map)[match(group_levels, srGroups(map))]
  colors[is.na(colors)] <- "black"
  names(colors) <- group_levels
  colors
}

# Create a facet labeller function for antigens
ag_facet_labeller <- function(map) {

  ag_nums <- seq_len(numAntigens(map))
  ag_names <- agNames(map)
  names(ag_names) <- ag_nums
  ggplot2::as_labeller(ag_names)

}

# Create a facet labeller function for sera
sr_facet_labeller <- function(map) {

  sr_nums <- seq_len(numSera(map))
  sr_names <- srNames(map)
  names(sr_names) <- sr_nums
  ggplot2::as_labeller(sr_names)

}

# Get group colors
group_cols <- function(map, groups, fill, outline) {

  # Get primary colors
  if (length(unique(fill)) > length(unique(outline))) {
    cols <- fill
  } else {
    cols <- outline
  }

  vapply(
    groups,
    function(group) {
      group_cols <- col2rgb(cols[groups == group])
      rgb(
        red   = mean(group_cols["red", ]),
        green = mean(group_cols["green", ]),
        blue  = mean(group_cols["blue", ]),
        maxColorValue = 255
      )
    },
    character(1)
  )

}

# Convert a matrix of antigen and sera data to the long format
agsr_matrix_to_long <- function(values) {

  ag_nums <- seq_len(nrow(values))
  sr_nums <- seq_len(ncol(values))

  colnames(values) <- sr_nums
  rownames(values) <- ag_nums

  values %>%
    as_tibble(
      rownames = "ag_num"
    ) %>%
    pivot_longer(
      cols = -ag_num,
      names_to = "sr_num",
      values_to = "value"
    ) %>%
    mutate(
      ag_num = factor(ag_num, levels = as.character(ag_nums)),
      sr_num = factor(sr_num, levels = as.character(sr_nums))
    )

}

# Get map antigen info in long form
#' @export
long_ag_info <- function(map) {

  ag_groups <- agGroups(map)
  if (is.null(ag_groups)) ag_groups <- factor(agNames(map))

  if (numOptimizations(map) > 0) {
    ag_reactivity_adjustments <- agReactivityAdjustments(map)
    ag_outliers <- rstatix::is_outlier(agStressPerTiter(map))
  } else {
    ag_reactivity_adjustments <- rep(0, numAntigens(map))
    ag_outliers <- rep(FALSE, numAntigens(map))
  }

  if (length(unique(agFill(map))) > length(unique(agOutline(map)))) {
    ag_cols <- agFill(map)
  } else {
    ag_cols <- agOutline(map)
  }

  ag_info <- tibble(
    ag_num = factor(seq_len(numAntigens(map))),
    ag_name = agNames(map),
    ag_name_adjusted = agNamesAdjusted(map),
    ag_group = ag_groups,
    ag_group_cols = group_cols(map, ag_groups, agFill(map), agOutline(map)),
    ag_fill = agFill(map),
    ag_outline = agOutline(map),
    ag_sequence = unlist(apply(agSequences(map), 1, list), recursive = F),
    ag_col = ag_cols,
    ag_reactivity_adjustment = ag_reactivity_adjustments,
    ag_outlier = ag_outliers
  )

  # Split groups if multiple groups are indicated
  if (grepl(";", ag_info$ag_group[1])) {
    ag_info <- split_groups(ag_info, "ag_group")
  }

  # Return info
  ag_info

}

# Get map sera info in long form
#' @export
long_sr_info <- function(map) {

  sr_groups <- srGroups(map)
  if (is.null(sr_groups)) sr_groups <- factor(srNames(map))

  if (numOptimizations(map) > 0) {
    sr_outliers <- rstatix::is_outlier(srStressPerTiter(map))
    sr_colbases <- colBases(map)
  } else {
    sr_outliers <- rep(FALSE, numSera(map))
    sr_colbases <- rep(NA_real_, numSera(map))
  }

  if (length(unique(srFill(map))) > length(unique(srOutline(map)))) {
    sr_cols <- srFill(map)
  } else {
    sr_cols <- srOutline(map)
  }

  sr_info <- tibble(
    sr_num = factor(seq_len(numSera(map))),
    sr_name = srNames(map),
    sr_group = sr_groups,
    sr_group_cols = group_cols(map, sr_groups, srFill(map), srOutline(map)),
    sr_fill = srFill(map),
    sr_outline = srOutline(map),
    sr_sequence = unlist(apply(srSequences(map), 1, list), recursive = F),
    sr_col = sr_cols,
    sr_outlier = sr_outliers,
    sr_colbase = sr_colbases
  )

  # Split groups if multiple groups are indicated
  if (grepl(";", sr_info$sr_group[1])) {
    sr_info <- split_groups(sr_info, "sr_group")
  }

  # Return info
  sr_info

}


# Get titer info in long form
#' @export
long_titer_info <- function(map) {

  # Get titer info
  titer_info <- titerTable(map)
  colnames(titer_info) <- seq_len(numSera(map))
  rownames(titer_info) <- seq_len(numAntigens(map))
  titer_info %>%
    as_tibble(
      rownames = "ag_num"
    ) %>%
    pivot_longer(
      cols = -ag_num,
      names_to = "sr_num",
      values_to = "titer"
    ) %>%
    mutate(
      logtiter = Racmacs:::log_titers(titer, 1),
      titertype = factor(as.vector(Racmacs:::titer_types_int(titer)), levels = -1:3),
      ag_num = as.factor(as.numeric(ag_num)),
      sr_num = as.factor(as.numeric(sr_num))
    )

}


# Get map info in long form
#' @export
long_map_info <- function(map) {

  # Get titer info
  titer_info <- long_titer_info(map)

  # Get antigen and sera info
  ag_info <- long_ag_info(map)
  sr_info <- long_sr_info(map)

  # Merge in information and map info
  titer_info %>%
    left_join(ag_info, by = "ag_num") %>%
    left_join(sr_info, by = "sr_num") %>%
    mutate(
      titer_adjusted = Racmacs:::reactivity_adjust_titers(titer, ag_reactivity_adjustment),
      logtiter_adjusted = logtiter + ag_reactivity_adjustment,
      dilution_stepsize = Racmacs::dilutionStepsize(map)
    )

}


#' @export
long_map_list_info <- function(maps) {

  map_info_list <- lapply(seq_along(maps), function(n) {
    long_map_info(maps[[n]]) %>%
      mutate(
        map = names(maps)[n]
      )
  })

  map_info <- do.call(bind_rows, map_info_list)
  # map_info$ag_name <- factor(map_info$ag_name, unique(map_info$ag_name))
  # map_info$sr_name <- factor(map_info$sr_name, unique(map_info$sr_name))
  select(map_info, -c("ag_num", "sr_num"))

}

