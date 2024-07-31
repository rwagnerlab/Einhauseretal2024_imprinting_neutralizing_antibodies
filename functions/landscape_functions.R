library(r3js)

source("./functions/remove_reactivity_bias.R")

# set plotting parameters
plot_zlim <- c(-1, 10)

plot_xlim <-  read.csv("./data/metadata/xlim_zoom.csv")$x
plot_ylim <-  read.csv("./data/metadata/ylim_zoom.csv")$x

# Set viewing angles
angle <- list(
  rotation = c(-1.377, 0, -0.1350),# rotation = c(-1.4370, 0.0062, -0.5350),
  translation = c(0, 0.05,0.1), #translation = c(0.0344, 0.0459, 0.1175),
  #zoom = 1.4342
  zoom = 1.2646
)


# Functions to remove buttons
addObject3js <- function(
  data3js,
  object,
  number_of_ids = 1
){
  
  # Generate an object ID
  if(is.null(data3js$lastID)){ data3js$lastID <- 0 }
  object$ID <- max(data3js$lastID) + seq_len(number_of_ids)
  
  # If object is interactive and highlighted add a reference to itself to
  # it's highlight group by default
  if(!is.null(object$properties$interactive)){
    object$group <- object$ID
  }
  
  # Add the object to the plot data
  data3js$plot[[length(data3js$plot)+1]] <- object
  
  # Update the ID of the last object added
  data3js$lastID <- object$ID
  
  # Return the new data
  data3js
  
}


remove_buttons <- function(data3js){
  
  new_data3js = data3js
  
  new_data3js = data3js
  
  new_data3js[['lastID']] = 0
  new_data3js[['plot']] = list()
  
  N = data3js[['lastID']] 
  
  
  
  
  for (i in 1:N)
  {
    obj = data3js[['plot']][[i]]
    
    
    
    if ('toggle' %in% names(obj[['properties']])){
      obj[['properties']][['toggle']] <- NULL
    }
    
    new_data3js = addObject3js(new_data3js,obj)
    
    
  }
  
  
  
  return (new_data3js)
  
}





# # Function to plot individual serum group landscape for fitted landscape
plot_sr_group_lndscp <- function(map, sr_group_logtiters, sr_group_cone_coords, sr_group_colbases, sr_group_cone_slopes, sr_group, 
                                 remove_buttons = TRUE, adjust_reactivity_bias = TRUE) {
  
  # # Work out which sera are out of bounds
  sr_x_coords <- srCoords(map)[,1]
  sr_y_coords <- srCoords(map)[,2]
  margin <- 0.2
  sr_out_of_bound <- sr_x_coords < plot_xlim[1] + margin |
    sr_x_coords > plot_xlim[2] - margin |
    sr_y_coords < plot_ylim[1] + margin |
    sr_y_coords > plot_ylim[2] - margin
  
  lndscp_xlim <- range(agCoords(map)[,1], na.rm = T)
  lndscp_ylim <- range(agCoords(map)[,2], na.rm = T)
  
  if(adjust_reactivity_bias) {
    sr_group_logtiters_adj <- remove_reactivity_bias_logtiter(sr_group_logtiters)
    
    sr_group_mean_logtiters <- rowMeans(sr_group_logtiters_adj, na.rm = T)
    sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
    
    
  } else {
    sr_group_mean_logtiters <- rowMeans(sr_group_logtiters, na.rm = T)
    sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
    
  }
  
  
  
  # Set map subset
  map_subset <- map
  srShown(map_subset)[sr_out_of_bound] <- FALSE
  srShown(map_subset)[srGroups(map_subset) != sr_group] <- FALSE
  
  srShown(map_subset) <- FALSE
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 3,
      cex.basemap.sr = 1.5,
      # lwd.grid             = 1.5,
      opacity.basemap.sr = 1
    )
  )
  
  # Get fitted surface
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  fit_lndscp_val <- function(x, y, sr_cone_coords, sr_colbases, sr_cone_slopes) {
    
    sr_distances <- as.matrix(dist(rbind(c(x, y), sr_cone_coords)))[1, -1]
    mean(sr_colbases - sr_distances*sr_cone_slopes)
    
  }
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        grid_x_matrix[n], grid_y_matrix[n],
        sr_group_cone_coords,
        sr_group_colbases,
        sr_group_cone_slopes
      )
    }, numeric(1)
  )
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = FALSE,
    doubleSide = TRUE
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = adjustcolor(
      sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    ),
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = TRUE
  )
  
  # Add individual landscapes
  for (i in seq_along(sr_group_colbases)) {
    
    ## Calculate the individual surface
    sr_cone_coord <- sr_group_cone_coords[i,]
    sr_colbase <- sr_group_colbases[i]
    sr_cone_slope <- sr_group_cone_slopes#[i]
    
    grid_dists <- as.matrix(dist(rbind(
      sr_cone_coord,
      cbind(
        as.vector(grid_x_matrix),
        as.vector(grid_y_matrix)
      )
    )))[1, -1]
    
    grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
    grid_z_matrix[] <- sr_colbase - grid_dists*sr_cone_slope
    
    # Add the individual surface
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = "grey70",
      opacity = 0.2,
      toggle = "Individual landscapes",
      wireframe = FALSE,
      doubleSide = TRUE
    )
    
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = "grey70",
      opacity = 0.4,
      toggle = "Individual landscapes",
      wireframe = TRUE
    )
    
  }
  #
  # # Add the titers
  data3js <- lndscp3d_titers(
    data3js = data3js,
    object = list(
      coords = agCoords(map)[!is.na(sr_group_mean_logtiters),],
      logtiters = sr_group_mean_logtiters[!is.na(sr_group_mean_logtiters)],
      indices = which(!is.na(sr_group_mean_logtiters)),
      acmap = map
    ),
    zlim = plot_zlim,
    options = list(
      cex.titer = 1,
      col.impulse = "grey60"
    )
  )
  
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(log2(5), length(x_grid), length(y_grid))
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey80",
    opacity = 0.2,
    toggle = "Titer 50"
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey40",
    opacity = 0.4,
    toggle = "Titer 50",
    wireframe = TRUE
  )
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  
  if(remove_buttons) {
    data3js <- remove_buttons(data3js)
    
  }
  
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom
  )
  
  htmlwidgets::onRender(
    widget,
    jsCode = paste0("function(el, x, data){
    el.style.outline = 'solid 2px #eeeeee';
    }")
  )
  
}

# Function to plot serum group landscape for fitted parameters
plot_sr_group_lndscp_gmt <- function(map, landscape_pars, sr_groups, remove_buttons =TRUE, adjust_reactivity_bias = TRUE) {
  
  # # Work out which sera are out of bounds
  sr_x_coords <- srCoords(map)[,1]
  sr_y_coords <- srCoords(map)[,2]
  margin <- 0.2
  sr_out_of_bound <- sr_x_coords < plot_xlim[1] + margin |
    sr_x_coords > plot_xlim[2] - margin |
    sr_y_coords < plot_ylim[1] + margin |
    sr_y_coords > plot_ylim[2] - margin
  
  lndscp_xlim <- range(agCoords(map)[,1], na.rm = T)
  lndscp_ylim <- range(agCoords(map)[,2], na.rm = T)
  
  sr_group <- sr_groups[1]
  
  if(adjust_reactivity_bias) {
    
  
    landscape_pars[[sr_group]]$log_titers <- remove_reactivity_bias_logtiter(landscape_pars[[sr_group]]$log_titers)
    landscape_pars[[sr_group]]$sr_colbases <- apply(landscape_pars[[sr_group]]$log_titers, 2, function(x) max(x, na.rm = T))
    
  }
  
 
  # Set map subset
  map_subset <- map
  srShown(map_subset)[sr_out_of_bound] <- FALSE
  srShown(map_subset)[srGroups(map_subset) != sr_group] <- FALSE

  # do not show any sera
  srShown(map_subset) <- FALSE
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 3,
      cex.basemap.sr = 1.5,
      # lwd.grid             = 1.5,
      opacity.basemap.sr = 1
    )
  )
  
  # Get fitted surface
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  fit_lndscp_val <- function(x, y, sr_cone_coords, sr_colbases, sr_cone_slopes) {
    
    sr_distances <- as.matrix(dist(rbind(c(x, y), sr_cone_coords)))[1, -1]
    mean(sr_colbases - sr_distances*sr_cone_slopes)
    
  }
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        grid_x_matrix[n], grid_y_matrix[n],
        landscape_pars[[sr_group]]$sr_cone_coords,
        landscape_pars[[sr_group]]$colbases,
        landscape_pars[[sr_group]]$slope
      )
    }, numeric(1)
  )
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
    opacity = 0.8,
    toggle = paste0(sr_group, " landscape"),
    wireframe = FALSE,
    doubleSide = TRUE
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = adjustcolor(
      sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    ),
    opacity = 0.8,
    toggle = paste0(sr_group, " landscape"),
    wireframe = TRUE
  )
  
  # Add other gmt_surfaces
  for (i in 2:length(sr_groups)) {
    
    sr_group <- sr_groups[[i]]
    
    if(adjust_reactivity_bias) {
      
      landscape_pars[[sr_group]]$log_titers <- remove_reactivity_bias_logtiter(landscape_pars[[sr_group]]$log_titers)
      landscape_pars[[sr_group]]$sr_colbases <- apply(landscape_pars[[sr_group]]$log_titers, 2, function(x) max(x, na.rm = T))
      
    }
    
    grid_z_matrix[] <- vapply(
      seq_len(length(grid_z_matrix)),
      \(n) {
        fit_lndscp_val(
          grid_x_matrix[n], grid_y_matrix[n],
          landscape_pars[[sr_group]]$sr_cone_coords,
          landscape_pars[[sr_group]]$colbases,
          landscape_pars[[sr_group]]$slope
        )
      }, numeric(1)
    )
    
    # Add the surface
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
      opacity = 0.8,
      toggle = paste0(sr_group, " landscape"),
      wireframe = FALSE,
      doubleSide = TRUE
    )
    
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = adjustcolor(
        sr_group_colors$Color[match(sr_group, sr_group_colors$Serum.group)],
        red.f = 0.25,
        green.f = 0.25,
        blue.f = 0.25
      ),
      opacity = 0.8,
      toggle = paste0(sr_group, " landscape"),
      wireframe = TRUE
    )
    
  }
  
  
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(log2(5), length(x_grid), length(y_grid))
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey80",
    opacity = 0.2,
    toggle = "Titer 50"
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey40",
    opacity = 0.4,
    toggle = "Titer 50",
    wireframe = TRUE
  )
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  if(remove_buttons) {
    data3js <- remove_buttons(data3js)
  }
  
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom
  )
  
  htmlwidgets::onRender(
    widget,
    jsCode = paste0("function(el, x, data){
    el.style.outline = 'solid 2px #eeeeee';
    }")
  )
  
}


# Function to plot serum group landscape for not fitted
plot_sr_group_lndscp_from_map <- function(map, sr_group, remove_buttons = FALSE, adjust_reactivity_bias = TRUE) {
  # Work out which sera are out of bounds
  sr_x_coords <- srCoords(map)[,1]
  sr_y_coords <- srCoords(map)[,2]
  margin <- 0.2
  sr_out_of_bound <- sr_x_coords < plot_xlim[1] + margin | 
    sr_x_coords > plot_xlim[2] - margin | 
    sr_y_coords < plot_ylim[1] + margin |
    sr_y_coords > plot_ylim[2] - margin
  
  lndscp_xlim <- range(agCoords(map)[,1], na.rm = T)
  lndscp_ylim <- range(agCoords(map)[,2], na.rm = T)
  
  # Set sera group
  sr_group_color <- unique(srOutline(map)[srGroups(map) == sr_group])
  sr_group_logtiters <- adjustedLogTiterTable(map)[ ,srGroups(map) == sr_group]
  sr_group_colbases <- colBases(map)[srGroups(map) == sr_group]
  
  if(adjust_reactivity_bias) {
    sr_group_logtiters_adj <- remove_reactivity_bias_logtiter(sr_group_logtiters)
    
    if(length(srGroups(map)[srGroups(map) == sr_group]) > 1) {
      sr_group_mean_logtiters <- rowMeans(sr_group_logtiters_adj, na.rm = T)
    } else {
      sr_group_mean_logtiters <- sr_group_logtiters_adj
    }
    
    
  } else {
    
    if(length(srGroups(map)[srGroups(map) == sr_group]) > 1) {
      sr_group_mean_logtiters <- rowMeans(sr_group_logtiters, na.rm = T)
    } else {
      sr_group_mean_logtiters <- sr_group_logtiters
    }
  }
  
  
  
   sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
  
  sr_group_coords <- srCoords(map)[srGroups(map) == sr_group,]
  
  # Set map subset
  map_subset <- map
 
  srShown(map_subset)[sr_out_of_bound] <- FALSE
  srShown(map_subset)[srGroups(map_subset) != sr_group] <- FALSE
  
  # do not show any sera for now
  srShown(map_subset) <- FALSE
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 3,
      cex.basemap.sr = 1.5,
      # lwd.grid             = 1.5,
      opacity.basemap.sr = 1
    )
  )
  
  # Get fitted surface
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  fit_lndscp_val <- function(x, y, sr_coords, sr_colbases) {
    
    sr_distances <- as.matrix(dist(rbind(c(x, y), sr_coords)))[1, -1]
    mean(sr_colbases - sr_distances, na.rm = T)
    
  }
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        grid_x_matrix[n], grid_y_matrix[n],
        sr_group_coords,
        sr_group_colbases
      )
    }, numeric(1)
  )
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = sr_group_color,
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = FALSE,
    doubleSide = TRUE
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = adjustcolor(
      sr_group_color,
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    ),
    opacity = 0.8,
    toggle = "Mean landscape",
    wireframe = TRUE
  )
  
  # Add individual landscapes
  if(length(srGroups(map)[srGroups(map) == sr_group]) > 1) {
    for (i in seq_along(sr_group_colbases)) {
      
      # Calculate the individual surface
      sr_coord <- sr_group_coords[i,]
      sr_colbase <- sr_group_colbases[i]
      
      grid_dists <- as.matrix(dist(rbind(
        sr_coord,
        cbind(
          as.vector(grid_x_matrix),
          as.vector(grid_y_matrix)
        )
      )))[1, -1]
      
      grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
      grid_z_matrix[] <- sr_colbase - grid_dists
      
      # Add the individual surface
      data3js <- r3js::surface3js(
        data3js,
        x = grid_x_matrix,
        y = grid_y_matrix,
        z = grid_z_matrix,
        col = "grey70",
        opacity = 0.2,
        toggle = "Individual landscapes",
        wireframe = FALSE,
        doubleSide = TRUE
      )
      
      data3js <- r3js::surface3js(
        data3js,
        x = grid_x_matrix,
        y = grid_y_matrix,
        z = grid_z_matrix,
        col = "grey70",
        opacity = 0.4,
        toggle = "Individual landscapes",
        wireframe = TRUE
      )
      
    }
    
  }
  
  
  # Add the titers
  data3js <- lndscp3d_titers(
    data3js = data3js,
    object = list(
      coords = agCoords(map)[!is.na(sr_group_mean_logtiters),],
      logtiters = sr_group_mean_logtiters[!is.na(sr_group_mean_logtiters)],
      indices = which(!is.na(sr_group_mean_logtiters)),
      acmap = map
    ),
    zlim = plot_zlim,
    options = list(
      cex.titer = 1,
      col.impulse = "grey60"
    )
  )
  
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(log2(5), length(x_grid), length(y_grid))
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey80",
    opacity = 0.2,
    toggle = "Titer 50"
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey40",
    opacity = 0.4,
    toggle = "Titer 50",
    wireframe = TRUE
  )
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  
  if(remove_buttons) {
    data3js <- remove_buttons(data3js)
  }
  
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom
  )
  
  htmlwidgets::onRender(
    widget,
    jsCode = paste0("function(el, x, data){
    el.style.outline = 'solid 2px #eeeeee';
    }")
  )
  
}


# Function to plot serum group landscape for not fitted
plot_gmt_lndscps_from_map <- function(map, sr_groups = unique(srGroups(map)), remove_buttons = FALSE, adjust_reactivity_bias = TRUE) {
  # Work out which sera are out of bounds
  sr_x_coords <- srCoords(map)[,1]
  sr_y_coords <- srCoords(map)[,2]
  margin <- 0.2
  sr_out_of_bound <- sr_x_coords < plot_xlim[1] + margin | 
    sr_x_coords > plot_xlim[2] - margin | 
    sr_y_coords < plot_ylim[1] + margin |
    sr_y_coords > plot_ylim[2] - margin
  
  lndscp_xlim <- range(agCoords(map)[,1], na.rm = T)
  lndscp_ylim <- range(agCoords(map)[,2], na.rm = T)
  
  sr_group <- sr_groups[1]
  # Set sera group
  sr_group_color <- unique(srOutline(map)[srGroups(map) == sr_group])
  sr_group_logtiters <- adjustedLogTiterTable(map)[ ,srGroups(map) == sr_group]
  
  sr_group_colbases <- colBases(map)[srGroups(map) == sr_group]
  
  if(adjust_reactivity_bias) {
      sr_group_logtiters <- remove_reactivity_bias_logtiter(sr_group_logtiters)
      sr_group_colbases <- apply(sr_group_logtiters, 2, function(x) max(x, na.rm = T))
  } 
  
  sr_group_mean_logtiters <- sr_group_logtiters
  sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
  
  sr_group_coords <- srCoords(map)[srGroups(map) == sr_group,]
 
  
  # Set map subset
  map_subset <- map
  srShown(map_subset)[sr_out_of_bound] <- FALSE
  srShown(map_subset)[srGroups(map_subset) != sr_group] <- FALSE
  
  # do not show any sera for now
  srShown(map_subset) <- FALSE
  # Plot the base plot
  data3js <- lndscp_r3jsplot(
    fit = list(acmap = map_subset),
    aspect.z = 0.5,
    show.surface = FALSE,
    show.titers = FALSE,
    output = "data3js",
    xlim = plot_xlim,
    ylim = plot_ylim,
    zlim = plot_zlim,
    show.sidegrid = TRUE,
    show.axis = FALSE,
    options = list(
      opacity.basemap.ags = 1,
      cex.basemap.ags = 3,
      cex.basemap.sr = 1.5,
      # lwd.grid             = 1.5,
      opacity.basemap.sr = 1
    )
  )
  
  # Get fitted surface
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  fit_lndscp_val <- function(x, y, sr_coords, sr_colbases) {
    
    sr_distances <- as.matrix(dist(rbind(c(x, y), sr_coords)))[1, -1]
    mean(sr_colbases - sr_distances, na.rm = T)
    
  }
  
  grid_z_matrix[] <- vapply(
    seq_len(length(grid_z_matrix)),
    \(n) {
      fit_lndscp_val(
        grid_x_matrix[n], grid_y_matrix[n],
        sr_group_coords,
        sr_group_colbases
      )
    }, numeric(1)
  )
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = sr_group_color,
    opacity = 0.8,
    toggle = sr_group,
    wireframe = FALSE,
    doubleSide = TRUE
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = grid_x_matrix,
    y = grid_y_matrix,
    z = grid_z_matrix,
    col = adjustcolor(
      sr_group_color,
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    ),
    opacity = 0.8,
    toggle = sr_group,
    wireframe = TRUE
  )
  
  # Add individual landscapes
  for (i in 2:length(sr_groups)) {
    
    sr_group <- sr_groups[i]
    
    # Set sera group
    sr_group_color <- unique(srOutline(map)[srGroups(map) == sr_group])
    sr_group_logtiters <- adjustedLogTiterTable(map)[ ,srGroups(map) == sr_group]
    
    sr_group_colbases <- colBases(map)[srGroups(map) == sr_group]
    
    if(adjust_reactivity_bias) {
      sr_group_logtiters <- remove_reactivity_bias_logtiter(sr_group_logtiters)
      sr_group_colbases <- apply(sr_group_logtiters, 2, function(x) max(x, na.rm = T))
    }
    
    sr_group_mean_logtiters <- sr_group_logtiters
    
    
    sr_group_mean_logtiters[sr_group_mean_logtiters < plot_zlim[1]] <- plot_zlim[1]
    
    sr_group_coords <- srCoords(map)[srGroups(map) == sr_group,]
   
    # Set map subset
    map_subset <- map
    srShown(map_subset)[sr_out_of_bound] <- FALSE
    srShown(map_subset)[srGroups(map_subset) != sr_group] <- FALSE
    
    
    # Get fitted surface
    grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
    
    grid_z_matrix[] <- vapply(
      seq_len(length(grid_z_matrix)),
      \(n) {
        fit_lndscp_val(
          grid_x_matrix[n], grid_y_matrix[n],
          sr_group_coords,
          sr_group_colbases
        )
      }, numeric(1)
    )
    
    # Add the individual surface
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = sr_group_color,
      opacity = 0.8,
      toggle = sr_group,
      wireframe = FALSE,
      doubleSide = TRUE
    )
    
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = sr_group_color,
      opacity = 0.8,
      toggle = sr_group,
      wireframe = TRUE
    )
    
  }
  
  # # Add the titers
  # data3js <- lndscp3d_titers(
  #   data3js = data3js,
  #   object = list(
  #     coords = agCoords(map)[!is.na(sr_group_mean_logtiters),],
  #     logtiters = sr_group_mean_logtiters[!is.na(sr_group_mean_logtiters)],
  #     indices = which(!is.na(sr_group_mean_logtiters)),
  #     acmap = map
  #   ),
  #   zlim = plot_zlim,
  #   options = list(
  #     cex.titer = 1.5,
  #     col.impulse = "grey60"
  #   )
  # )
  # 
  # Add the titer 50 plane
  x_grid <- seq(from = plot_xlim[1], to = plot_xlim[2], by = 0.5)
  y_grid <- seq(from = plot_ylim[1], to = plot_ylim[2], by = 0.5)
  z_grid <- matrix(log2(5), length(x_grid), length(y_grid))
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey80",
    opacity = 0.2,
    toggle = "Titer 50"
  )
  
  data3js <- r3js::surface3js(
    data3js,
    x = x_grid,
    y = y_grid,
    z = z_grid,
    col = "grey40",
    opacity = 0.4,
    toggle = "Titer 50",
    wireframe = TRUE
  )
  
  # Draw border
  data3js <- r3js::lines3js(
    data3js,
    x = c(plot_xlim[1], plot_xlim[1], plot_xlim[2], plot_xlim[2], plot_xlim[1]),
    y = c(plot_ylim[1], plot_ylim[2], plot_ylim[2], plot_ylim[1], plot_ylim[1]),
    z = rep(plot_zlim[1], 5),
    lwd = 2,
    col = "grey70"
  )
  
  if(remove_buttons) {
    data3js <- remove_buttons(data3js)
  }
  
  # Create html widget
  widget <- r3js(
    data3js = data3js,
    rotation = angle$rotation,
    translation = angle$translation,
    zoom = angle$zoom
  )
  
  htmlwidgets::onRender(
    widget,
    jsCode = paste0("function(el, x, data){
    el.style.outline = 'solid 2px #eeeeee';
    }")
  )
  
}



plot_single_landscape_panel <- function(landscape, label, label_size = 10, label_x_pos = 2, label_y_pos = 9,
                                        sr_group_label = "", sr_group_y_pos = 0, sr_group_size = 3, show_border = FALSE){
  
  
  to_save <- file.path("temp.html")
  png_save <- gsub(".html", ".png", to_save)
  saveWidget(landscape, to_save, selfcontained = FALSE)
  webshot(url=to_save,file = png_save)
  temp_plot <- readPNG(png_save)
  
  qplot(c(1:10),c(1:10), geom="blank") +
    annotation_custom(rasterGrob(temp_plot, height = unit(1, "npc")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    annotate(geom="text", x=label_x_pos, y=label_y_pos, label=label,size= label_size, hjust = 0) + 
    annotate(geom="text", x=label_x_pos, y=sr_group_y_pos, label=sr_group_label,size= sr_group_size, hjust = 0) +
    theme_void() -> p

  if(show_border) {
    p + theme(panel.border = element_rect(color = "grey50",
                                      fill = NA,
                                      size = 0.5))-> p
  }
  
  if (file.exists(to_save)) {
    #Delete file if it exists
    file.remove(to_save)
  }
  if (file.exists(png_save)) {
    #Delete file if it exists
    file.remove(png_save)
  }
  
 return(p) 
}
