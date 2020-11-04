branded_colors <- list(
  "blue"   = "#00798c",
  "red"    = "#d1495b",
  "yellow" = "#edae49",
  "green"  = "#66a182",
  "navy"   = "#2e4057", 
  "grey"   = "#8d96a3"
)

branded_pal <- function(
  primary = "blue", 
  other = "grey", 
  direction = 1
) {
  stopifnot(primary %in% names(branded_colors))
  
  function(n) {
    if (n > 6) warning("Branded Color Palette only has 6 colors.")
    
    if (n == 2) {
      other <- if (!other %in% names(branded_colors)) {
        other
      } else {
        branded_colors[other]
      }
      color_list <- c(other, branded_colors[primary])
    } else {
      color_list <- branded_colors[1:n]
    }
    
    color_list <- unname(unlist(color_list))
    if (direction >= 0) color_list else rev(color_list)
  }
}

scale_colour_branded <- function(
  primary = "blue", 
  other = "grey", 
  direction = 1, 
  ...
) {
  ggplot2::discrete_scale(
    "colour", "branded", 
    branded_pal(primary, other, direction), 
    ...
  )
}

scale_color_branded <- scale_colour_branded