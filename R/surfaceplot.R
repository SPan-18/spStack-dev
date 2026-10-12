#' Make a surface plot
#'
#' @param tab a data-frame containing spatial co-ordinates and the variable to
#' plot
#' @param coords_name name of the two columns that contains the co-ordinates of
#' the points, also used as the axis labels
#' @param var_name name of the column containing the variable to be plotted
#' @param h integer; (optional) controls smoothness of the spatial interpolation
#' as appearing in the [MBA::mba.surf()] function. Default is 8.
#' @param col.pal Optional; color palette, preferably divergent, use
#' \code{colorRampPalette} function from \code{grDevices}. Default is the
#' colorblind-friendly diverging palette 'RdBu' from ColorBrewer.
#' @param mark_points Logical; if \code{TRUE}, the input points are marked.
#' Default is \code{FALSE}.
#' @return a `ggplot` object containing the surface plot
#' @importFrom MBA mba.surf
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_distiller
#'   geom_point scale_fill_gradientn
#' @importFrom ggplot2 theme_bw theme element_line element_blank element_text labs
#' @importFrom stats na.omit
#' @importFrom rlang .data
#' @examples
#' data(simSpatial)
#' plot1 <- surfaceplot(simSpatial, coords_name = c("s1", "s2"),
#'                      var_name = "z_true")
#' plot1
#'
#' # try your favourite color palette
#' plot2 <- surfaceplot(simSpatial, coords_name = c("s1", "s2"),
#'                      var_name = "z_true",
#'                      col.pal = hcl.colors(100, "PuOr", rev = TRUE))
#' plot2
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
surfaceplot <- function(tab, coords_name, var_name, h = 8,
                        col.pal, mark_points = FALSE){

  surf_df <- surface_grid(tab, coords_name, var_name, h)
  plot <- surface_ggplot(surf_df, tab, coords_name, range(surf_df$z), col.pal,
                         mark_points)
  plot

}

#' Make two side-by-side surface plots
#'
#' @description Make two side-by-side surface plots, particularly useful towards
#' a comparative study of two spatial surfaces.
#' @param tab a data-frame containing spatial co-ordinates and the variables to
#' plot
#' @param coords_name name of the two columns that contains the co-ordinates of
#' the points, also used as the axis labels
#' @param var1_name name of the column containing the first variable to be
#' plotted
#' @param var2_name name of the column containing the second variable to be
#' plotted
#' @param h integer; (optional) controls smoothness of the spatial interpolation
#' as appearing in the [MBA::mba.surf()] function. Default is 8.
#' @param col.pal Optional; color palette, preferably divergent, use
#' \code{colorRampPalette} function from \code{grDevices}. Default is the
#' colorblind-friendly diverging palette 'RdBu' from ColorBrewer.
#' @param mark_points Logical; if \code{TRUE}, the input points are marked.
#' Default is \code{FALSE}.
#' @return a list containing two `ggplot` objects
#' @importFrom MBA mba.surf
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_distiller
#'   geom_point scale_fill_gradientn
#' @importFrom ggplot2 theme_bw theme element_line element_blank element_text labs
#' @importFrom stats na.omit
#' @importFrom rlang .data
#' @examples
#' data(simSpatial)
#' plots_2 <- surfaceplot2(simSpatial, coords_name = c("s1", "s2"),
#'                         var1_name = "z_true", var2_name = "y_gauss")
#' plots_2
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
surfaceplot2 <- function(tab, coords_name, var1_name, var2_name,
                        h = 8, col.pal, mark_points = FALSE){

  surf_df1 <- surface_grid(tab, coords_name, var1_name, h)
  surf_df2 <- surface_grid(tab, coords_name, var2_name, h)

  # one colour scale for both surfaces, so that a colour means the same value
  common_limits <- range(c(surf_df1$z, surf_df2$z))

  plot1 <- surface_ggplot(surf_df1, tab, coords_name, common_limits, col.pal,
                          mark_points)
  plot2 <- surface_ggplot(surf_df2, tab, coords_name, common_limits, col.pal,
                          mark_points)

  return(list(plot1, plot2))

}

# internal: surface of var_name interpolated by MBA on a 250 x 250 grid over
# the extent of the coordinates, as a data frame with columns x, y, z
surface_grid <- function(tab, coords_name, var_name, h){

  surf <- mba.surf(tab[, c(coords_name, var_name)],
                   no.X = 250, no.Y = 250, h = h, m = 1, n = 1,
                   extend = FALSE)$xyz.est
  surf_df <- data.frame(expand.grid(surf$x, surf$y), z = as.vector(surf$z))
  surf_df <- na.omit(surf_df)
  names(surf_df) <- c("x", "y", "z")
  surf_df

}

# internal: raster plot of a surface from surface_grid(), with the fill scale on
# the given limits, axes labelled by coords_name, and optionally the observed
# points of tab marked
surface_ggplot <- function(surf_df, tab, coords_name, limits, col.pal,
                           mark_points){

  plot <- ggplot(surf_df, aes(x = .data$x, y = .data$y)) +
    geom_raster(aes(fill = .data$z)) +
    labs(x = coords_name[1], y = coords_name[2]) +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = 0.25),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.box.just = "center",
          aspect.ratio = 1)

  if(missing(col.pal)){
    plot <- plot + scale_fill_distiller(palette = "RdBu", direction = -1,
                                        labels = function(x) sprintf("%.1f", x),
                                        limits = limits)
  }else{
    plot <- plot + scale_fill_gradientn(colours = col.pal, limits = limits)
  }

  if(mark_points){
    pts <- data.frame(x = tab[[coords_name[1]]], y = tab[[coords_name[2]]])
    plot <- plot + geom_point(aes(x = .data$x, y = .data$y), data = pts,
                              color = "black", fill = NA, shape = 21,
                              stroke = 0.5, alpha = 0.5)
  }

  plot

}
