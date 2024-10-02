#' Make a surface plot
#'
#' @param tab a data-frame containing spatial co-ordinates and the variable to
#' plot
#' @param coords_name name of the two columns that contains the co-ordinates of
#' the points
#' @param var_name name of the column containing the variable to be plotted
#' @param h integer; (optional) controls smoothness of the spatial interpolation
#' as appearing in the [MBA::mba.surf()] function. Default is 8.
#' @param col.pal Optional; color palette, preferably divergent, use
#' \code{colorRampPalette} function from \code{grDevices}. Default is 'RdYlBu'.
#' @param mark_points Logical; if \code{TRUE}, the input points are marked.
#' Default is \code{FALSE}.
#' @return a `ggplot` object containing the surface plot
#' @importFrom MBA mba.surf
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_distiller
#' geom_point scale_fill_gradientn
#' @importFrom ggplot2 theme_bw theme element_line element_blank element_text
#' @importFrom stats na.omit
#' @examples
#' data(simGaussian)
#' plot1 <- surfaceplot(simGaussian, coords_name = c("s1", "s2"),
#'                      var_name = "z_true")
#' plot1
#'
#' # try your favourite color palette
#' col.br <- colorRampPalette(c("blue", "white", "red"))
#' col.br.pal <- col.br(100)
#' plot2 <- surfaceplot(simGaussian, coords_name = c("s1", "s2"),
#'                      var_name = "z_true", col.pal = col.br.pal)
#' plot2
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
surfaceplot <- function(tab, coords_name, var_name, h = 8,
                        col.pal, mark_points = FALSE){

  surf <- mba.surf(tab[,c(coords_name, var_name)],
                   no.X = 250, no.Y = 250, h = h, m = 1, n = 1,
                   extend=FALSE)$xyz.est

  surf_df <- data.frame(expand.grid(surf$x, surf$y), z = as.vector(surf$z))
  surf_df <- na.omit(surf_df)
  names(surf_df) <- c("x", "y", "z")

  plot <- ggplot(surf_df, aes_string(x = 'x', y = 'y')) +
    geom_raster(aes_string(fill = 'z')) +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = 0.25),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.box.just = "center",
          aspect.ratio = 1)

  if(missing(col.pal)){
    plot <- plot + scale_fill_distiller(palette = "RdYlBu", direction = -1,
                                        label = function(x) sprintf("%.1f", x))
  }else{
    plot <- plot + scale_fill_gradientn(colours = col.pal)
  }

  if(mark_points){
    plot <- plot + geom_point(aes_string(x = coords_name[1],
                                         y = coords_name[2]),
                              data = tab, color = "black", fill = NA,
                              shape = 21, stroke = 0.5, alpha = 0.5)
  }

  plot

}

#' Make two side-by-side surface plots
#'
#' @description Make two side-by-side surface plots, particularly useful towards
#' a comparative study of two spatial surfaces.
#' @param tab a data-frame containing spatial co-ordinates and the variables to
#' plot
#' @param coords_name name of the two columns that contains the co-ordinates of
#' the points
#' @param var1_name name of the column containing the first variable to be
#' plotted
#' @param var2_name name of the column containing the second variable to be
#' plotted
#' @param h integer; (optional) controls smoothness of the spatial interpolation
#' as appearing in the [MBA::mba.surf()] function. Default is 8.
#' @param col.pal Optional; color palette, preferably divergent, use
#' \code{colorRampPalette} function from \code{grDevices}. Default is 'RdYlBu'.
#' @param mark_points Logical; if \code{TRUE}, the input points are marked.
#' Default is \code{FALSE}.
#' @return a list containing two `ggplot` objects
#' @importFrom MBA mba.surf
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_distiller
#' geom_point scale_fill_gradientn
#' @importFrom ggplot2 theme_bw theme element_line element_blank element_text
#' @importFrom stats na.omit
#' @examples
#' data(simGaussian)
#' plots_2 <- surfaceplot2(simGaussian, coords_name = c("s1", "s2"),
#'                         var1_name = "z_true", var2_name = "y")
#' plots_2
#' @author Soumyakanti Pan <span18@ucla.edu>,\cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
surfaceplot2 <- function(tab, coords_name, var1_name, var2_name,
                        h = 8, col.pal, mark_points = FALSE){

  surf1 <- mba.surf(tab[,c(coords_name, var1_name)],
                   no.X = 250, no.Y = 250, h = h, m = 1, n = 1,
                   extend=FALSE)$xyz.est

  surf2 <- mba.surf(tab[,c(coords_name, var2_name)],
                    no.X = 250, no.Y = 250, h = h, m = 1, n = 1,
                    extend=FALSE)$xyz.est

  surf_df1 <- data.frame(expand.grid(surf1$x, surf1$y), z = as.vector(surf1$z))
  surf_df2 <- data.frame(expand.grid(surf2$x, surf2$y), z = as.vector(surf2$z))
  surf_df1 <- na.omit(surf_df1)
  surf_df2 <- na.omit(surf_df2)
  names(surf_df1) <- c("x", "y", "z")
  names(surf_df2) <- c("x", "y", "z")

  plot1 <- ggplot(surf_df1, aes_string(x = 'x', y = 'y')) +
    geom_raster(aes_string(fill = 'z')) +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = 0.25),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.box.just = "center",
          aspect.ratio = 1)

  plot2 <- ggplot(surf_df2, aes_string(x = 'x', y = 'y')) +
    geom_raster(aes_string(fill = 'z')) +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = 0.25),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.box.just = "center",
          aspect.ratio = 1)

  common_limits <- range(c(surf_df1$z, surf_df2$z))

  if(missing(col.pal)){
    plot1 <- plot1 + scale_fill_distiller(palette = "RdYlBu", direction = -1,
                                          label = function(x) sprintf("%.1f", x),
                                          limits = common_limits)

    plot2 <- plot2 + scale_fill_distiller(palette = "RdYlBu", direction = -1,
                                          label = function(x) sprintf("%.1f", x),
                                          limits = common_limits)

  }else{
    plot1 <- plot1 + scale_fill_gradientn(colours = col.pal,
                                          limits = common_limits)

    plot2 <- plot2 + scale_fill_gradientn(colours = col.pal,
                                          limits = common_limits)
  }

  if(mark_points){
    plot1 <- plot1 + geom_point(aes_string(x = coords_name[1],
                                           y = coords_name[2]),
                                data = tab, color = "black", fill = NA,
                                shape = 21, stroke = 0.5, alpha = 0.5)

    plot2 <- plot2 + geom_point(aes_string(x = coords_name[1],
                                           y = coords_name[2]),
                                data = tab, color = "black", fill = NA,
                                shape = 21, stroke = 0.5, alpha = 0.5)
  }

  return(list(plot1, plot2))

}