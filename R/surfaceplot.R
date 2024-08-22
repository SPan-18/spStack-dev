#' Make a surface plot
#'
#' @param tab a data-frame containing spatial co-ordinates and the variable to
#' plot
#' @param coords_name name of the two columns that contains the co-ordinates of
#' the points
#' @param var_name name of the column containing the variable to be plotted
#' @param h integer; controls smoothness of the spatial interpolation as
#' appearing in the \code{mba.surf} function of the \code{MBA} package.
#' Default is 8.
#' @param col.pal Optional; color palette, preferably divergent, use
#' \code{colorRampPalette} function from \code{grDevices}. Default is 'RdYlBu'.
#' @param mark_points Logical; if \code{TRUE}, the input points are marked.
#' Default is \code{FALSE}.
#' @importFrom MBA mba.surf
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_distiller
#' geom_point scale_fill_gradientn
#' @importFrom ggplot2 theme_bw theme element_line element_blank element_text
#' @importFrom stats na.omit
#' @examples
#' data(simLMdat)
#' plot1 <- surfaceplot(simLMdat, coords_name = c("s1", "s2"),
#'                      var_name = "z_true")
#' plot1
#'
#' # try your favourite color palette
#' col.br <- colorRampPalette(c("blue", "white", "red"))
#' col.br.pal <- col.br(100)
#' plot2 <- surfaceplot(simLMdat, coords_name = c("s1", "s2"),
#'                      var_name = "z_true", col.pal = col.br.pal)
#' plot2
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
