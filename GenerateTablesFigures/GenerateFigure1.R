# Black-and-white schematic. Uses only grid (included with R).
# Run: source("polished_schematic.R")
# Outputs: polished_schematic.pdf and polished_schematic.png
library(grid)

draw_schematic <- function(text_scale = 1.0) {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, 12), yscale = c(0, 8.5)))
  on.exit(popViewport())
  used <- "#BFBFBF"
  line <- gpar(col = "black", lwd = 1.2, fill = NA)
  txt <- function(label, x, y, size = 16, face = "plain",
                  just = "centre", family = "sans") {
    grid.text(
      label,
      x = unit(x, "native"),
      y = unit(y, "native"),
      just = just,
      gp = gpar(
        fontfamily = family,
        fontsize = size * text_scale,
        fontface = face,
        col = "black"
      )
    )
  }
  box <- function(x, y, w, h, dashed = TRUE) {
    grid.rect(unit(x, "native"), unit(y, "native"),
              unit(w, "native"), unit(h, "native"),
              gp = gpar(col = "black", fill = "white", lwd = 1.2,
                        lty = if (dashed) "dashed" else "solid"))
  }
  # Cubic Bezier paths sampled into a polyline for portable grid rendering.
  curve_arrow <- function(p) {
    t <- seq(0, 1, length.out = 150)
    b <- cbind((1-t)^3, 3*(1-t)^2*t, 3*(1-t)*t^2, t^3) %*% p
    grid.lines(unit(b[, 1], "native"), unit(b[, 2], "native"),
               arrow = arrow(length = unit(0.13, "inches"), type = "open"),
               gp = gpar(col = "black", lwd = 1.5, lineend = "round"))
  }
  node <- function(x, y, radius, label, is_used) {
    # snpc units preserve circular nodes when the device is resized.
    grid.circle(unit(x, "native"), unit(y, "native"),
                r = unit(radius / 8.5, "snpc"),
                gp = gpar(col = "black", fill = if (is_used) used else "white",
                          lwd = 1.2))
    txt(label, x, y + radius + 0.28, 20)
  }
  
  # Compact legend, with the original used / not-used encoding.
  for (i in 1:2) {
    yy <- c(7.96, 7.57)[i]
    grid.rect(unit(0.64, "native"), unit(yy, "native"),
              unit(0.16, "inches"), unit(0.16, "inches"),
              gp = gpar(col = "black", fill = c(used, "white")[i], lwd = 1))
    txt(c("Used", "Not used")[i], 0.90, yy, 13, just = "left")
  }
  
  # Matched panels, without the redundant outer frame.
  box(2.90, 3.15, 4.70, 4.90)
  box(9.10, 3.15, 4.70, 4.90)
  txt("External Data", 2.90, 5.88, 21)
  txt("Internal Data", 9.10, 5.88, 21)
  
  # Model learned externally; transport assumptions shown to its right.
  box(4.52, 7.45, 4.05, 0.78)
  txt(expression(rank(f[E](bold(Z))) ==
                   rank(E[E]*"{"*tilde(Y)*"|"*bold(Z)*"}")),
      4.52, 7.45, 15)
  txt(expression(E[E]*"{"*tilde(Y)*"|"*bold(Z)*"}" !=
                   E[plain(I)]*"{"*Y*"|"*bold(Z)*"}"), 9.03, 8.00, 14, family="serif")
  txt(expression(E[E]*"{"*Y*"|"*bold(Z)*"}" !=
                   E[plain(I)]*"{"*Y*"|"*bold(Z)*"}"), 9.03, 7.57, 14, family="serif")
  txt(
    expression(
      plain(cor) * bgroup(
        "[",
        rank(E[E] * "{" * tilde(Y) * "|" * bold(Z) * "}") * "," ~
          rank(E[I] * "{" * Y * "|" * bold(Z) * "}"),
        "]"
      ) > 0
    ),
    9.03, 7.10, 12, family="serif"
  )
  curve_arrow(rbind(c(1.18, 5.60), c(1.28, 6.70),
                    c(1.60, 7.45), c(2.43, 7.45)))
  curve_arrow(rbind(c(6.60, 7.30), c(7.12, 6.42),
                    c(10.80, 6.60), c(10.90, 5.65)))
  
  for (left in c(TRUE, FALSE)) {
    xs <- if (left) c(1.70, 4.10) else c(7.90, 10.30)
    node(xs[1], 4.50, 0.44, expression(tilde(Y)), left)
    node(xs[2], 4.50, 0.44, expression(Y), !left)
    node(xs[1], 1.85, 0.78, expression(bold(Z)), TRUE)
    node(xs[2], 1.85, 0.78, expression(bold(B)), !left)
  }
  invisible(NULL)
}

## Save in a 5.5 x 4.0 aspect ratio
draw_schematic(text_scale=0.67)
