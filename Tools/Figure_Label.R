 # R add figure label
 line2user <- function(line, side) {
   lh <- par('cin')[2] * par('cex') * par('lheight')
   x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
   y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
   switch(side,
          `1` = grconvertY(-line * y_off, 'npc', 'user'),
          `2` = grconvertX(-line * x_off, 'npc', 'user'),
          `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
          `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
          stop("Side must be 1, 2, 3, or 4", call.=FALSE))
 }

 addfiglab <- function(lab, xl = par()$mar[2], yl = par()$mar[3]) {

   text(x = line2user(xl, 2), y = line2user(yl, 3),
        lab, xpd = NA, font = 2, cex = 1.5, adj = c(0, 1))

 }

 par(mfrow = c(1, 2))
 plot(0)
 addfiglab("A")
 plot(1000)
 addfiglab("B")
