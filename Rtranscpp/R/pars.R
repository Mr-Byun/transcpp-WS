# just some stuf for setting user parameters

require(compiler)
enableJIT(3)

pars <- list()

pars$default <- list()                                          
pars$default$xlog      = F                                      # plot x on log scale
pars$default$ylog      = F                                      # plot y on log scale
pars$default$adj       = 0.5                                    # string justification (0=left,1=right)
pars$default$ann       = T                                      # annotate plots with titles?
pars$default$ask       = F                                      # ask user for input
pars$default$bg        = "transparent"                          # background color
pars$default$bty       = "o"                                    # box type ("l", "7", "c", "u", or "]" or "n" for none)
pars$default$cex       = 0.83                                   # text and plotting magnification
pars$default$cex.axis  = 1                                      # axis magnification
pars$default$cex.lab   = 1                                      # label magnification
pars$default$cex.main  = 1.2                                    # main magnification
pars$default$cex.sub   = 1                                      # subtitle magnification
pars$default$col       = "black"                                # line and piont color
pars$default$col.axis  = "black"                                # axis color
pars$default$col.lab   = "black"                                # label color
pars$default$col.main  = "black"                                # main color
pars$default$col.sub   = "black"                                # subtitle color
pars$default$crt       = 0                                      # letter rotation only supported by text
pars$default$err       = 0                                      # not implemented
pars$default$family    = ""                                     # font family name (i.e "serif")
pars$default$fg        = "black"                                # foreground color (lines, axes...)
pars$default$fig       = c(0.0, 0.5, 0.0, 0.5)                  # figure display coordinates
pars$default$fin       = c(5.635622, 5.681742)                  # figure in inches
pars$default$font      = 1                                      # font (1=normal, 2=bold,3=italic,4=both,5=symbol)
pars$default$font.axis = 1                                      # axis font
pars$default$font.lab  = 1                                      # label font
pars$default$font.main = 2                                      # main font
pars$default$font.sub  = 1                                      # sub font
pars$default$lab       = c(5, 5, 7)                             # number of x and y ticks (last value unimplemented)
pars$default$las       = 0                                      # axis label orientation (0=parallel,1=horizontal,2=perpendiculat to axis,3=always vertical)
pars$default$lend      = "round"                                # style of line endings ("round","butt","square")
pars$default$lheight   = 1                                      # height of text line
pars$default$ljoin     = "round"                                # style of line joinings ("round","mitre","bevel")
pars$default$lmitre    = 10                                     # The line mitre limit. This controls when mitred line joins are automatically converted into bevelled line joins. The value must be larger than 1 and the default is 10. Not all devices will honour this setting.
pars$default$lty       = "solid"                                # The line type ("blank", "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash")
pars$default$lwd       = 1                                      # line width
pars$default$mai       = c(0.8466, 0.6806, 0.6806, 0.3486)      # margin size in inches c(bottom, left, top, right)
pars$default$mar       = c(5.1, 4.1, 4.1, 2.1)                  # number of lines of margin c(bottom, left, top, right)
pars$default$mex       = 1                                      # character size expansion 
pars$default$mfcol     = c(1, 1)                                # nrows, nfolcs
pars$default$mfg       = c(1, 1, 1, 1)                          # which to figure to draw next
pars$default$mfrow     = c(1, 1)                                # nrow, ncols
pars$default$mgp       = c(3, 1, 0)                             # The margin line (in mex units) for the axis title, axis labels and axis line
pars$default$mkh       = 0.001                                  # does nothing
pars$default$new       = F                                      # if T, do not clear frame on next plotting command 
pars$default$oma       = c(0, 0, 0, 0)                          # A vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text.
pars$default$omd       = c(0, 1, 0, 1)                          # A vector of the form c(x1, x2, y1, y2) giving the region inside outer margins in NDC (= normalized device coordinates), i.e., as a fraction (in [0, 1]) of the device region.
pars$default$omi       = c(0, 0, 0, 0)                          # A vector of the form c(bottom, left, top, right) giving the size of the outer margins in inches.
pars$default$pch       = F                                      # default plotting point
pars$default$pin       = c(5, 5)                                # current plot dimensions
pars$default$plt       = c(.1207675,.9381435,.1490036,.8802128) # A vector of the form c(x1, x2, y1, y2) giving the coordinates of the plot region as fractions of the current figure region.
pars$default$ps        = 11                                     # point size of text
pars$default$pty       = "m"                                    # "s" square or "m" max plotting region
pars$default$smo       = 1                                      # unimplemented 
pars$default$srt       = 0                                      # string rotation (only supported by text)
pars$default$tck       = NA                                     # The length of tick marks as a fraction of the smaller of the width or height of the plotting region
pars$default$tcl       = -0.5                                   # The length of tick marks as a fraction of the height of a line of text
pars$default$usr       = c(-280, 12680, -34.62356, 900.21256)   # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region.
pars$default$xaxp      = c(0, 12000, 6)                         # A vector of the form c(x1, x2, n) giving the coordinates of the extreme tick marks and the number of intervals between tick-marks
pars$default$xaxs      = "r"                                    # The style of axis interval calculation to be used for the x-axis. Possible values are "r", "i"
pars$default$xaxt      = "s"                                    # set to "n" to suppress axis
pars$default$xpd       = FALSE                                  # clip to plot (F) fig (T) or device(NA)
pars$default$yaxp      = c(0, 800, 4)                           # A vector of the form c(x1, x2, n) giving the coordinates of the extreme tick marks and the number of intervals between tick-marks
pars$default$yaxs      = "r"                                    # The style of axis interval calculation to be used for the y-axis. Possible values are "r", "i"
pars$default$yaxt      = "s"                                    # set to "n" to suppress axis
pars$default$ylbias    = 0.2                                    # A positive real value used in the positioning of text in the margins
                                                                


pars$kenneth <- list()

pars$kenneth$pch      = 19
pars$kenneth$cex      = 1.2
pars$kenneth$cex.axis = 1.1
pars$kenneth$cex.main = 1.4
pars$kenneth$cex.lab  = 1.4
pars$kenneth$cex.sub  = 1.2
pars$kenneth$lwd      = 2
pars$kenneth$tcl      = 0.4
pars$kenneth$xaxs     ='i'
pars$kenneth$yaxs     ='i'
pars$kenneth$cex      = 1.2
pars$kenneth$mar      = c(3,2,2,.5)+.5
pars$kenneth$mgp      = c(1.5,0.5,0)

pars$black <- list()
pars$black$pch      = 19
pars$black$cex      = 1.2
pars$black$cex.axis = 1.1
pars$black$cex.main = 1.4
pars$black$cex.lab  = 1.4
pars$black$cex.sub  = 1.2
pars$black$lwd      = 2
pars$black$tcl      = 0.4
pars$black$xaxs     ='i'
pars$black$yaxs     ='i'
pars$black$cex      = 1.2
pars$black$mar      = c(3,2,2,.5)+.5
pars$black$mgp      = c(1.5,0.5,0)
pars$black$bg        = "black"
pars$black$fg        = "white"
pars$black$col       = "white"
pars$black$col.axis  = "white"
pars$black$col.lab   = "white"
pars$black$col.main  = "white"
pars$black$col.sub   = "white"

