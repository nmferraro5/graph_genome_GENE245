# TEST CODE FOR R SANKEY DIAGRAM

#install.packages("googleVis")
library(googleVis)
options(gvis.plot.tag='chart')

popA <- data.frame(source = c("consensus1","consensus1","rs1-ref","rs1-alt","consensus2","consensus2","rs2-ref","rs2-alt"),
      target = c("rs1-ref","rs1-alt","consensus2","consensus2","rs2-ref","rs2-alt","consensus3","consensus3"),
      freq = c(.3,.7,1,1,0.5,0.5,1,1))
s=gvisSankey(popA[,c('source','target','freq')])
plot(s)

#Copy output of plot(s) into html file and open
#pbpaste > sankeytestX.html
#open sankeytestX.html
