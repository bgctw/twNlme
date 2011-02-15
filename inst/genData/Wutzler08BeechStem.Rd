     \name{Wutzler08BeechStem}
     \docType{data}
     \alias{Wutzler08BeechStem}
     \title{Wutzler08BeechStem}
     \description{Tree weights by tree and stand properties.}
     \usage{Wutzler08BeechStem}
     \format{data.frame with columns \describe{
		\item{author}{Original author of the data}
		\item{stand}{Identifier of the forest plot}
		\item{alt}{Altitude above sea level in m}
		\item{si}{site index of the plot, i.e. expected tree height at age 100, in m }
		\item{age}{tree age in yr}
		\item{stockden}{stocking density}
		\item{dbh}{tree diameter at breastheight (1.3m) in cm}
		\item{height}{tree height in m}
		\item{stem}{stem weight in kg}
		}
     }
     \source{http://pubs.nrc-cnrc.gc.ca/rp/rpsuppl/x07-194e.RData}
     \references{
     	Wutzler, T.; Wirth, C. & Schumacher, J. (2008) Generic biomass functions for Common beech (Fagus sylvatica L.) in Central Europe - predictions and components of uncertainty. Canadian Journal of Forest Research, 38, 1661-1675
     }
     \keyword{datasets}
     \examples{
     	data(Wutzler08BeechStem)
     	plot( log(stem)~dbh, col=author, data=Wutzler08BeechStem)
     }
      
