\name{PNAS2015}
\alias{PNAS2015}
\encoding{latin1}

\title{Small dataset from a 2015 PNAS paper}

\description{
The data was acquired by digitization of a graph from a 2015 PNAS paper. Contains three datapoints that exert significance reversal.
}

\usage{
data(PNAS2015) 
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## See examples in 'lmInfl' and 'lmThresh'.
LM <- lm(y ~ x, data = PNAS2015)
lmInfl(LM)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

