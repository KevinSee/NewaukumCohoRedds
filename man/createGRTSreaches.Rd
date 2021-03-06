% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createGRTSreaches.R
\name{createGRTSreaches}
\alias{createGRTSreaches}
\title{Create GRTS reaches}
\usage{
createGRTSreaches(grts_pts = NULL, strm_sf = NULL, length_buffer_mi = 0.5)
}
\arguments{
\item{grts_pts}{\code{sf} object of GRTS points. Could be the output of \code{spsurvey::grts}.}

\item{strm_sf}{\code{sf} object containing the line layer of the GRTS frame. Should be in a projection with feet as the distant unit, and contain a column \code{llid} with unique codes for each stream.}

\item{length_buffer_mi}{distance (in miles) that the reach should extend on either side of the GRTS point.}
}
\value{
sf
}
\description{
Based on GRTS points and a stream layer, determines the boundaries of a reach that will be sampled as part of the study design. Strives to keep that reach on the same stream, and extend for a specified distance upstream and downstream of the GRTS point.
}
\examples{
createGRTSreaches()
}
\author{
Kevin See
}
