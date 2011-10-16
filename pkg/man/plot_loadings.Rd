\name{plot_loadings}
\alias{plot_loadings}
\title{Display a heatmap of factor loadings (note that this function requires ggplot2 package)...}
\usage{plot_loadings(x, xlabel=NA, ylabel=NA, color=NA, sorting=NA, scale.name="Loading", type="color")
plot_loadings(x, xlabel=NA, ylabel=NA, color=NA, sorting=NA, scale.name="Loading", type="color")
}
\description{Display a heatmap of factor loadings (note that this function requires ggplot2 package)}
\alias{plot_loadings}
\arguments{\item{x}{an sbfac model object or a (loadings) matrix}
\item{xlabel}{x-axis label}
\item{ylabel}{y-axis label}
\item{color}{character vector of colors (default Bl-Wh-Rd)}
\item{sorting}{a permutaion of 1:P (where P is the number of variables) providing sort order}
\item{type}{Either "color" (heatmap) or "line" (B&W alternative)
for the rows of the loadings matrix}
\item{scale.name}{label for the legend}
}
\alias{plot_loadings}
\alias{plot_loadings}
