\name{cna-deprecated}
\alias{cna-deprecated}
\alias{truthTab}
\alias{cstt}
\alias{fstt}
\alias{mvtt}
\alias{full.tt}
\alias{tt2df}


\title{Deprecated functions in the cna package}
\description{
These functions are provided for compatibility with older versions of the \pkg{cna} package only, and may be removed eventually. Commands that worked in versions of the \pkg{cna} package prior to version 3.0.0 will not necessarily work in version 3.0.0 and beyond, or may not work in the same manner.
}

\usage{
truthTab(...)
cstt(...)
fstt(...)
mvtt(...)

full.tt(...)
tt2df(...)
}

\arguments{
  \item{...}{Arguments passed to renamed functions. See the corresponding help pages for help.}
}

\details{
\code{truthTab} has been replaced by \code{\link{configTable}}.
\code{cstt} has been replaced by \code{\link{csct}}.
\code{fstt} has been replaced by \code{\link{fsct}}.
\code{mvtt} has been replaced by \code{\link{mvct}}.

\code{full.tt} has been replaced by \code{\link{full.ct}}.
\code{tt2df} has been replaced by \code{\link{ct2df}}.
}

\seealso{\code{\link{configTable}}, \code{\link{full.ct}}, \code{\link{ct2df}}}

