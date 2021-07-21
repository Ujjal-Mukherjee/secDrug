\name{AssignGroup}
\alias{AssignGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  AssignGroup
}
\description{
%%  Assigns the Cell Lines to Resistant and Responsive Groups
%%
%% This function takes the IC50 values for each cell lines corresponding to each drug.
%% The function assigns each of the drugs into responsive and resistant using the threshold l.
%% The threshold l can be set by the user to achieve different levels of classification.
%% The main combination function automatically chooses l to maximize the out-of-sample classification accuracy.
%%

}
\usage{
AssignGroup(x, l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     Vector of IC50 Values
}
  \item{l}{
%%     Threshold for classification.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, l) 
{
    ina = which(is.na(x))
    x[ina] = mean(x[-ina])
    n = length(x)
    n1 = ceiling(l * n)
    n2 = floor((1 - l) * n)
    class = rep("N", n)
    ind = order(x)
    class[ind[1:n1]] = rep("S", n1)
    class[ind[n2:n]] = rep("R", ((n - n2) + 1))
    return(class)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
