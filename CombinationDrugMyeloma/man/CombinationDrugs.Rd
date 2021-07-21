\name{CombinationDrugs}
\alias{CombinationDrugs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CombinationDrugs(dir, file, dest = "TopCombinations_For_Bortazomib_Non_Sensitive_CellLInes.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{dest}{
%%     ~~Describe \code{dest} here~~
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
function (dir, file, dest = "TopCombinations_For_Bortazomib_Non_Sensitive_CellLInes.txt") 
{
    setwd(dir)
    d = read.csv(file, header = TRUE, as.is = TRUE)
    d1 = d[which(d$TargetCell == "B-Cell"), ]
    DrugClassCellLine = as.data.frame(matrix(nrow = dim(d1)[1], 
        ncol = dim(d1)[2]))
    rownames(DrugClassCellLine) = d1$Cell.Line
    colnames(DrugClassCellLine) = colnames(d1)
    DrugClassCellLine[, c(1:4)] = d1[, c(1:4)]
    for (k in 5:(dim(d1)[2])) {
        DrugClassCellLine[, k] = AssignGroup(d1[, k], 0.3)
    }
    frac = function(x) {
        return(length(x[which(x == "S")])/length(x))
    }
    d3 = DrugClassCellLine[-which(DrugClassCellLine$Bortezomib == 
        "S"), -c(1:4, 36)]
    IndexN = c()
    TopDrugEffect = 1
    sink(dest)
    while (TopDrugEffect >= 0.25) {
        if (length(IndexN) > 0) 
            d2 = d3[, setdiff(colnames(d3), IndexN)]
        else d2 = d3
        topCombinations = c()
        n = dim(d2)[1]
        n1 = n
        Coverage = c(n)
        pVal = 1
        while ((n1 * pVal) > 0) {
            FracS = c()
            for (i in 1:dim(d2)[2]) {
                FracS = c(FracS, frac(d2[, i]))
            }
            TopDrugs = order(FracS, decreasing = TRUE)
            pV = FracS[TopDrugs[1]]
            if (length(topCombinations) == 0) {
                IndexN = c(IndexN, colnames(d2)[TopDrugs[1]])
                TopDrugEffect = pV
            }
            if ((n1/n) <= 0.2) {
                TestP = PropTest(pV, n = dim(d2)[1])
                if (TestP <= 0.01) {
                  topCombinations = c(topCombinations, colnames(d2)[TopDrugs[1]])
                }
                else {
                  pVal = 0
                }
            }
            else {
                topCombinations = c(topCombinations, colnames(d2)[TopDrugs[1]])
            }
            ThisDrug = d2[, TopDrugs[1]]
            ind = which(ThisDrug == "S")
            d2 = d2[-ind, -TopDrugs[1]]
            n1 = dim(d2)[1]
            Coverage = c(Coverage, n1)
        }
        print("---------------------------------------")
        print(topCombinations)
        print(Coverage)
        print("---------------------------------------")
    }
    sink()
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
