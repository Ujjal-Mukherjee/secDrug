PropTest <-
function(p1,p2=0,n){

Stdev=sqrt((p1*(1-p1)/n)+(p2*(1-p2)/n))
Stat=(p1-p2)/Stdev

p=2*(1-pnorm(abs(Stat),Stdev))

return(p);

}
