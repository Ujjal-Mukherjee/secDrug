AssignGroup <-
function(x,l){

ina=which(is.na(x));
x[ina]=mean(x[-ina]);

n=length(x);
n1=ceiling(l*n);
n2=floor((1-l)*n);

class=rep("N",n);
ind=order(x);
class[ind[1:n1]]=rep("S",n1);
class[ind[n2:n]]=rep("R",((n-n2)+1));

return(class);

}
