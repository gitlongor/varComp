const.polyqlist=function(coefs)
{
	do.call('polyqlist', lapply(coefs, polynomialq)) 
}
linear.polyqlist=function(const.coefs=0, linear.coefs=rep(0, length(const.coefs)), linear.name='x')
{
	L1=length(const.coefs); L2= length(linear.coefs)
	if(L1!=L2){
		L1=L2=max(L1,L2); 
		const.coefs=rep(const.coefs,length.out=L1); 
		linear.coefs=rep(linear.coefs,length.out=L2)
	}
	
	const.polyqlist(const.coefs) + const.polyqlist(linear.coefs) * rep(polynomialq(0:1), L1)
}
