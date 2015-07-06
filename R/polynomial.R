`*.polylist` = function(e1, e2) ## needed for polynom package
{
	if(is.numeric(e2) || is.polynomial(e2)) {
		e2 = as.polynomial(e2)
		ans = structure(lapply(e1, "*", e2), class='polylist')
	}else if(inherits(e2, 'polylist')){
		l1=length(e1); l2=length(e2)
		if(l1!=l2){
			L=max(c(l1,l2))
			e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
		}
		ans = structure(mapply("*", e1, e2, SIMPLIFY=FALSE), class='polylist')
	}else stop('unknown e2 class')
    for(i in seq_along(ans)) class(ans[[i]])=c('bigq','polynomial')
    ans
}
`/.polylist` = function(e1, e2) ## not used by varComp, just for completeness
{
	if(is.numeric(e2) || is.polynomial(e2)) {
		e2 = as.polynomial(e2)
		structure(lapply(e1, "/", e2), class='polylist')
	}else if(inherits(e2, 'polylist')){
		l1=length(e1); l2=length(e2)
		if(l1!=l2){
			L=max(c(l1,l2))
			e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
		}
		structure(mapply("/", e1, e2, SIMPLIFY=FALSE), class='polylist')
	}else stop('unknown e2 class')
}
`%%.polylist` = function(e1, e2) ## not used by varComp, just for completeness
{
	if(is.numeric(e2) || is.polynomial(e2)) {
		e2 = as.polynomial(e2)
		structure(lapply(e1, "%%", e2), class='polylist')
	}else if(inherits(e2, 'polylist')){
		l1=length(e1); l2=length(e2)
		if(l1!=l2){
			L=max(c(l1,l2))
			e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
		}
		structure(mapply("%%", e1, e2, SIMPLIFY=FALSE), class='polylist')
	}else stop('unknown e2 class')
}
`+.polylist` = function(e1, e2) ## needed for polynom package
{
	if(missing(e2)) return(e1)
	if(is.numeric(e2) || is.polynomial(e2)) {
		e2 = as.polynomial(e2)
		ans = structure(lapply(e1, "+", e2), class='polylist')
	}else if(inherits(e2, 'polylist')){
		l1=length(e1); l2=length(e2)
		if(l1!=l2){
			L=max(c(l1,l2))
			e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
		}
		ans = structure(mapply("+", e1, e2, SIMPLIFY=FALSE), class='polylist')
	}else stop('unknown e2 class')
	for(i in seq_along(ans)) class(ans[[i]])=c('bigq','polynomial')
	ans
}
`-.polylist` = function(e1, e2) ## not used by varComp, just for completeness
{
	if(missing(e2)) return(e1 * (-1))
	if(is.numeric(e2) || is.polynomial(e2)) {
		e2 = as.polynomial(e2)
		structure(lapply(e1, "-", e2), class='polylist')
	}else if(inherits(e2, 'polylist')){
		l1=length(e1); l2=length(e2)
		if(l1!=l2){
			L=max(c(l1,l2))
			e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
		}
		structure(mapply("-", e1, e2, SIMPLIFY=FALSE), class='polylist')
	}else stop('unknown e2 class')
}
`^.polylist` = function(e1, e2)
{
	stopifnot(mode(e2)=='numeric' && length(e2)==1L && e2==as.integer(e2+.5))
	structure(lapply(e1, "^", e2), class='polylist')
}

rep.polynomial=function(x, ...)
{
	structure(rep(list(x), ...), class='polylist')
}
const.polylist=function(coefs)
{
	do.call(polylist, lapply(coefs, polynomial)) 
}
linear.polylist=function(const.coefs=0, linear.coefs=rep(0, length(const.coefs)), linear.name='x')
{
	L1=length(const.coefs); L2= length(linear.coefs)
	if(L1!=L2){L1=L2=max(L1,L2); const.coefs=rep(const.coefs,length.out=L1); linear.coefs=rep(linear.coefs,length.out=L2)}
	
	const.polylist(const.coefs) + const.polylist(linear.coefs) * rep.polynomial(polynomial(0:1), L1)
}

rational=function(numer, denom=polynomial(1))rationalfun::rationalfun(numer, denom)
if(FALSE){
rational=function(num.polynomial, denom.polynomial=polynomial(1))
{
	structure(list(numerator=num.polynomial, denominator=denom.polynomial), class=c('rational','polylist'))
}
`*.rational`=function(e1, e2)
{
	if(is.numeric(e2) && length(e2)==1L) {
		e2=rational(polynomial(e2))
	}else if(is.polynomial(e2)) e2=rational(e2)
	stopifnot(inherits(e2, 'rational'))
	
	rational(e1$numerator * e2$numerator, e1$denominator * e2$denominator)
}
`+.rational`=function(e1, e2)
{
	if(missing(e2)) return(e1)
	.commonDenom(
		polylist(e1$numerator, e2$numerator), 
		polylist(e1$denominator, e2$denominator)
	)
}
`-.rational`=function(e1, e2)
{
	if(missing(e2)) return(rational(-e1$numerator, e2$denominator))
	.commonDenom(
		polylist(e1$numerator, -e2$numerator), 
		polylist(e1$denominator, e2$denominator)
	)
}
}

ratlist=function(numer.polylist, denom.polylist=rep.polynomial(polynomial(1), length.out=length(numer.polylist)))
{
	structure(mapply(rational, numer.polylist, denom.polylist, SIMPLIFY=FALSE), class='ratlist')
}
getDenominator=function(x)UseMethod('getDenominator')
getDenominator.ratlist=function(x)
{
	do.call('polylist', lapply(x, '[[', 'denominator'))
}

if(FALSE){
.commonDenom=function(num.polylist, denom.polylist) ## finds n[1]/d[1] + n[2]/d[2] + ... ### need to optimize for speed
{  
	Ln=length(num.polylist); Ld=length(denom.polylist)
	if(Ln!=Ld){Ld=Ln=max(Ld,Ln); num.polylist=rep(num.polylist, length.out=Ln); denom.polylist=rep(denom.polylist, length.out=Ld); }
	if(Ln==1) return(rational(num.polylist[[1L]], denom.polylist[[1L]]))
	denom.ans = LCM(denom.polylist)
	num.ans = polynomial(0)
	for(i in seq(Ln)){
		num.ans = num.ans + num.polylist[[i]] * (denom.ans / denom.polylist[[i]])
	}
	rational(num.ans, denom.ans)
}
}

.sum1ratlist=function(x) {
	ans0 = x[[1L]]
	if(length(x)==1L) return(ans0)
	Reduce("+", x[-1L], ans0)
}	
sum.ratlist = function(..., na.rm = FALSE)
{
	all.rat=lapply(list(...), .sum1ratlist)
	
	ans0 = all.rat[[1]]
	if(length(all.rat)==1L) return(ans0)
	Reduce("+", all.rat[-1L], ans0)	
}

if(FALSE){
`+.ratlist`=function(e1, e2)  ## not elementwise addition, but sum(e1) + sum(e2)
{
	if(missing(e2)) return(.commonDenom(e1$numerator, e1$denominator))
	if(inherits(e2, 'polynomial')){
		e2 = ratlist(polylist(e2))
	}else if(inherits(e2, 'polylist')){
		e2 = ratlist(e2)
	}else if(mode(e2)=='numeric') {
		e2 = ratlist(do.call('polylist',as.list(e2)))
	}
	
	.commonDenom(c(e1$numerator, e2$numerator), 
				c(e1$denominator, e2$denominator))
}
`-.ratlist`=function(e1, e2)
{
	if(missing(e2)) return(.commonDenom(e1$numerator * (-1), e1$denominator))
	.commonDenom(c(e1$numerator, e2$numerator * -1 ), 
				 c(e1$denominator, e2$denominator))
}
}

degree=function(x,...)UseMethod('degree')
degree.polynomial=function(x,...)
{
	seq_along(x)-1L
}
#coef.polynomial is already defined in polynom
