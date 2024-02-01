#!/usr/bin/mawk -f

function porf(gc,k){
	s = gc/2
	w = 0.5 -s
	atg = s*w^2
	stp = 2*s*w^2 + w^3
	return atg*stp*(1-stp)^(k-2)
}

function porfx(atg,stp,k){
	return atg*stp*(1-stp)^(k-2)
}

function norfs(orflen,loc){
	return loc-3*orflen+1
}

BEGIN{
	FS = OFS = "\t"
	atgx = 0.016850003370000675
	stpx = 0.05985721197144239
}

NR==FNR{
	gc = gsub(/[GCgc]/,"&",$2)
	l = length($2)
	gc = gc/l
	for(i=10;i<=l/3;i++){
		n_orfs = 2*norfs(i,l)
		xp1+= n_orfs*porf(gc,i)
		xp2+= n_orfs*porfx(atgx,stpx,i)
		nloci+= n_orfs
	}
	next
}

{
	ob++ 
	for(i=1;i<=length($2)-30;i+=3){
		if(substr($2,i,3) == "ATG") 
			ows++
	}

}

END{

	print "IG", nloci, sprintf("%d",xp1), sprintf("%d",xp2), ob, ows
		
}
 
 
