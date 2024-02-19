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

function stpf(gc,f){
	s = gc/2
	w = 0.5 -s
	if(f ==1){
		exc = 4*w^6 + 16*s*w^5 + 24*(s^2)*(w^4) + 16*(w^3)*(s^3) + 4*(s^4)*(w^2)
		return 2*s*w^2 + w^3 - exc
	} else 
		return 2*s*w^2 + w^3
}

BEGIN{
	FS = OFS = "\t"
	str = 1
	atgx = 0.01724240000000001
	stpx[0] = 0.02155216103924185
	stpx[1] = 0.03192973986845678
	stpx[2] = 0.04233394696515669
}

NR==FNR{
	gc = gsub(/[GCgc]/,"&",$2)
	l = length($2)
	gc = gc/l
	for(i=10;i<=l/3;i++){
		n_orfs = norfs(i,l)
		for(j in stpx){
			atg = (0.5*gc)*(0.5-0.5*gc)^2
			xp1[j]+= n_orfs*porfx(atg,stpf(gc,j),i)/3
			xp2[j]+= n_orfs*porfx(atgx,stpx[j],i)/3
		}
		nloci+= sprintf("%d",n_orfs/3)
	}
	next
}

{
	match($1,":")
	f = 0+substr($1,RSTART-1,1)
	ob[f]++ 
	for(i=1;i<=length($2)-30;i+=3){
		if(substr($2,i,3) == "ATG")
			ows[f]++
	} 
}

END{
	for(i=0;i<=2;i++)
		print substr(FILENAME,1,2), i, nloci, sprintf("%d",xp1[i]), sprintf("%d",xp2[i]), ob[i], ows[i]
		
}
 
