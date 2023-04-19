#!/usr/bin/gawk -f

# NC_001133.9	RefSeq	CDS	1807	2169
BEGIN{
	FS = OFS ="\t"
}

NR==FNR{
	parent = $9
	sub(/.*ID=cds-/,"",parent)
	sub(/;.*/,"",parent)
	b[$1][FNR] = $4
	e[$1][FNR] = $5
	o[$1][FNR] = $7
	p[$1][FNR] = parent
	next
}

function min(x,y){
	if(x<y)
		return x
	else
		return y
}

function max(x,y){
	if(x>y)
		return x
	else
		return y
}


function abs(x){
	if(x<0)
		return 0-x
	else
		return x
}

function mod(x,y, k){
	k = x%y
	if(k<0)
		return y+k
	else
		return k
}

{	as = 0
	for(i in b[$1]){
		if( (($4>=b[$1][i] && $4<=e[$1][i]) || ($5>=b[$1][i] && $5<=e[$1][i])) && $7 != o[$1][i]){
			as = 1
			y = b[$1][i]
			z = e[$1][i]
			if($7 == "+")
				f = mod(($4-b[$1][i]),3)
			else
				f = mod((e[$1][i]-$5),3)
			ol = sprintf("%0.1f",100*(min(e[$1][i],$5) - max(b[$1][i],$4)+1)/$6)
			ss = p[$1][i]
			break
		}
	}
	if(as == 1){
		print $0, f, ol, ss
	}
}
