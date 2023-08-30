#!/usr/bin/gawk -f

BEGIN{
	FS = OFS ="\t"
}

NR==FNR{
	if($3 == "CDS"){
		parent = $9
		sub(/;.*/,"",parent)
		sub(/.*=/,"",parent)
		b[$1][parent] = $4
		e[$1][parent] = $5
		o[$1][parent] = $7
		f[$1][parent] = $6
	}
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

$3=="transcript"{
	for(i in b[$1]){
		if($7!=o[$1][i] && (($4>=b[$1][i] && $4<=e[$1][i]) || ($5>=b[$1][i] && $5<=e[$1][i])))
			print $1, $2, $3, max($4,b[$1][i]), min($5,e[$1][i]), $6, $7, $8, $9
	}
}
