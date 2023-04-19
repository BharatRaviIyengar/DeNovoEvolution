#!/usr/bin/gawk -f

# NC_001133.9	RefSeq	CDS	1807	2169
BEGIN{
	FS = OFS ="\t"
}

NR==FNR{
	if($3 == "CDS"){
		parent = $9
		sub(/.*=rna-/,"",parent)
		sub(/;.*/,"",parent)
		b[$1][parent] = $4
		e[$1][parent] = $5
		o[$1][parent] = $7
		f[$1][parent] = $8
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

{	
	for(i in b[$1]){
		if( (($4>=b[$1][i] && $4<=e[$1][i]) || ($5>=b[$1][i] && $5<=e[$1][i])) && $7 != o[$1][i]){
			if($7 == "+")
				fx = mod( ($4 - b[$1][i] + f[$1][i]),3)
			else
				fx = mod((e[$1][i]-f[$1][i]-$5),3)
			len = $5-$4+1
			$8 = fx
			$6 = sprintf("%0.1f",100*(min(e[$1][i],$5) - max(b[$1][i],$4)+1)/len)
			$9 = $9 "; " ss "anti-" i
			print $0
			next
		}
	}
}
