#!/usr/bin/mawk -f

NR==FNR{
	if($0~/>/){
		h = substr($0,2)
		next
	}
	lseq[h]+=length($0)
	next
}

/>/{
	htrue = 0
	split($0,x," ")
	parent = substr(x[1],2)
	sub(/_[0-9]+$/,"",parent)
	orfpos1 = 0 + substr(x[2],2)
	orfpos2 = 3 + x[4]
	if(orfpos2>lseq[parent])
		next
	head = $0
	htrue = 1
	next
}

htrue==1{
	orfs[head] = orfs[head] $0
}

END{
	for(i in orfs)
		print i "\n" orfs[i]
}

