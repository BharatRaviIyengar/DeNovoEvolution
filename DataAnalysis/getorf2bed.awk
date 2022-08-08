#!/usr/bin/mawk -f

BEGIN{
	OFS="\t"
}
/>/{
	split($0,x," ")
	match(x[1],/_[0-9]+$/)
	parent = substr(x[1],2,RSTART-2)
	num = substr(x[1],RSTART+1,RLENGTH)
	pos1 = 0 + substr(x[2],2)
	pos2 = 0 + x[4]
	f1 = 1; f2 = pos1-f1
	if(pos1<pos2){
		s = "+"
		start = pos1+0
		stop = pos2+3
	}else{
		s = "-"
		start = pos2-3
		stop = pos1+0
	}
#	if(start>f1 && start<f2 && stop>f1 && stop<f2)
		print parent, start, stop, "ORF-"num, 0, s
}
