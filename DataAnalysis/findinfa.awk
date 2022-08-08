#!/usr/bin/mawk -f

BEGIN{
	FS = "\n"
	RS = ">"
	OFS = "\t"
	ORS = "\n"
}

NR == FNR{
	if(FNR==1)
		next
	h1 = $1
	seq = toupper(substr($0,length(h1)+2))
#	seq = toupper(substr($0,length(h1)+2,20))
	gsub("\n","",seq)
#	print h1, seq
	x[seq] = h1
	next
}

FNR>1{
	h2 = $1
	h2x = h2
	sub(/\.[0-9]*/,"",h2x)
	sq2 = toupper(substr($0,length(h2)+2))
#	sq2 = toupper(substr($0,length(h2)+2,20))
	gsub("\n","",sq2)
	if(sq2 in x)
		print h2x, x[sq2]
}


