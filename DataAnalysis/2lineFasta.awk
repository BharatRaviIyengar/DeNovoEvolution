#!/usr/bin/mawk -f

BEGIN{
	FS = "\n"
	RS = ">"
	OFS = "\n"
	ORS = "\n"
}

{
	if(FNR==1)
		next
	h1 = $1
	seq = toupper(substr($0,length(h1)+2))
	gsub("\n","",seq)
	print ">"h1, seq
}
