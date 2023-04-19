#!/usr/bin/mawk -f

BEGIN{
	FS = OFS = "\t"
	while((getline < "novel.tsq") > 0){
		rna[$1] = $2
	}
	stpcdn["TAA"]
	stpcdn["TGA"]
	stpcdn["TAG"]
}

/>/{
	split($0,x," ")
	orfname = parent = substr(x[1],2)
	sub(/_[1-9]+$/,"",parent)
	split(substr(x[1],2),y,/[_)(]/)
	chr = y[1]"_"y[2]
	strand = y[4]
	orfnum = y[6]
	split(y[3],gpos,"-")
	orfpos1 = 0 + substr(x[2],2)
	orfpos2 = 3 + x[4]

	last3 = toupper(substr(rna[parent],x[4]+1,3))
	frst3 = toupper(substr(rna[parent],orfpos1,3))
	if(!(last3 in stpcdn && frst3 == "ATG"))
		next
	
	if(strand=="+"){
		strt = gpos[1]+orfpos1
		stop = gpos[1]+orfpos2
	}
	if(strand == "-"){
		strt = gpos[2]-orfpos2+1
		stop = gpos[2]-orfpos1+1
	}
	if(stop >= gpos[2])
		next
	print chr, "getorf", "CDS", strt, stop, stop-strt+1, strand, orfname
}