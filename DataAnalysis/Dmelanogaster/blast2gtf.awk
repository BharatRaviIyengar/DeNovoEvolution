#!/usr/bin/mawk -f

BEGIN{
	FS=OFS="\t"
}

$5==100{
	if($3>$4){
		orie = "-"
		strt = $4-3
		stop = $3
	}else{
	 	orie = "+"
        strt = $3
        stop = $4+3
	}
	xx[$2 FS "BLASTN" FS "CDS" FS strt FS stop FS "." FS orie FS 0] = $1
}

END{
	for(i in xx)
		print i,xx[i]
}

