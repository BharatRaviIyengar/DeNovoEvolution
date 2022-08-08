#!/usr/bin/mawk -f

BEGIN{
	FS = OFS = "\t"
}

FNR>1{
	start = $7
	stop = $8
	strand = "+"
	if($7>$8){
		start = $8
		stop = $7
		strand = "-"
	}
	print $2, "BLASTN", "transcript", start, stop, 0, strand, 0, $1 > $3"_dmel-6.44-ncRNAs.gtf"
}
