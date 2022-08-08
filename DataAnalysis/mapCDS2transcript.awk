#!/usr/bin/mawk -f

BEGIN{
	OFS ="\t"
}

NR==FNR{
	if($0~/>/)
		h = substr($1,2)
	else
		RNA[h] = RNA[h] $0
	next
}

/>/{
	h2 = substr($1,2)
	parent[h2] = $0
	sub(/.*parent=/,"", parent[h2])
	sub(/; .*/,"", parent[h2])
	sub(/.*,/,"", parent[h2])
#	print h2, substr(parent,RSTART+1)
	next
}

{
	cds[h2] = cds[h2] $0
}

END{
	for(i in parent){
#		match(RNA[parent[i]],cds[i])
#		print parent[i], RSTART, RLENGTH, i
		if(parent[i] in RNA)
			print ">" parent[i] "\n" RNA[parent[i]]
	}
}

