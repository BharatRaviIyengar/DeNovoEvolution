#!/usr/bin/gawk -f

BEGIN{
    FS = OFS = "\t"
}

{
    for(s=50;s<=100;s+=10){
        if($10>=s){
            counts[s][$9]++
            xolen = $6*$10/100
            olen[s][$9] += xolen
            tlen[s] += xolen
            total[s]++
        }
    }
    
} 
    
END{
    print "MinOverlap%","Frame","TotalORF","CountORF","OverlapLen","ORF%Frm", "AvgLen", "Overlap%Frm"
    for(s=50;s<=100;s+=10){
        for(i in counts[s])
            print s, i, total[s], counts[s][i], sprintf("%0.0f",olen[s][i]),  sprintf("%0.2f",100*counts[s][i]/total[s]),  sprintf("%0.2f",olen[s][i]/counts[s][i]),sprintf("%0.2f",100*olen[s][i]/tlen[s])
    }
}