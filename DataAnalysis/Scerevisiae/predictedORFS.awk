#!/usr/bin/mawk -f

BEGIN{
    FS = OFS = "\t"
    patg = 0.016850003370000675
    pstp[0] = 0.059240958727919374
    pstp[1] = 0.03987753248113172
    pstp[2] = 0.048208999659451635
}

{
    len = $5-$4
    for(i=90;i<=len-3;i+=3){
        norf = len - i -2
        k = i/3
        for(f in pstp){
            porf = pstp[f]*patg*(1-pstp[f])^(k-2)
            eorf[f]+= norf*porf
        }
        
    }
}

END{
    for(i in eorf)
        print i, eorf[i]
}

