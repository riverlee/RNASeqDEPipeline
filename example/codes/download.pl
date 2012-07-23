$base="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/litesra/SRP/SRP000/SRP000225";

while($id=<DATA>){
    chomp($id);
    print "Download $id ...\n";
    $url=$base."/$id/${id}.lite.sra";
    if (! -e "${id}.lite.sra"){
#        print $url,"\n";
        `wget $url`;
    }
}

__DATA__
SRR002320
SRR002321
SRR002322
SRR002323
SRR002324
SRR002325
