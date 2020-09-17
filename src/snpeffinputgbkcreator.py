import sys
from Bio import SeqIO

try:
    gi=SeqIO.parse("../data/annotation/patho.ids.gbk","genbank")
    go=open("../external/snpEff/data/ipd1060/genes.gbk",'w')
    for rec in gi:
        if rec.features:
            newfeatures=[]
            loctags=[]
            for feature in rec.features:
                if feature.type == "mRNA":
                    try:
                        if not feature.qualifiers["locus_tag"] in loctags:
                            newfeatures.append(feature)
                            loctags.append(feature.qualifiers["locus_tag"])
                        else:
                            print(feature.qualifiers["locus_tag"])
                    except Exception:
                        newfeatures.append(feature)
                else:
                    newfeatures.append(feature)
            rec.features=newfeatures
        SeqIO.write(rec, go, "genbank")
    go.close()
except Exception as e:
    print(e)
    sys.exit(0)
