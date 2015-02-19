ofh = None
lastFname = None
for line in open("orig/nbt.3117-S2.filter.txt"):
    fields = line.rstrip("\n").split("\t")
    fname = fields[0]
    if fname=="Cells":
        continue
    if fname!=lastFname:
        ofh = open(fname+".bed", "w")
    ofh.write("\t".join(fields[1:]))
    ofh.write("\n")
    lastFname = fname

