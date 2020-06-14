def readFastqFile (path):
    file = open(path,'r')
    reads = []
    qualities = [] #quality of a read located in the same index as the read
    count = 0
    while True:
        count+=1
        line = file.readline()
        if not line:
                break
        line = line[:-1]
        if((count%2==0) and (count%4)):
            reads.append(line)
        if(count%4==0):
            qualities.append(line)
    file.close()
    return (reads,qualities)

r,q= readFastqFile('SP1.fq')
print(r)

