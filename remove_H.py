#-*-coding:utf-8-*-
import sys
if len(sys.argv)<2:
    print('python remove_H -h for help')
    sys.exit(1)
remove_water=False
for i in sys.argv:
    if "-h" in i:
        print("'python remove_H yourpdb.pdb' to remove H from PDB, with '-water' for removing water from PDB")    
        sys.exit(0)
    if "-water" in i:
        remove_water=True
s=[]
ids=1
with open(sys.argv[1]) as reader, open("new.pdb",'w') as writer:
    for index,line in enumerate(reader):
        if  line.startswith("ATOM") or line.startswith("TER"):
            if line[11:].lstrip().startswith('H'):
                continue
            tmpname=line.split()[0]
            s.append("%-6s%5d%s" %(tmpname,ids,line[11:]))
            ids+=1
        elif  line.startswith("HETATM"):
            if line[17:].lstrip().startswith('HOH') and remove_water==True:
                continue
            tmpname=line.split()[0]
            s.append("%-6s%5d%s" %(tmpname,ids,line[11:]))
            ids+=1            
        elif "CONECT" in line:
            continue
        elif "MASTER" in line:
            continue        
        else:
            s.append(line)
    writer.writelines(s)