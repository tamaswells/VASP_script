#!/usr/bin/python

with open("new",'r') as reader,open('INCAR.txt','w') as writer:
    read_lines=reader.readlines()
    for i in range(0,len(read_lines),2):
        read_lines[i+1]=read_lines[i+1].replace('false','.FALSE.').replace('true','.TRUE.')
        writer.write(read_lines[i].strip('\r\n').strip('\n')+" = "+read_lines[i+1])
