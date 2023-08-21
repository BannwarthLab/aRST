import re
keyw = 'FINAL SINGLE'

with open(f'orca.out', 'r') as f:
    lines = f.readlines()
    line = re.findall(r'\-?\d+\.?\d*', [k for k in lines if keyw in k][0])
    if line and line[0].isdigit():
        print("yeyeee")
        e = line[0]
    else:
        if line:
            print("a")
        if str.isnumeric(line[0]):
            print('b')
        print("nonono")
        e = 'inf'
print(line,line[0])
print('e= ',e)
with open('ENERGY', 'w') as f:
    f.write(e)
