
N = 10000
with open('batch.txt', 'r') as f:
    j = int(f.read())
with open('batch.txt', 'w') as w:
    w.write(str(j+1))

with open("ast.txt", "w") as ast:
    for i in range(j*N+1, (j+1)*N+1):
        ast.write("'"+f'{i:07}'+"'"+"\n")
