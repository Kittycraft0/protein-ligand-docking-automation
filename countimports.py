# 4/17/2026

import numpy as np

list=[]

exit=0
while not exit:
    inputt=input("input next")
    if inputt=="exit":
        break
    list.append(inputt)

print(list)
print(".")
print(np.unique(list, return_counts=True))
uniquelist=np.unique(list, return_counts=True)
i=0
while i<len(uniquelist[0]):
    print(f"{uniquelist[1][i]}: {uniquelist[0][i]}")
    i+=1