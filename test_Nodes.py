import re
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

f = open('Nodes.txt')
lines = f.readlines()
i = 1
points = []
for line in lines:
    if i==1:
        R = float(line)
        i = i+1
    elif i==2:
        N = int(line)
        i = i+1
    else:
        nums = re.findall(r"\-?\d+\.?\d*",line)
        points.append([float(num) for num in nums])
        i = i+1
pp = np.array(points)
print(np.shape(pp))
# ax = plt.subplot(111,projection = '3d')
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(pp[:,0],pp[:,1],pp[:,2])
plt.show()