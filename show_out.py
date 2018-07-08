import re
import matplotlib.pyplot as plt

file_name = 'Cd_collidetype1_bbtype1_100_0.8.out'
f = open(file_name,'rU')
num = 0.0
time = []
Re = []
Cd1 = []
Cd2 = []
Cd = []
for line in f:
    if num>0.5:
        line_str = ''.join(line)
        temp = re.split(r' +',line_str)
        time.append(temp[1])
        Re.append(temp[2])
        Cd1.append(temp[3])
        Cd2.append(temp[4])
        Cd.append(temp[5])
        
    num =num + 1 
plt.loglog(Re,Cd1)
plt.show()
