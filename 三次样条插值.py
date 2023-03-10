import numpy as np
import sys
import json
import matplotlib.pyplot as plt

#程序从对应json文件读取参数
#'nodes'代表插值节点
#'boundtype'代表边界条件类型，其中2代表已知边界二阶导数值，1代表已知边界一阶导数值，0代表周期边界条件
#'boundval'代表边界条件的值，当'boundtype'为0时无效
file = open('sanciyangtiaochazhi.json','r')
data = json.load(file)
file.close()

nodes = data['nodes']
type = data['boundtype']
n = len(nodes) - 1
nodes.append(nodes[1].copy())
nodes[n + 1][0] = nodes[n][0] + nodes[1][0] - nodes[0][0]
miu = []
lamda = []
d = []
M = []

if (n <= 0): #仅有一个插值节点
    print("statistic not enough")
    sys.exit(0)
elif (type == 0 and nodes[0][1] != nodes [n][1]): #是周期边界条件但边界值不符合要求
    print("statistic not available")
    sys.exit(0)

for i in range(n):
    x1 = nodes[i][0]
    y1 = nodes[i][1]
    x2 = nodes[i + 1][0]
    y2 = nodes[i + 1][1]
    x3 = nodes[i + 2][0]
    y3 = nodes[i + 2][1]
    h1 = x2 - x1
    h2 = x3 - x2
    miu.append(h1 / (h1 + h2))
    lamda.append(h2 / (h1 + h2))
    if (i < n - 1):
        d.append(((y3 - y2) / h2 - (y2 - y1) / h1) / (h1 + h2) * 6)

if (type == 2): #第一类边界条件
    if (n > 1):
        d[0] = d[0] - miu[0] * data['boundval'][0]
        d[-1] = d[-1] - lamda[-1] * data['boundval'][1]
        B = np.zeros([n - 1,n - 1])
        for i in range(n - 1):
            if (i == 0):
                B[0][0] = 2
                if (n > 2):
                    B[0][1] = lamda[0]
            elif (i == n - 2):
                if (n > 2):
                    B[i][n - 3] = miu[n - 2]
                B[i][n - 2] = 2
            else:
                B[i][i - 1] = miu[i]
                B[i][i] = 2
                B[i][i + 1] = lamda[i]
        D = np.asarray(d,dtype = float)
        Binv = np.linalg.inv(B)
        M = np.dot(Binv, D)
    M = np.insert(M,0,data['boundval'][0])
    M = np.append(M,data['boundval'][1])
elif (type == 1): #第二类边界条件
    d.insert(0,6 * ((nodes[1][1] - nodes[0][1]) / (nodes[1][0] - nodes[0][0]) - data['boundval'][0]) / (nodes[1][0] - nodes[0][0]))
    d.append(6 * (data['boundval'][1] - (nodes[n][1] - nodes[n - 1][1]) / (nodes[n][0] - nodes[n - 1][0])) / (nodes[n][0] - nodes[n - 1][0]))
    B = np.zeros([n + 1,n + 1])
    B[0][0] = 2
    B[0][1] = 1
    B[n][n - 1] = 1
    B[n][n] = 2
    for i in range(n - 1):
        B[i + 1][i] = miu[i]
        B[i + 1][i + 1] = 2
        B[i + 1][i + 2] = lamda[i]
    D = np.asarray(d, dtype=float)
    Binv = np.linalg.inv(B)
    M = np.dot(Binv, D)
else: #周期性边界条件
    h1 = nodes[1][0] - nodes[0][0]
    hn = nodes[n][0] - nodes[n-1][0]
    y1 = nodes[1][1]
    y0 = nodes[0][1]
    ynm = nodes[n - 1][1]
    d.append(6 / (h1 + hn) * ((y1 - y0) / h1 - (y0 - ynm) / hn))
    B = np.zeros([n,n])
    for i in range(n):
        if (i == 0):
            B[0][0] = 2
            if (n > 2):
                B[0][1] = lamda[0]
        elif (i == n - 1):
            if (n > 2):
                B[i][n - 2] = miu[n - 1]
            B[i][n - 1] = 2
        else:
            B[i][i - 1] = miu[i]
            B[i][i] = 2
            B[i][i + 1] = lamda[i]
    B[0][-1] = B[0][-1] + miu[0]
    B[-1][0] = B[-1][0] + lamda[-1]
    print(B)
    D = np.asarray(d, dtype=float)
    Binv = np.linalg.inv(B)
    M = np.dot(Binv, D)
    M = np.insert(M, 0, M[-1])

plt.grid()
for i in range(n + 1):
    plt.scatter(nodes[i][0],nodes[i][1],color = 'red')
    plt.text(nodes[i][0],nodes[i][1],(round(nodes[i][0],5),round(nodes[i][1],5)),fontsize = 16)

for i in range(n):
    xi = nodes[i + 1][0]
    xim = nodes[i][0]
    yi = nodes[i + 1][1]
    yim = nodes[i][1]
    Mi = M[i + 1]
    Mim = M[i]
    hi = xi - xim
    X = np.linspace(xim,xi,50)
    Y = np.zeros([50])
    a3 = round((Mi - Mim) / 6 / hi,5) #记录系数,便于查看保留五位小数
    a2 = round((Mim * xi - Mi * xim) / 2 / hi,5)
    a1 = round((Mi * xim * xim / 2 - Mim * xi * xi / 2 - (yim - hi * hi * Mim / 6) + (yi - hi * hi * Mi / 6)) / hi,5)
    a0 = round((Mim * (xi ** 3) / 6 - Mi * (xim ** 3) / 6 + xi * (yim - hi * hi * Mim / 6) - xim * (yi - hi * hi * Mi /6)) / hi,5)
    print([xim,xi],':',a3,'* x^3 +',a2,'* x^2 +',a1,'* x +',a0)
    for j in range(50): #M数组求解完毕，带入三次样条插值函数
        item1 = (xi - X[j]) * (xi - X[j]) * (xi - X[j]) / 6 / hi * Mim
        item2 = (X[j] - xim) * (X[j] - xim) * (X[j] - xim) / 6 / hi * Mi
        item3 = (yim - hi * hi * Mim / 6) * (xi - X[j]) / hi
        item4 = (yi - hi * hi * Mi / 6) * (X[j] - xim) / hi
        Y[j] = item1 + item2 + item3 + item4
    plt.plot(X,Y,color = 'black')
plt.show()