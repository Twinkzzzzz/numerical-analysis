import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.abc import x,y
yp = sp.Symbol('z')

yp = x * y * y #导数表达式
low = 0
high = 1 #求解范围
initco = 0 #初始条件(坐标)
initval = 1 #初始条件(值)
h = 0.1 #步长(带方向,负数代表沿数轴负方向计算)
X = [initco]
Y = [initval]

xi = initco
yi = initval
plt.scatter(xi,yi,color = 'red')
plt.text(xi,yi,(round(xi,5),round(yi,5)),fontsize = 16)
i = 0

while((xi + h) >= low and (xi + h) <= high): #R-K法主程序
    i = i + 1
    K1 = h * yp.evalf(subs={x: xi,y: yi})
    K2 = h * yp.evalf(subs={x: xi + h / 2,y: yi + K1 / 2})
    K3 = h * yp.evalf(subs={x: xi + h / 2, y: yi + K2 / 2})
    K4 = h * yp.evalf(subs={x: xi + h, y: yi + K3})
    yi = yi + (K1 + 2 * K2 + 2 * K3 + K4) / 6
    xi = i * h + initco
    print(round(K1,3), round(K2,3), round(K3,3), round(K4,3), round(xi,3), round(yi,3))
    X.append(xi)
    Y.append(yi)

plt.scatter(xi, yi, color='red')
plt.text(xi,yi,(round(xi,5),round(yi,5)),fontsize = 16)
h = int((len(X) - 1) / 4)
if (h == 0):
    h = 1
for i in range(3): #除边界外再取至多三个点画图
    if ((i + 1) * h >= len(X)-1):
        break
    xi = X[(i + 1) * h]
    yi = Y[(i + 1) * h]
    plt.scatter(xi, yi, color='red')
    plt.text(xi,yi,(round(xi,5),round(yi,5)),fontsize = 16)

for i in range(len(X)):
    print((X[i],Y[i]))
plt.plot(X,Y)
plt.show()