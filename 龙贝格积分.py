import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.abc import x,y

y = sp.sin(x) / x #在此处输入函数表达式(sympy函数库)
low = 0
high = 1 #积分上下界
h = high - low
ep = 0.002 #误差界要求
#画函数图像
X = np.linspace(low,high,500)
Y = np.empty(500)
for i in range(498):
    Y[i+1] = y.evalf(subs = {x:X[i+1]})
Y[0] = sp.limit(y,x,X[0])
Y[499] = sp.limit(y,x,X[499])
plt.plot(X,Y)
plt.scatter(X[0],Y[0],color = 'red')
plt.scatter(X[-1],Y[-1],color = 'red')
plt.text(X[0],Y[0],(round(X[0],5),round(Y[0],5)),fontsize = 16)
plt.text(X[-1],Y[-1],(round(X[-1],5),round(Y[-1],5)),fontsize = 16)

T = [h * (Y[0] + Y[-1]) / 2]
S = []
C = []
R = []

for i in range(30): #龙贝格积分
    TT = 0
    for j in range(2 ** i):
        TT = TT + y.evalf(subs = {x:(low + (2 * j + 1) * h / (2 ** (i+1)))})
    TT = TT / (2 ** (i+1)) * h + T[i] / 2
    T.append(TT)
    S.append(T[i + 1] + (T[i + 1] - T[i]) / 3)
    if (i <= 0):
        continue
    C.append(S[i] + (S[i] - S[i - 1]) / 15)
    if (i <= 1):
        continue
    R.append(C[i - 1] + (C[i - 1] - C[i - 2]) / 63)
    if (i >= 3):
        if (abs(R[i - 2] - R[i - 3]) < ep):
            break;

print('integral =',R[-1])
plt.text((X[0] + X[-1]) / 2,(Y[0] + Y[-1]) / 2,'integral='+str(R[-1]),fontsize = 16)
plt.show()