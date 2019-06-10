#(.pdf).eps

#未知領域に侵入していく生物種の挙動について
#離散版拡散ロジスティック方程式

import math
import copy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

a = 1#自然増殖率
b = 1#bu個体数の増加に伴う死滅率
d = 1#拡散率
drealx = 1#刻み幅
drealt = 0.15#刻み時間

L = 37#L*Lマス
k = 0.2#境界付近の死滅率
dx = math.sqrt(a / d) * drealx
dt = a * drealt
p = dt/(dx ** 2)
#セルを作る、w = b/a*u
w = np.zeros((L, L))
#セルを作る
uh = np.zeros((L, L))
#初期分布の設定
w[L // 2][L//2] = 0.3
w[L // 2, L // 2 + 1] = 0.3
w[L // 2, L // 2 - 1] = 0.3
w[L // 2, L // 2 + 2] = 0.3
w[L // 2, L // 2 - 2] = 0.3
w[L // 2, L // 2 + 3] = 0.3
w[L // 2, L // 2 - 3] = 0.3
w[L // 2, L // 2 + 4] = 0.3
w[L // 2, L // 2 - 4] = 0.3
#w[L // 8 - 1][L // 2 + 1] = 0.8
#w[L // 8 + 1][L // 2 - 1] = 0.8
#w[L // 8 - 1][L // 2 - 1] = 0.8
#w[L // 8 + 1][L // 2 + 1] = 0.8
#移動ベクトル
dd = ((0, -1), (-1, 0), (0, 1), (1, 0))
#ステップ数の設定
stepcount = 1
#指定したstepcountの回数まで生物種の移動を実行
while stepcount <= 200:
    wc = copy.deepcopy(w)#同時更新を実現するため前の状態を保持しておく
    for i in range(L):
        for j in range(L):
            if i == 0 or i == L - 1 or j == 0 or j == L - 1:#端の部分は0とする
                w[i][j] = 0
            else:
                if w[i][j] > 0:#生物種が存在しているとき
                    w[i][j] = wc[i][j] * (1 - 4 * p + (1 - wc[i][j]) * dt)
                for dx, dy in dd:#下、左、上、右の四方向から生物種が流れこむ、もしくは流れ出る（濃度の違いによって決定される、濃度が高いほうから低いほうに流れる）
                    nx = i + dx
                    ny = j + dy
                    if 0 <= nx < L and 0 <= ny < L:#L*Lの幅に収まっているときのみ
                        w[i][j] += p * wc[nx][ny]
                if wc[i][j] == 0:#生物種が存在していないときはkの死滅率を受ける（自由境界以外のところではw[i][j]=0となり、w[i][j] = max(0, - k)で結局0となるので、自由境界のところとそれ以外の、種が存在しないところを一緒に考えても良い）
                    if (a / b) * w[i][j] - k <= 0:
                        uh[i][j] += (a / b) * w[i][j]
                        if uh[i][j] > k:
                            w[i][j] = uh[i][j] - k
                            uh[i][j] = 0
                        else:
                            w[i][j] = 0
                    else:
                        w[i][j] = (a / b) * w[i][j] - k

    if stepcount % 40 == 0 or stepcount == 1:
        plt.figure()
        sns.heatmap(w, square=True, cmap='Reds', vmax = 1, vmin = 0)
        #sns.heatmap(uh, square=True, cmap='Reds', vmax=1, vmin=0)
        plt.show()

        #plt.draw()
        #plt.show()
        #plt.pause(1)
        #plt.savefig("sffa{stp:03}".format(stp=stepcount))
        #plt.close
        plt.close('all')
    stepcount += 1  # 試行回数を１増やす