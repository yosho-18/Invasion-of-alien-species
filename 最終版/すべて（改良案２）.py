#(.pdf).eps

#未知領域に侵入していく生物種の挙動について
#離散版拡散ロジスティック方程式

from initialdistribution import *
from otherfanction import *
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

L = 51#L*Lマス
k1 = 0.001#もうひとつの死滅率（履歴なし）（絶滅率）
k2 = 0.138#境界付近の死滅率0.4～0.2（履歴率）
dx = math.sqrt(a / d) * drealx
dt = a * drealt
p = dt/(dx ** 2)
#セルを作る、w = b/a*u
w = np.zeros((L, L))
#セルを作る
uh = np.zeros((L, L))
#初期分布の設定
#center_cross(w, L)#
#cross(w, L)
#square(w,L)#1,72
#rectangle(w, L)#4,120
#interval_square(w, L)#6,144
#I(w, L)#2,104
S(w, L)#5,128
#T(w, L)#2,104
#overdensity(w, L, a, b)
title = "S"#cross, square, rectangle, interval_square, I, S, T, overdensity

#隣接度
#print(adjacent_degree(w, L))

#移動ベクトル
dd = ((0, -1), (-1, 0), (0, 1), (1, 0))

stepcount = 0
heatmap(w, stepcount, a, b, title)
#ステップ数の設定
stepcount += 1
ascount1 = 0
ascount2 = 0
vanicount = 0
distcount = 0
fn = 200

iok = 0

scli = []
hli = []
uli = []
scli.append(0)
hli.append(L // 2 + 2)#L // 2 + 2

uli.append((a / b) * w[L // 2][L // 2])

allarea = (L - 2) ** 2
redarea = 0
areali = [9 / allarea]

arearootli = []
arearootli.append(math.sqrt(9 / allarea))

denallarea = (L - 2) ** 2 * (a / b)
denredarea = 0
denareali = [9 / allarea * 0.3]

denarearootli = [math.sqrt(9 / allarea * 0.3)]


#指定したstepcountの回数まで生物種の移動を実行
while stepcount <= fn:
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
                    w[i][j] = max(0, (a / b) * w[i][j] - k1)
                    if (a / b) * (w[i][j] + uh[i][j]) - k2 <= 0:#密度がkに達しないとき
                        uh[i][j] += (a / b) * w[i][j]#uhにuの履歴を保存しておく
                        if uh[i][j] > k2:#uhがkを越したら
                            w[i][j] = uh[i][j] - k2#wに密度を入れる
                            uh[i][j] = 0#uhをゼロにしておく
                        else:#
                            w[i][j] = 0
                    else:#
                        w[i][j] = (a / b) * (w[i][j] + uh[i][j]) - k2
                        uh[i][j] = 0  # uhをゼロにしておく
    """"#拡散しきったら終了
    for i in range(L):
        for j in range(L):
            if w[i][j] > 0:
                iok += 1
    if iok == (L - 2) * (L - 2):
        print("拡散", stepcount, min(uli), scli[uli.index(min(uli))])
        heatmap(w, stepcount, a, b, title)
        graph(scli, hli, 'h position', stepcount)
        graph(scli, uli, 'u density', stepcount)
        graph(scli, areali, 'area velocity', stepcount)
        graph(scli, arearootli, 'area velocity root', stepcount)
        graph(scli, denareali, 'area velocity density', stepcount)
        graph(scli, denarearootli, 'area velocity root density', stepcount)
        exit()
    else:
        iok = 0"""

    #自由境界速度
    for i in range(1, L):
        if w[i - 1][L // 2] > 0 and w[i][L // 2] == 0:
            #scli.append(stepcount)
            hli.append(i)
    #scli.append(stepcount)
    #密度
    scli.append(stepcount)
    uli.append((a / b) * w[L // 2][L // 2])

    if stepcount % (fn // 5) == 0 or stepcount == 200:
        heatmap(w, stepcount, a, b, title)
    #面積速度を求める，（面積だけのと，濃度も考慮するのの二つ）
    redarea = 0
    denredarea = 0
    for i in range(L):
        for j in range(L):
            if (a / b) * w[i][j] > 10 ** (-4):
                redarea += 1
                denredarea += (a / b) * w[i][j]
    areali.append(redarea / allarea)
    arearootli.append(math.sqrt(redarea / allarea))
    denareali.append(denredarea / denallarea)
    denarearootli.append(math.sqrt(denredarea / denallarea))

    #100回前に密度のある面積を数えておく
    if stepcount == fn // 10 * 9:
        for i in range(L):
            for j in range(L):
                if (a / b) * w[i][j] > 10 ** (-4):
                    ascount1 += 1
    """"#stepcount=400,1400のときのhの位置
    if stepcount == 100 or stepcount == 350:
        print(stepcount, hli[-1])"""
    #終了時に密度のある面積を数えておく
    if stepcount == fn:
        for i in range(L):
            for j in range(L):
                if (a / b) * w[i][j] > 10 ** (-4):
                    ascount2 += 1
                    distcount += 1
                    vanicount += 1
        sprevanistay(vanicount, ascount1, ascount2, distcount, L)
        print(stepcount, min(uli), scli[uli.index(min(uli))])
        print(stepcount, max(uli), scli[uli.index(max(uli))])

        heatmap(w, stepcount, a, b, title)
        graph(scli, hli, 'h position', stepcount)
        graph(scli, uli, 'u density', stepcount)
        """graph(scli, areali, 'area velocity', stepcount)#redarea / (L - 2) * (L - 2)
        graph(scli, arearootli, 'area velocity root', stepcount)
        graph(scli, denareali, 'area velocity density', stepcount)
        graph(scli, denarearootli, 'area velocity root density', stepcount)"""
        #plt.draw()
        #plt.show()
        #plt.pause(1)
        #plt.savefig("sffa{stp:03}".format(stp=stepcount))
        #plt.close
        plt.close('all')
        """for i in range(L):
        for j in range(L):
            if i == 1 or i == L - 2 or j == 1 or j == L - 2:#端の部分は0とする
                if w[i][j] > 0:  # 壁に当たったら終了
                    print(stepcount)
                    heatmap(w, stepcount, a, b, title)
                    graph(scli, hli, 'h position', stepcount)
                    graph(scli, uli, 'u density', stepcount)
                    graph(scli, areali, 'area velocity', stepcount)
                    graph(scli, arearootli, 'area velocity root', stepcount)
                    graph(scli, denareali, 'area velocity density', stepcount)
                    graph(scli, denarearootli, 'area velocity root density', stepcount)
                    exit()"""
    """if w[L - 2][L // 2] > 0:#壁に当たったら終了
        print(stepcount)
        heatmap(w, stepcount, a, b, title)
        #graph(scli, hli, 'h position', stepcount)
        #graph(scli, uli, 'u density', stepcount)
        graph(scli, areali, 'area velocity', stepcount)
        graph(scli, arearootli, 'area velocity root', stepcount)
        graph(scli, denareali, 'area velocity density', stepcount)
        graph(scli, denarearootli, 'area velocity root density', stepcount)
        exit()"""
    stepcount += 1  # 試行回数を１増やす