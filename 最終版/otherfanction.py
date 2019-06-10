import matplotlib.pyplot as plt
import seaborn as sns
#関数

#ヒートマップ
def heatmap(w, stepcount, a, b, title):
    plt.figure()
    ax = plt.axes()
    sns.heatmap(w, square=True, cmap='Reds', vmax=a / b, vmin=0, ax=ax)
    ax.set_title(title)
    # sns.heatmap(uh, square=True, cmap='Reds', vmax = 1, vmin = 0)
    #plt.show()
    plt.savefig(title + "_heatmap{stp:03}".format(stp=stepcount))

#グラフ
def graph(scli, hli, ysym, stepcount):
    plt.figure()
    plt.plot(scli, hli)
    plt.xlabel('stepcount')
    plt.ylabel(ysym)
    # plt.xlim(1300,)
    #plt.show()
    plt.savefig(ysym + "{stp:03}".format(stp=stepcount))

#拡散消滅停留　判定
def sprevanistay(vanicount, ascount1, ascount2, distcount, L):
    if vanicount == 0:
        return print("消滅")
    elif ascount1 < ascount2:
        return print("拡散")
    elif ascount1 == ascount2 and distcount == (L - 2) * (L - 2):
        return print("拡散")
    elif ascount1 == ascount2 and distcount != (L - 2) * (L - 2) and distcount != 0:
        return print("停留")
    elif ascount1 > ascount2:
        return print("消滅")

#隣接度
def adjacent_degree(w, L):
    adj = 0
    for i in range(L):
        for j in range(L):
            if w[i][j] > 0:
                for k in range(10):
                    for l in range(k):
                        x = i + (k - l) - 1
                        y = j + l + 1
                        #abs(x - i) + abs(y - j) = k
                        if w[x][y] > 0:
                            adj += k
                        x = i - (k - l) + 1
                        y = j - l - 1
                        # abs(x - i) + abs(y - j) = k
                        if w[x][y] > 0:
                            adj += k
                        x = i + (k - l)
                        y = j - l
                        # abs(x - i) + abs(y - j) = k
                        if w[x][y] > 0:
                            adj += k
                        x = i - (k - l)
                        y = j + l
                        # abs(x - i) + abs(y - j) = k
                        if w[x][y] > 0:
                            adj += k
    return adj // 2
