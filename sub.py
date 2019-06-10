import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
# 正方行列と X および Y のラベルの行列を渡す
"""def draw_heatmap(data, row_labels, column_labels):
    # 描画する
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data, cmap=plt.cm.Blues)

    ax.set_xticks(np.arange(data.shape[0]) + 0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1]) + 0.5, minor=False)

    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)
    plt.show()
    #plt.savefig('image.png')

    return heatmap"""
#draw_heatmap(data, row_labels, column_labels)
""""# generate data
x = np.random.rand(100)
y = np.random.rand(100)

fig = plt.figure()

ax = fig.add_subplot(1,1,1)

ax.scatter(x,y)

ax.set_title('first scatter plot')
ax.set_xlabel('x')
ax.set_ylabel('y')
for i in range(4):
    plt.show()

fig = plt.figure()
#fig.add_subplot(1,1,1)
fig.add_subplot(3,4,10)
#plt.show()
df_flights = sns.load_dataset('flights')
df_flights.head(5)
df_flights_pivot = pd.pivot_table(data=df_flights, values='passengers',
                                  columns='year', index='month', aggfunc=np.mean)
plt.figure(figsize=(12, 9))
sns.heatmap(df_flights_pivot, annot=True, fmt='g', cmap='Blues')"""
"""list_2d = [[0, 6, 2], [3, 4, 5]]
plt.figure()
sns.heatmap(list_2d)
plt.show()
plt.close('all')"""
