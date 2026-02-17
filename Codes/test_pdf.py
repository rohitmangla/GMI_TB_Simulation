import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
np.random.seed(0)
x= np.random.randint(1,10,30)
y= x+np.random.normal(0,1,30)

ax = sns.regplot(x,y)

plt.show()
