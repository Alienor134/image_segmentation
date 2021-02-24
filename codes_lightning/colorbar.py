import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(20,12))
for i, image in enumerate(ims):
    plt.subplot(2, 4, i + 1)
    plt.axis("off")
    plt.imshow(image, aspect='auto')
    plt.subplots_adjust(hspace=0.05, wspace=0.05)
    plt.show()


folder = "F://image_analysis/"
cop = np.copy(ims[0])
R = cop[:,:,0]
G = cop[:,:,1]
B = cop[:,:,2]

R[:, :] = -2.2
R[0, 0] = -0.5

G[:,:] = -3.2
G[0,0] = -1.3

B[:,:] = -2.8
B[0,0] = -0.5

plt.close('all')
fig = plt.figure()
ax1 = fig.add_subplot(111)

imax1 = ax1.imshow(R,cmap="Reds")#plot
cbar = plt.colorbar(imax1, extend='neither', spacing='proportional',
                orientation='horizontal', shrink=1, format="%.01f")
cbar.set_label(r"tau 2 high", size=15)
cbar.ax.tick_params(labelsize=15) 
plt.savefig(folder + "/RGB/R.tiff")
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)

imax1 = ax1.imshow(G,cmap="Greens")#plot
cbar = plt.colorbar(imax1, extend='neither', spacing='proportional',
                orientation='horizontal', shrink=1, format="%.01f")
cbar.set_label(r"tau 1 high", size=15)
cbar.ax.tick_params(labelsize=15) 
plt.savefig(folder + "/RGB/G.tiff")
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)

imax1 = ax1.imshow(B,cmap="Blues")#plot
cbar = plt.colorbar(imax1, extend='neither', spacing='proportional',
                orientation='horizontal', shrink=1, format="%.01f")
cbar.set_label(r"tau 1 low", size=15)
cbar.ax.tick_params(labelsize=15) 
plt.savefig(folder + "/RGB/B.tiff")
plt.show()