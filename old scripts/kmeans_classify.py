c = np.zeros((650,19))
for i in range(650):
    hist,bi = np.histogram(p0[i][4],bins=fixed_bin)
    bin_indices = np.digitize(p0[i][4],bi)-1
    if bin_indices[-1] < 18:
        bin_indices[-1] = 18
    c[i] = np.bincount(bin_indices,weights=p0[i][5])

kmeans = KMeans(n_clusters=3, random_state=0, n_init="auto").fit(c1)

e = kmeans.labels_

plt.bar(x=np.where(e==1),height=1,width = 0.9,label='bad scans',color='r')
plt.bar(x=np.where(e!=1),height=1,width = 0.9,label='good scans',color = 'g')
plt.legend()
plt.xlabel("scan number")
plt.show()

plt.bar(x=np.where(int_log_label==0)[0],height=1,width = 0.9,label='bad scans',color='r')
plt.bar(x=np.where(int_log_label!=0)[0],height=1,width = 0.9,label='good scans',color = 'g')
plt.legend()
plt.xlabel("scan number")
plt.show()
plt.savefig('int_log_peaks_filtered')
