fig6, ax = plt.subplots(figsize=(8,8))
ax.loglog(np.abs(ccos_cf[1::2]), 'yo', label=r"$a_{n}$ using Integration")
ax.loglog(np.abs(ccos_cf[2::2]), 'ro', label=r"$b_{n}$ using Integration")
ax.loglog(np.abs(ccos_cf_fit[1::2]), 'go', label=r"$a_{n}$ using Least Squares")
ax.loglog(np.abs(ccos_cf_fit[2::2]), 'bo', label=r"$b_{n}$ using Least Squares")
ax.legend(loc='upper right')
ax.set_title("Figure 5 : Fourier coefficients of $\cos(\cos{x})$ (Log-Log)")
ax.set_xlabel(r'n$\longrightarrow$')
ax.set_ylabel(r'Magnitude of coeffients$\longrightarrow$')
ax.grid(True)
fig6.savefig(PATH + 'Figure6.png')
