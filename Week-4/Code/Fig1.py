fig1, ax = plt.subplots(figsize=(8,8))
ax.semilogy(x, exp(x), 'k', label="Original Function")
ax.semilogy(x, exp(x % (2*np.pi)), 'b--', label="Expected Function from fourier series")
ax.semilogy(x_, exp_fit, 'go', label="Function from Least Square fit Co-efficients")
ax.legend()
ax.set_title(r'Figure 1 : Plot of $e^{x}$')
ax.set_xlabel(r'x$\longrightarrow$')
ax.set_ylabel(r'$e^{x}\longrightarrow$')
ax.grid(True)
fig1.savefig(PATH + "Figure1.png")
