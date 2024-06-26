fig2, ax = plt.subplots(figsize=(8,8))
ax.plot(x, ccos(x), 'k', linewidth=2, label="Original Function")
ax.plot(x, ccos(x % (2*np.pi)), 'r--', label="Expected Function from fourier series")
ax.plot(x_, ccos_fit, 'go', label="Functions from Least Square fit Co-efficients")
ax.legend(loc='upper right')
ax.set_title(r"Figure 2 : Plot of $\cos(\cos(x))$")
ax.set_xlabel(r'x$\longrightarrow$')
ax.set_ylabel(r'$\cos(\cos(x)\longrightarrow$')
ax.grid(True)
fig2.savefig(PATH + "Figure2.png")
