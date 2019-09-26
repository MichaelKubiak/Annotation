import matplotlib.pyplot as plt

def freq_bars(brackets, labels, xpos, yscale, xlabel, ylabel):
    plt.bar(xpos, brackets, edgecolor="black")
    plt.xticks(xpos, labels, rotation=10)
    plt.yscale(yscale)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
