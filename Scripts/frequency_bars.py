# ------------------------------------------------------------------------------------------------------
# A module containing a function for use in plotting frequency bar charts
# ------------------------------------------------------------------------------------------------------
# Import

import matplotlib.pyplot as plt


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to plot a bar chart with necessary elements as arguments

def freq_bars(brackets, labels, xpos, yscale, xlabel, ylabel):
    # Set up bar chart with correct bars
    plt.bar(xpos, brackets, edgecolor="black")
    # Set bar labels
    plt.xticks(xpos, labels, rotation=10)
    # Set y scale type
    plt.yscale(yscale)
    # Set x and y labels
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # Draw the plot
    plt.show()
