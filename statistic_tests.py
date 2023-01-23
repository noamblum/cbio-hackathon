import numpy as np
from scipy.stats import chi2_contingency

# Create the contingency table
data = np.array([[10, 5, 3], [8, 4, 2], [12, 6, 4]])

# Perform the chi-squared test
chi2, p, dof, expected = chi2_contingency(data)

# Print the results
print("Chi-squared test statistic:", chi2)
print("p-value:", p)
