# Rosalee Hayes, Skye Russ, Simon Buan, Anders Stall

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import OneHotEncoder


# Read the data
data = pd.read_csv('condensed_cvr.csv', dtype=str).sample(n=20000, random_state=42)

# Set aside some data for testing
train_set, test_set = train_test_split(data, test_size=0.2, shuffle=True, random_state=42)
# Use 42 as the random number seed.


# Separate the target column as the label for supervised learning
# Uncomment one of the two lines below to predict a specific race

# TARGET_COLUMN = 'City of Portland, Mayor'
TARGET_COLUMN = 'President and Vice President'

train_X = train_set.drop(TARGET_COLUMN, axis=1)
train_y = train_set[TARGET_COLUMN].copy()


# Use one-hot encoding for the input attributes, as sklearn needs these
hot_encoder = OneHotEncoder(handle_unknown='ignore')
train_X_1hot = hot_encoder.fit_transform(train_X)
# may be votes that appear in the test set but not in the training set. Call the augmented train_X "train_X_1hot".


# Train a decision tree
tree = DecisionTreeClassifier(max_depth=6)
tree = tree.fit(train_X_1hot, train_y)

# Test it on the training set
print(f'Accuracy on training set: {tree.score(train_X_1hot, train_y)}')


# Test it on the test set
test_X = test_set.drop(TARGET_COLUMN, axis=1)
test_y = test_set[TARGET_COLUMN].copy()
test_X_1hot = hot_encoder.transform(test_X)

# test data.
print(f'Accuracy on test set: {tree.score(test_X_1hot, test_y)}')
