import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import BernoulliNB
from sklearn.naive_bayes import MultinomialNB
from sklearn import svm



import pickle

#reading files
df_train = pd.read_csv('plos1sequence.csv')
df_test = pd.read_csv('newextra.csv')


#merging together
frames = [df_train, df_test]
merged_df = pd.concat(frames)

aa = ['A','R','N','D','B','C','E','Q','Z','G','H','I','L','K','M','F','P','S','T','W','Y','V']
encoded_aa = np.arange(len(aa))

subclasses = ['PIP', 'TIP', 'NIP', 'SIP', 'XIP', 'GIP', 'HIP']
encoded_subclasses = np.arange(len(subclasses))


merged_df = merged_df.replace(aa, encoded_aa)
merged_df = merged_df.replace(subclasses, encoded_subclasses)

#separating feature and labels
features = merged_df.drop(['SubClass'], axis=1)
labels = merged_df['SubClass']

features_train, features_test = features[:176], features[176:]
labels_train, labels_test = labels[:176], labels[176:]

# #splitting training and testing data
# features_train, features_test, labels_train, labels_test = train_test_split(features, labels, 
# test_size=0.2)
# i=2
# for l in labels:
#     print(i, l)
#     i+=1

clf = GaussianNB()
clf.fit(features_train, labels_train)
prediction = clf.predict(features_test)
print("Gausian score",accuracy_score(labels_test, prediction))

# #dumping the mode for future use
# pickle.dump(clf, open( "gausian_model.pkl", "wb" ) )

clf = KNeighborsClassifier()
clf.fit(features_train, labels_train)
prediction = clf.predict(features_test)
print("Nearest neighbour",accuracy_score(labels_test, prediction))

cl = MLPClassifier()
clf.fit(features_train, labels_train)
prediction = clf.predict(features_test)
print("MLP",accuracy_score(labels_test, prediction))

clf = RandomForestClassifier()
clf.fit(features_train, labels_train)
prediction = clf.predict(features_test)
print("RF",accuracy_score(labels_test, prediction))
pickle.dump(clf, open( "randomForest_model.pkl", "wb" ) )

clf = svm.SVC()
clf.fit(features_train, labels_train)
prediction = clf.predict(features_test)
print("SVM",accuracy_score(labels_test, prediction))

clf = MultinomialNB()
clf.fit(features_train, labels_train)
prediction = clf.predict(features_test)
print("Multinomial score",accuracy_score(labels_test, prediction))

clf = BernoulliNB()
clf.fit(features_train, labels_train)
prediction = clf.predict(features_test)
print("Bernoulli score",accuracy_score(labels_test, prediction))
