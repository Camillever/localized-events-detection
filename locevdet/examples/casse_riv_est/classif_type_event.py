""" Script to classify the type of event - RandomForest """

import os
import numpy as np
import pickle


################################### Create X, Y train and test ###################################
## X (for train) Already labelled
rockfall_eventlist_path = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', 'rockfall', \
    'eventlist')
rockfalls_eventlist = "eventlist_2020-02-01T00-00-00_2020-02-11T00-00-00.pkl"

start_period = rockfalls_eventlist.split('_')[1]
end_period = (rockfalls_eventlist.split('_')[2]).split('.')[0]
print(f"start and end period : {start_period} -{end_period}")
eventlist_filepath = os.path.join(rockfall_eventlist_path, rockfalls_eventlist)
with open(eventlist_filepath, "rb") as content:
    eventlist = pickle.load(content)

    X_train = eventlist.to_array()
    print('X_train :', X_train)
    print("X_train shape :", np.shape(X_train))

# X Unknown (for test)
continuous_eventlist_path = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
    'continuous', 'eventlist')
all_eventlists = os.listdir(os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
            'continuous', 'eventlist'))
X_test = None
for eventlist_file in all_eventlists:
    start_period = eventlist_file.split('_')[2]
    end_period = (eventlist_file.split('_')[3]).split('.')[0]
    print(f"start and end period : {start_period} -{end_period}")
    eventlist_filepath = os.path.join(continuous_eventlist_path, eventlist_file)
    with open(eventlist_filepath, "rb") as content:
        eventlist = pickle.load(content)

        X_test_per_eventlist = eventlist.to_array()

        if X_test is None:
            X_test = X_test_per_eventlist
        else:
            X_test = np.concatenate((X_test, X_test_per_eventlist), axis=0)
        print("X_test shape (in progress):", np.shape(X_test))
        
print("X_test :", X_test)
print("X_test shape :", np.shape(X_test))


################################ BINARY CLASSIFICATION ##########################################
################################### by Random Forest ############################################


 


