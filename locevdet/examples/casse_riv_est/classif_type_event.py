""" Script to classify the type of event - RandomForest """

import os
import pickle
import numpy as np


################################### Create X, Y train and test ###################################
## X (for train) Already labelled
rockfall_eventlist_path = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', 'rockfall', \
    'eventlist')
rockfalls_eventlist = "eventlist_2020-02-01T00-00-00_2020-02-11T00-00-00.pkl"

volcano_eventlist_path = os.path.join('dictionaries_data', '01-02-2020_11-02-2020', \
    'volcanotectonic', 'eventlist')

start_period = rockfalls_eventlist.split('_')[1]
end_period = (rockfalls_eventlist.split('_')[2]).split('.')[0]
print(f"start and end period : {start_period} -{end_period}")
eventlist_filepath = os.path.join(rockfall_eventlist_path, rockfalls_eventlist)
with open(eventlist_filepath, "rb") as content:
    eventlist = pickle.load(content)

    X_train = eventlist.to_array()
    print('X_train :', X_train)
    print("X_train shape :", np.shape(X_train))
    events_rows = X_train[:,:1]

    rockfall_id = np.ones((events_rows.shape[0], 1))
    Y_train = np.concatenate((events_rows, rockfall_id), axis=1)
    # print('Y_train :', Y_train)
    print("Y_train shape :", np.shape(Y_train))


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

events_rows = X_train[:,:1]

rockfall_id = np.zeros((events_rows.shape[0], 1))
Y_test = np.concatenate((events_rows, rockfall_id), axis=1)
# print('Y_train :', Y_train)
print("Y_test shape :", np.shape(Y_test))







# We saved the numpy.ndarray to the .npz format
# This will allow us to load the dataset faster later on
np.savez(
    'dataset_split',
    X_train=X_train,
    Y_train=Y_train,
    X_test=X_test,
    Y_test=Y_test
)
print(f"Dataset sauvegardé dans {os.getcwd()}")
dataset_split = np.load('dataset_split.npz')
print('Dataset :',dataset_split.files)
print('Exemple - taille de X_train :',dataset_split['X_train'].shape)

################################ BINARY CLASSIFICATION ##########################################
################################### by Random Forest ############################################
def accuracy(y_pred, y_true):
    return np.mean(y_pred == y_true)

def get_confusion_data(y_pred, y_true):
    y_true_is_positive = y_true == 1
    true_positives = np.mean(y_pred[y_true_is_positive] == y_true[y_true_is_positive])
    false_negatives = np.mean(y_pred[y_true_is_positive] != y_true[y_true_is_positive])

    y_true_is_negative = np.logical_not(y_true_is_positive)
    true_negatives = np.mean(y_pred[y_true_is_negative] == y_true[y_true_is_negative])
    false_positives = np.mean(y_pred[y_true_is_negative] != y_true[y_true_is_negative])

    return true_positives, false_negatives, true_negatives, false_positives

def recall(y_pred, y_true):
    TP, FN, TN, FP = get_confusion_data(y_pred, y_true)
    return TP / (TP + FN)

def precision(y_pred, y_true):
    TP, FN, TN, FP = get_confusion_data(y_pred, y_true)
    return TP / (TP + FP)

def f1_score(y_pred, y_true):
    TP, FN, TN, FP = get_confusion_data(y_pred, y_true)
    return TP / (TP + (FP + FN)/2)



###########################################

# from sklearn.ensemble import RandomForestClassifier

# # Load the data 
# dataset = np.load('dataset_split.npz')

# X_train = dataset['X_train']
# Y_train = dataset['Y_train']
# X_test = dataset['X_test']
# Y_test = dataset['Y_test']

# import wandb

# sweep_config = {
#   "name": "My Sweep",
#   "method": "random",
#   "metric": "f1_score",
#   "goal": "maximize",
#   "parameters": {
#         "model": {
#             "value": "RandomForestClassifier"
#         },
#         "n_estimators": {
#             "distribution": "q_log_uniform",
#             "min": 0,
#             "max": 5,
#             "q": 1,
#         },
#         "max_depth": {
#             "distribution": "int_uniform",
#             "min": 2,
#             "max": 15,
#         },
#         "rockfall_weight": {
#             "distribution": "q_log_uniform",
#             "min": 0,
#             "max": 7.6,
#             "q": 1,
#         }
#     }
# }

# sweep_id = wandb.sweep(sweep_config, project="rockfalls")

# import wandb
# import time

# def train(n_fold=5, metrics=('accuracy', 'recall', 'precision', 'f1_score')):
#     run = wandb.init(project="rockfalls")
#     config = dict(run.config)

#     rockfall_weight = config.pop('rockfall_weight', 1)
#     class_weight = {0:1, 1:rockfall_weight}
#     print("Config:", config)

#     all_index = np.arange(len(Y_train_pixels))
#     indexes_permuted = np.random.permutation(all_index)
#     val_indexes_list = np.array_split(indexes_permuted, n_fold)

#     train_metrics = {metric_name: [] for metric_name in metrics}
#     val_metrics = {metric_name: [] for metric_name in metrics}

#     model_class = eval(config.pop('model'))

#     for i in range(n_fold):
#         print(f"Fold {i}:")
#         #Separate train et val
#         val_index = val_indexes_list[i]
#         train_index = np.delete(indexes_permuted, val_index)

#         x_train, y_train = X_train_pixels[train_index], Y_train_pixels[train_index]
#         x_val, y_val = X_train_pixels[val_index], Y_train_pixels[val_index]

#         model = model_class(class_weight=class_weight, **config)

#         t0 = time.time()
#         model.fit(x_train, y_train)
#         training_time = time.time() - t0
        
#         y_pred_train = model.predict(x_train)
#         y_pred_val = model.predict(x_val)
#         for metric_name in metrics:
#             metric = eval(metric_name)
#             train_value = metric(y_pred_train, y_train)
#             val_value = metric(y_pred_val, y_val)
#             print(f"\t{metric_name} - train:{train_value:.4f} val:{val_value:.4f}")
#             train_metrics[metric_name].append(train_value)
#             val_metrics[metric_name].append(val_value)
        
#         for f1_val in val_metrics['f1_score']:
#             if f1_val < 0.9:
#                 return

#     print(metrics)
#     for metric_name in metrics:
#         mean_value_train = np.mean(train_metrics[metric_name])
#         mean_value_val = np.mean(val_metrics[metric_name])
    
#         std_value_train = np.std(train_metrics[metric_name])
#         std_value_val = np.std(val_metrics[metric_name])

#         print(f"{metric_name}:  \t {mean_value_train:.4f}±{std_value_train:.4f}")
#         print(f"{'val_' + metric_name}:  \t {mean_value_val:.4f}±{std_value_val:.4f}")

#         wandb.log({metric_name: mean_value_train}, commit=False)
#         wandb.log({'val_' + metric_name: mean_value_val}, commit=False)

#     wandb.log({'training_time': training_time})
#     wandb.sklearn.plot_feature_importances(model, features_names)

# wandb.agent(sweep_id, function=train)


 


