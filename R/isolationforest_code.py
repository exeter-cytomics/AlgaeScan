#####################  import libraries #####

# import package
import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
import sys

##################### define functions #####


# function to perform training with isoforest model

def iso_forest_train(Xtrain,cont_val=0.1,random_state=42):
  random_state=int(random_state)
  # Initialize the Isolation Forest model
  model = IsolationForest(n_estimators=100, contamination=cont_val, random_state=random_state)

  # Fit the model to the data
  model_isoforest=model.fit(Xtrain)

  return model_isoforest

# function to perform prediction using trained model
def iso_forest_predict(isof_model,test_df):
  y_pred = isof_model.predict(test_df)
  y_pred = y_pred.flatten()
  return y_pred

# Probability-like anomaly scores
def iso_forest_predict_prob(isof_model, test_df):
    # decision_function gives scores (higher = more normal)
    scores = isof_model.decision_function(test_df)
    # Normalize scores to [0,1] range as "probabilities"
    probs = (scores - scores.min()) / (scores.max() - scores.min())
    return probs


# Note: if you are defininig a function inside a python script within R studio, make sure you highline also one line below the end of the function when you press the run command.
