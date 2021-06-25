#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to run classifiers to evaluate the goodness of feature extraction
Created on Wed Mar  3 21:39:00 2021

@author: givanna
"""

import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score, make_scorer
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC, LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score



def normalise_to_0_1(train_dat, test_dat):
    """
    Transform data to range between 0 and 1.
    Transformer is trained on training data, and applied on test data

    Parameters
    ----------
    train_dat : pandas.DataFrame
        Training data (features only).
    test_dat : pandas.DataFrame
        Test data (features only)

    Returns
    -------
    train_dat_copy : pandas.DataFrame.
        Training data normalised.
    test_dat_copy : pandas.DataFrame.
        Test data normalised

    """
    
    # copy data frame to avoid changing parameters
    train_dat_copy = train_dat.copy()
    test_dat_copy = test_dat.copy()
    
    scaler = MinMaxScaler()
    scaler.fit(train_dat_copy)
    
    # have to store the column names
    col_names = train_dat_copy.columns
    train_dat_copy = pd.DataFrame(scaler.transform(train_dat_copy), columns=col_names)
    col_names = test_dat_copy.columns
    test_dat_copy = pd.DataFrame(scaler.transform(test_dat_copy), columns=col_names)
    
    return train_dat_copy, test_dat_copy


def run_knn(train_dat, train_label, test_dat, test_label, 
            k_vals = range(1, 11), num_folds=None):
    """
    
    Train KNN classifier using Grid Search Cross Validation, and test
    the `optimal` tuned classifier on the test data.

    Parameters
    ----------
    train_dat : pandas.DataFrame
        Training data (features only) used to train KNN.
    train_label : array
        Label for each training data (1 label per row of train_dat).
    test_dat : pandas.DataFrame
        Test data (features only) to run the trained KNN on.
    test_label : array
        Label for each test data (1 label per row of test_dat).
    k_vals : array, optional
        Array of integer containing the list of K-values to train. 
        The default is range(1, 11).
    num_folds : int, optional
        Number of folds to be created for Grid Search CV.
        The default is None.
        If None, it will be determined based on the smallest training label.

    Returns
    -------
    res : dict
        Dictionary containing:
            1. best_param: best parameter as determined by Grid Search CV.
            2. best_param_train_score: best parameter's performance as evaluated
            using cohen kappa metric.
            3. predicted_test_label: predicted label of the test data produced
            by the `optimal` KNN classifier.
            4. proportion_correct_test_label: proportion of test label correctly
            predicted.
            5. num_folds: number of folds used in performing CV. Unless specified,
            this is automatically computed based on the smallest training label.

    """
    
    print("Running K-Nearest Neighbours")
    
    # Find the optimal parameter using grid search on training data.
    # Report the optimal Cohen Kappa metric (because we have imbalance classes).
    # Run the optimal KNN on the test_dat.
    # Cross check the KNN result with test_label, and return proportion 
    # of patients classified correctly.
    knn = KNeighborsClassifier(metric='euclidean')
        
    grid_search = GridSearchCV(
        estimator=knn,
        param_grid={"n_neighbors": k_vals},
        # scoring = make_scorer(cohen_kappa_score),
        scoring='accuracy',
        n_jobs=16,
        cv=LeaveOneOut(),
        return_train_score=True,
        verbose=0
        )
    grid_search.fit(X=train_dat, y=train_label)
    
    # Test the KNN on test data
    predicted_label = grid_search.predict(X=test_dat)
    cnt_correct = np.count_nonzero(predicted_label == np.array(test_label))
    prop_correct = cnt_correct / test_label.size
    
    res = {
        "best_param": grid_search.best_params_,
        "best_param_train_score": grid_search.best_score_,
        "predicted_test_label": predicted_label,
        "proportion_correct_test_label": round(prop_correct, 3)
        }
    
    return(res)
        


def run_nb(train_dat, train_label, test_dat, test_label):
    """
    
    Train Naive Bayes classifier on train_dat, and test it on test_dat.

    Parameters
    ----------
    train_dat : pandas.DataFrame
        Training data (features only) used to train KNN.
    train_label : array
        Label for each training data (1 label per row of train_dat).
    test_dat : pandas.DataFrame
        Test data (features only) to run the trained KNN on.
    test_label : array
        Label for each test data (1 label per row of test_dat).
    

    Returns
    -------
    res : dict
        Dictionary containing:
            1. predicted_test_label: predicted label of the test data produced
            by the `optimal` Naive Bayes classifier.
            2. proportion_correct_test_label: proportion of test label correctly
            predicted.
            
    """
    
    print("Running Naive Bayes")
    
    # NB has no parameter, so just run normal loocv for the training data
    cv = LeaveOneOut()
    
    lab_true, lab_pred = list(), list()
    for train_ix, test_ix in cv.split(train_dat):
    	# split data
    	dat_train, dat_test = train_dat.to_numpy()[train_ix, :], train_dat.to_numpy()[test_ix, :]
    	lab_train, lab_test = train_label[train_ix], train_label[test_ix]
    	# fit model
    	model = GaussianNB()
    	model.fit(dat_train, lab_train)
    	# evaluate model
    	pred = model.predict(dat_test)
    	# store
    	lab_true.append(lab_test[0])
    	lab_pred.append(pred[0])
    
    train_acc = accuracy_score(lab_true, lab_pred)
    
    
    gnb = GaussianNB()
    
    gnb.fit(X=train_dat, y=train_label)
    
    predicted_label = gnb.predict(X=test_dat)
    cnt_correct = np.count_nonzero(predicted_label == np.array(test_label))
    prop_correct = cnt_correct / test_label.size
    
    res = {
        "predicted_test_label": predicted_label,
        "proportion_correct_test_label": round(prop_correct, 3),
        "best_param_train_score": train_acc
        }
    
    return(res)


def run_dt(train_dat, train_label, test_dat, test_label, 
           max_depth_params = [None, 2, 4, 8, 16, 32, 64, 128, 256, 512],
           min_samples_split_params = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
           min_samples_leaf_params = [0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50], 
           num_folds=None):
    """
    
    Train Decision Tree classifier using Grid Search Cross Validation, and test
    the `optimal` tuned classifier on the test data.
    
    By default, the min_samples_leaf_params and min_samples_split params are 
    2 fold increase from default value of 1 and 2 respectively, keeping
    the number of elements to 10.

    Parameters
    ----------
    train_dat : pandas.DataFrame
        Training data (features only) used to train KNN.
    train_label : array
        Label for each training data (1 label per row of train_dat).
    test_dat : pandas.DataFrame
        Test data (features only) to run the trained KNN on.
    test_label : array
        Label for each test data (1 label per row of test_dat).
    min_samples_split_params : array, optional
        Array of integer or float containing the list of min_samples_split to train. 
        The default is [1, 2, 4, 8, 16, 32, 64, 128, 256, 512].
    min_samples_leaf_params : array, optional
        Array of integer containing the list of min_samples_leaf to train. 
        The default is [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024].
    num_folds : int, optional
        Number of folds to be created for Grid Search CV.
        The default is None.
        If None, it will be determined based on the smallest training label.

    Returns
    -------
    res : dict
        Dictionary containing:
            1. best_param: best parameter as determined by Grid Search CV.
            2. best_param_train_score: best parameter's performance as evaluated
            using cohen kappa metric.
            3. predicted_test_label: predicted label of the test data produced
            by the `optimal` Decision Tree classifier.
            4. proportion_correct_test_label: proportion of test label correctly
            predicted.
            5. num_folds: number of folds used in performing CV. Unless specified,
            this is automatically computed based on the smallest training label.

    """
    
    print("Running Decision Tree")
    
    # Find the optimal parameter using grid search on training data.
    # Report the optimal Cohen Kappa metric (because we have imbalance classes).
    # Run the optimal DT on the test_dat.
    # Cross check the DT result with test_label, and return proportion 
    # of patients classified correctly.
    # Random state is the seed.
    dt = DecisionTreeClassifier(random_state=42)
        
    grid_search = GridSearchCV(
        estimator=dt,
        param_grid={
            "min_samples_split": min_samples_split_params,
            "min_samples_leaf": min_samples_leaf_params,
            "max_depth": max_depth_params
            },
        # scoring = make_scorer(cohen_kappa_score),
        scoring='accuracy',
        n_jobs=16,
        cv=LeaveOneOut(),
        return_train_score=True,
        verbose=0
        )
    grid_search.fit(X=train_dat, y=train_label)
    
    # Test on test data
    predicted_label = grid_search.predict(X=test_dat)
    cnt_correct = np.count_nonzero(predicted_label == np.array(test_label))
    prop_correct = cnt_correct / test_label.size
    
    res = {
        "best_param": grid_search.best_params_,
        "best_param_train_score": grid_search.best_score_,
        "predicted_test_label": predicted_label,
        "proportion_correct_test_label": round(prop_correct, 3)
        }
    
    return(res)


def run_svm(train_dat, train_label, test_dat, test_label, 
            c_params = [2e-5, 2e-3, 2e-1, 2e1, 2e3, 2e5, 2e7, 2e9, 2e11, 2e13],
            gamma_params = [2e-15, 2e-13, 2e-11, 2e-9, 2e-7, 2e-5, 2e-3, 2e-1, 2e1, 2e3],
            max_iters = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000],
            num_folds = None):
    """
    
    Train Support Vector Machine classifier using Grid Search Cross Validation, and test
    the `optimal` tuned classifier on the test data.
    Grid search will use linear and rbf kernel only. 
    For Rbf kernel, gamma_params will be investigated.
    Same C values specified in c_params will be used for both Rbf and linear kernel.
    Default Gamma and C params are obtained from
    https://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf.
    However, we kept the number of element to 10, consistent with others.

    Parameters
    ----------
    train_dat : pandas.DataFrame
        Training data (features only) used to train KNN.
    train_label : array
        Label for each training data (1 label per row of train_dat).
    test_dat : pandas.DataFrame
        Test data (features only) to run the trained KNN on.
    test_label : array
        Label for each test data (1 label per row of test_dat).
    c_params : array, optional
        Array of float containing the list of parameter `C` to train. 
        The default is [2e-5, 2e-3, 2e-1, 2e1, 2e3, 2e5, 2e7, 2e9, 2e11, 2e13].
    gamma_params : array, optional
        Array of integer containing the list of parameter `gamma` to train. 
        The default is [2e-15, 2e-13, 2e-11, 2e-9, 2e-7, 2e-5, 2e-3, 2e-1, 2e1, 2e3].
    num_folds : int, optional
        Number of folds to be created for Grid Search CV.
        The default is None.
        If None, it will be determined based on the smallest training label.

    Returns
    -------
    res : dict
        Dictionary containing:
            1. best_param: best parameter as determined by Grid Search CV.
            2. best_param_train_score: best parameter's performance as evaluated
            using cohen kappa metric.
            3. predicted_test_label: predicted label of the test data produced
            by the `optimal` SMV classifier.
            4. proportion_correct_test_label: proportion of test label correctly
            predicted.
            5. num_folds: number of folds used in performing CV. Unless specified,
            this is automatically computed based on the smallest training label.

    """
    
    print("Running Support Vector Machine")
    
    # Find the optimal parameter using grid search on training data.
    # Report the optimal Cohen Kappa metric (because we have imbalance classes).
    # Run the optimal SVM on the test_dat.
    # Cross check the SVM result with test_label, and return proportion 
    # of patients classified correctly.
    # Random state is the seed.
    svm = SVC(random_state=42)
    
    param_grid=[
        {'kernel': ['rbf'], 'gamma': gamma_params, 'C': c_params, 'max_iter': max_iters},
        {'kernel': ['linear'], 'C': c_params, 'max_iter': max_iters}
        ]
        
    grid_search = GridSearchCV(
        estimator=svm,
        param_grid=param_grid,
        # scoring = make_scorer(cohen_kappa_score),
        scoring='accuracy',
        n_jobs=16,
        cv=LeaveOneOut(),
        return_train_score=True,
        verbose=0
        )
    grid_search.fit(X=train_dat, y=train_label)
    
    # Test on test data
    predicted_label = grid_search.predict(X=test_dat)
    cnt_correct = np.count_nonzero(predicted_label == np.array(test_label))
    prop_correct = cnt_correct / test_label.size
    
    res = {
        "best_param": grid_search.best_params_,
        "best_param_train_score": grid_search.best_score_,
        "predicted_test_label": predicted_label,
        "proportion_correct_test_label": round(prop_correct, 3)
        }
    
    return(res)


def run_rf(train_dat, train_label, test_dat, test_label, 
           bootstrap_params = [True, False],
           # max_depth_params = [None, 2, 4, 8, 16, 32, 64, 128, 256, 512],
           # min_samples_split_params = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
           # min_samples_leaf_params = [0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50], 
           n_estimators_params = [3, 6, 12, 25, 50, 100, 200, 400, 800, 1600],
           num_folds=None):
    """
    
    Train Random Forest classifier using Grid Search Cross Validation, and test
    the `optimal` tuned classifier on the test data.
    
    On the default parameter values varied by Grid search, 
    vary everything by 2 fold-change from default value, and keep the size to 10.
    max_depth default is None. We just start with 2.
    n_estimators (number of trees) default is 100, we just do 2 fold increase
    from 100
    Same for min_samples_split and min_samples_leaf. 2 fold-change from default,
    but only increase, not decrease as the default values are the minimum.
    bootstrap is kept to True and False

    Parameters
    ----------
    train_dat : pandas.DataFrame
        Training data (features only) used to train KNN.
    train_label : array
        Label for each training data (1 label per row of train_dat).
    test_dat : pandas.DataFrame
        Test data (features only) to run the trained KNN on.
    test_label : array
        Label for each test data (1 label per row of test_dat).
    bootstrap_params : array, optional
        Whether bottstrap samples are used when building trees. 
        See sklearn.ensemble.RandomForestClassifier.
        The default is [True, False].
    max_depth_params : array, optional
        Maximum depth of tree. 
        See sklearn.ensemble.RandomForestClassifier.
        The default is [None, 2, 4, 8, 16, 32, 64, 128, 256],.
    n_estimators_params : TYPE, optional
        Number of trees in the forest. 
        See sklearn.ensemble.RandomForestClassifier.
        The default is [25, 50, 100, 200, 400].
    num_folds : int, optional
        Number of folds to be created for Grid Search CV.
        The default is None.
        If None, it will be determined based on the smallest training label.

    Returns
    -------
    res : dict
        Dictionary containing:
            1. best_param: best parameter as determined by Grid Search CV.
            2. best_param_train_score: best parameter's performance as evaluated
            using cohen kappa metric.
            3. predicted_test_label: predicted label of the test data produced
            by the `optimal` Random Forest classifier.
            4. proportion_correct_test_label: proportion of test label correctly
            predicted.
            5. num_folds: number of folds used in performing CV. Unless specified,
            this is automatically computed based on the smallest training label.

    """
    
    print("Running Random Forest")
    
    # Find the optimal parameter using grid search on training data.
    # Report the optimal Cohen Kappa metric (because we have imbalance classes).
    # Run the optimal RF on the test_dat.
    # Cross check the RF result with test_label, and return proportion 
    # of patients classified correctly.
    # Random state is the seed.
    rf = RandomForestClassifier(random_state=42)
    
    param_grid = {
        'bootstrap': bootstrap_params,
        # 'max_depth': max_depth_params,
        # 'min_samples_leaf': min_samples_leaf_params,
        # 'min_samples_split': min_samples_split_params,
        'n_estimators': n_estimators_params
        }
        
    grid_search = GridSearchCV(
        estimator=rf,
        param_grid=param_grid,
        # scoring = make_scorer(cohen_kappa_score),
        scoring='accuracy',
        n_jobs=16,
        cv=LeaveOneOut(),
        return_train_score=True,
        verbose=0
        )
    grid_search.fit(X=train_dat, y=train_label)
    
    # Test on test data
    predicted_label = grid_search.predict(X=test_dat)
    cnt_correct = np.count_nonzero(predicted_label == np.array(test_label))
    prop_correct = cnt_correct / test_label.size
    
    res = {
        "best_param": grid_search.best_params_,
        "best_param_train_score": grid_search.best_score_,
        "predicted_test_label": predicted_label,
        "proportion_correct_test_label": round(prop_correct, 3)
        }
    
    return(res)


def run_logistic_regression(train_dat, train_label, test_dat, test_label,
                            Cs = 10,
                            solver_params = ['lbfgs', 'sag', 'newton-cg', 'liblinear'],
                            tol_params = np.logspace(-5, 5, 10, base=2),
                            max_iter_params = np.logspace(0, 9, 10, base=2) * 100,
                            num_folds = None):
    """
    
    Train Logistic Regression classifier using Grid Search Cross Validation, and test
    the `optimal` tuned classifier on the test data.
    
    By default, C parameters varied are a grid in logartihmic scale between 1e-4
    and 1e4 of size Cs.
    Refer to sklearn.linear_model.LogisticRegressionCV documentation.
    
    I2 penalty is used as it's the default in sklearn.
    
    solver_params are set defaulted for small datasets and multiclass prolems.

    tol_params (increase and decrease) is 2 fold 
    change from default, kept to maximum array size of 5.
    
    Cannot keep max_iter as 2 fold change from default as it failed to converge.
    That was the original [100, 200, 400, 800, 1600]
    Thus we start higher than 2000, and linearly increase it by 2 folds
    
    Parameters
    ----------
    train_dat : pandas.DataFrame
        Training data (features only) used to train KNN.
    train_label : array
        Label for each training data (1 label per row of train_dat).
    test_dat : pandas.DataFrame
        Test data (features only) to run the trained KNN on.
    test_label : array
        Label for each test data (1 label per row of test_dat).
    Cs : int
        The number of parameter C to investigate. 
        The default is np.logspace(-4, 4, Cs) to generate Cs evenly spread 
        number of logarithmic scale.
    tol_params : array
        The tolerance for stopping criteria.
        The default is np.logspace(-5, 5, 10, base=2).
    max_iter_params : array
        Maximum number of iterations of the optimisation algorithm.
        If method returns warning not able to converge, increase the start
        (and maybe end) of this array.
        The default is np.logspace(0, 9, 10, base=2) * 100.
    num_folds : int, optional
        Number of folds to be created for Grid Search CV.
        The default is None.
        If None, it will be determined based on the smallest training label.

    Returns
    -------
    res : dict
        Dictionary containing:
            1. best_param: best parameter as determined by Grid Search CV.
            2. best_param_train_score: best parameter's performance as evaluated
            using cohen kappa metric.
            3. predicted_test_label: predicted label of the test data produced
            by the `optimal` SMV classifier.
            4. proportion_correct_test_label: proportion of test label correctly
            predicted.
            5. num_folds: number of folds used in performing CV. Unless specified,
            this is automatically computed based on the smallest training label.

    """
    
    print("Running Logistic Regression")
    
    # Find the optimal parameter using grid search on training data.
    # Report the optimal Cohen Kappa metric (because we have imbalance classes).
    # Run the optimal Logistic Regression on the test_dat.
    # Cross check the Logistic Regression result with test_label, and return proportion 
    # of patients classified correctly.
    # Random state is the seed.
    log_reg = LogisticRegression(random_state=42)
    
    param_grid = {
        'C': np.logspace(-4, 4, num=Cs),
        'solver': solver_params,
        'max_iter': max_iter_params,
        'tol': tol_params
        }
    
    grid_search = GridSearchCV(
        estimator=log_reg,
        param_grid=param_grid,
        # scoring = make_scorer(cohen_kappa_score),
        scoring='accuracy',
        n_jobs=16,
        cv=LeaveOneOut(),
        return_train_score=True,
        verbose=0
        )
    grid_search.fit(X=train_dat, y=train_label)
    
    # Test on test data
    predicted_label = grid_search.predict(X=test_dat)
    cnt_correct = np.count_nonzero(predicted_label == np.array(test_label))
    prop_correct = cnt_correct / test_label.size
    
    res = {
        "best_param": grid_search.best_params_,
        "best_param_train_score": grid_search.best_score_,
        "predicted_test_label": predicted_label,
        "proportion_correct_test_label": round(prop_correct, 3)
        }
    
    return(res)


