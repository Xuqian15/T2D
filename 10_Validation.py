#!/usr/bin/env python3

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC,LinearSVC
from sklearn.metrics import roc_curve,auc,recall_score,precision_score,precision_recall_curve,f1_score,accuracy_score,roc_auc_score,matthews_corrcoef
from bayes_opt import BayesianOptimization, UtilityFunction 
import seaborn as sns
from numpy import interp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import argparse

#import data
parser = argparse.ArgumentParser(description = "Internal validation")
parser.add_argument('--Workplace','-W',help = 'Workplace : Input and output work place')
parser.add_argument('--metadata','-m',help = 'input file : metadata')
parser.add_argument('--profile','-p',help = 'input file : microbial profile')
parser.add_argument('--exposure','-e',help = 'input param : the experiment group(exposure) of interest')
parser.add_argument('--group','-g',help = 'input param : the column name of experimental interest(group) in metadata')
parser.add_argument('--batch','-b',help = 'input param : column name of batch(cohort)')
parser.add_argument('--classifier','-c',help = 'input param : selected classifier')
parser.add_argument('--seed','-s',help = 'input param : random seed')
parser.add_argument('--output','-o',help = 'output file prefix: Internal validation result')
args = parser.parse_args()

metadata = pd.read_table(args.Workplace+args.metadata,sep = ',',index_col = 0)
opt_biomarker = pd.read_table(args.Workplace+args.profile,sep = ',',index_col=0)
data_group = np.array([1 if i== str(args.exposure) else 0 for i in metadata[str(args.group)]])
RANDOM_SEED = int(args.seed)
opt_clf = args.classifier


#def function
class machine_learning:   
    def __init__(self):
        self.Method = {'LRl1':LogisticRegression(penalty='l1', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'LRl2':LogisticRegression(penalty='l2', random_state=RANDOM_SEED, solver='liblinear', class_weight='balanced'),
                  'DT':DecisionTreeClassifier(class_weight='balanced', random_state=RANDOM_SEED),
                  'RF':RandomForestClassifier(oob_score=True, class_weight='balanced', random_state=RANDOM_SEED),
                  'GB':GradientBoostingClassifier(random_state=RANDOM_SEED),
                  'KNN':KNeighborsClassifier(n_neighbors=3),
                  'SVC':SVC(class_weight='balanced',random_state=RANDOM_SEED,probability = True)
                  }

    
    def tune_parameter_cv(self,data, data_group,k_fold,**params):
        aucs = []
        tprs = []
        mean_fpr = np.linspace(0, 1, 100)
        i = 0
        splitor = StratifiedKFold(n_splits=k_fold, shuffle=True,random_state=RANDOM_SEED) 
        clf = self.Method[opt_clf].set_params(**params)
        
        for train_index, test_index in splitor.split(data, data_group):
            y_train, y_test = data_group[train_index], data_group[test_index]
            X_train, X_test = np.array(data)[train_index], np.array(data)[test_index]
            
            probas = clf.fit(X_train, y_train).predict_proba(X_test)
            fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            i += 1
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        return mean_auc

    def bayesian_optimise_rf(self,X, y, clf_kfold,k_fold, n_iter = 50, init_points = 5):
        def rf_crossval(n_estimators, max_features,max_depth,max_samples):
            return clf_kfold(
                data = X,
                data_group = y,
                k_fold = k_fold,
                n_estimators = int(n_estimators),
                max_samples = max(min(max_samples,0.999),1e-3),
                max_features = max(min(max_features, 0.999), 1e-3),
                max_depth = int(max_depth),
                bootstrap = True
            )
        
        optimizer = BayesianOptimization(
            random_state = RANDOM_SEED,
            f = rf_crossval,
            pbounds = {
                "n_estimators" : (10, 500),
                "max_features" : (0.1, 0.999),
                "max_samples" : (0.1,0.999),
                "max_depth" : (1,5)
            }
        )
        optimizer.maximize(n_iter = n_iter , init_points = init_points)
        print("Final result:", optimizer.max)
        tune_result = optimizer.max
        tune_result['params']['n_estimators'] = int(tune_result['params']['n_estimators'])
        tune_result['params']['max_depth'] = int(tune_result['params']['max_depth'])
        return tune_result

    def bayesian_optimise_l1(self,X, y, clf_kfold,k_fold, n_iter = 100, init_points = 5):
        def l1_crossval(tol,C):
            return clf_kfold(
                data = X,
                data_group = y,
                k_fold = k_fold,
                tol = max(min(tol,0.1),1e-3),
                C = max(min(C, 0.999), 1e-3)
            )
        
        optimizer = BayesianOptimization(
            random_state = RANDOM_SEED,
            f = l1_crossval,
            pbounds = {
                "tol" : (0.00000001, 0.1),
                "C" : (0,0.999)
            }
        )
        optimizer.maximize(n_iter = n_iter , init_points = init_points)
        print("Final result:", optimizer.max)
        return optimizer.max

    def bayesian_optimise_l2(self,X, y, clf_kfold,k_fold, n_iter = 100, init_points = 5):
        def l2_crossval(tol,C):
            return clf_kfold(
                data = X,
                data_group = y,
                k_fold = k_fold,
                tol = max(min(tol,0.1),1e-3),
                C = max(min(C, 0.999), 1e-3)
            )
        
        optimizer = BayesianOptimization(
            random_state = RANDOM_SEED,
            f = l2_crossval,
            pbounds = {
                "tol" : (0.00000001, 0.1),
                "C" : (0,0.999)
            }
        )
        optimizer.maximize(n_iter = n_iter , init_points = init_points)
        print("Final result:", optimizer.max)
        return optimizer.max

    def bayesian_optimise_dt(self,X, y, clf_kfold,k_fold, n_iter = 50, init_points = 5):
        def dt_crossval(min_samples_leaf,max_depth,min_samples_split):
            return clf_kfold(
                data = X,
                data_group = y,
                k_fold = k_fold,
                min_samples_leaf = max(min(min_samples_leaf, 0.999), 1e-3),
                min_samples_split = max(min(min_samples_split, 0.999), 1e-3),
                max_depth = int(max_depth)
            )
        
        optimizer = BayesianOptimization(
            random_state = RANDOM_SEED,
            f = dt_crossval,
            pbounds = {
                "min_samples_split" : (0.1, 0.999),
                "min_samples_leaf" : (0.00001,0.5),
                "max_depth" : (1,5)
            }
        )
        optimizer.maximize(n_iter = n_iter , init_points = init_points)
        print("Final result:", optimizer.max)
        return optimizer.max

    def bayesian_optimise_gb(self,X, y, clf_kfold,k_fold, n_iter = 50, init_points = 5):
        def gb_crossval(n_estimators,learning_rate,subsample,max_depth,max_features):
            return clf_kfold(
                data = X,
                data_group = y,
                k_fold = k_fold,
                n_estimators = int(n_estimators),
                learning_rate = max(min(learning_rate,0.999),1e-3),
                subsample = max(min(subsample, 0.999), 1e-3),
                max_depth = int(max_depth),
                max_features = max(min(subsample, 0.999), 1e-3)
            )
        
        optimizer = BayesianOptimization(
            random_state = RANDOM_SEED,
            f = gb_crossval,
            pbounds = {
                "n_estimators" : (10, 500),
                "learning_rate" : (0, 1),
                "subsample" : (0.4,1),
                "max_depth" : (1,5),
                "max_features" : (0.1, 0.999)
            }
        )
        optimizer.maximize(n_iter = n_iter , init_points = init_points)
        print("Final result:", optimizer.max)
        tune_result = optimizer.max
        tune_result['params']['n_estimators'] = int(tune_result['params']['n_estimators'])
        tune_result['params']['max_depth'] = int(tune_result['params']['max_depth'])
        return tune_result

    def bayesian_optimise_svc(self,X, y, clf_kfold,k_fold, n_iter = 100, init_points = 5):
        def svc_crossval(C):
            return clf_kfold(
                data = X,
                data_group = y,
                k_fold = k_fold,
                C = max(min(C,0.999),1e-3)
            )
        
        optimizer = BayesianOptimization(
            random_state = RANDOM_SEED,
            f = svc_crossval,
            pbounds = {
                "C" : (0, 0.999)
            }
        )
        optimizer.maximize(n_iter = n_iter , init_points = init_points)
        print("Final result:", optimizer.max)
        return optimizer.max

    def bayesian_optimise_kn(self,X, y, clf_kfold,k_fold, n_iter = 50, init_points = 5):
        def kn_crossval(n_neighbors):
            return clf_kfold(
                data = X,
                data_group = y,
                k_fold = k_fold,
                n_neighbors= int(n_neighbors)
            )
        
        optimizer = BayesianOptimization(
            random_state = RANDOM_SEED,
            f = kn_crossval,
            pbounds = {
                "n_neighbors" : (1,6)
            }
        )
        optimizer.maximize(n_iter = n_iter, init_points = init_points)
        print("Final result:", optimizer.max)
        tune_result = optimizer.max
        tune_result['params']['n_neighbors'] = int(tune_result['params']['n_neighbors'])
        tune_result['params']['max_depth'] = int(tune_result['params']['max_depth'])
        return tune_result

    def model_construction(self,data, data_group, params, k_fold):
        aucs = []
        tprs = []
        mean_fpr = np.linspace(0, 1, 100)
        plot_data = []
        i = 0
        sens = []
        spes = []
        pres = []
        f1s = []
        accus = []
        mccs = []
        auprs = []        
        splitor = StratifiedKFold(n_splits=k_fold, shuffle=True,random_state=RANDOM_SEED) 
        clf = self.Method[opt_clf].set_params(**params)
        
        for train_index, test_index in splitor.split(data, data_group):
            y_train, y_test = data_group[train_index], data_group[test_index]
            X_train, X_test = np.array(data)[train_index], np.array(data)[test_index]
            clf.fit(X_train,y_train)
            
            probas = clf.predict_proba(X_test)
            pred = clf.predict(X_test)
            sen = recall_score(y_test,pred)
            sens.append(sen)
            TP = ((pred==1) & (y_test==1)).sum()
            FP = ((pred==1) & (y_test==0)).sum()
            TN = ((pred==0) & (y_test==0)).sum()
            FN = ((pred==0) & (y_test==1)).sum()
            spe = TN / float(FP + TN)
            spes.append(spe)

            precision, recall, _ = precision_recall_curve(y_test,pred)
            aupr = auc(recall,precision)
            auprs.append(aupr)            
            pre = precision_score(y_test, pred)
            pres.append(pre)
            f1 = f1_score(y_test, pred)
            f1s.append(f1)
            fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
            roc_auc = auc(fpr, tpr)
            spes.append(spe)
            mcc = matthews_corrcoef(y_test,pred)
            mccs.append(mcc)
            aucs.append(roc_auc)
            accu = accuracy_score(y_test, pred)
            accus.append(accu)
            
            
            ### plot data
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            plot_data.append([fpr, tpr, 'ROC Fold %d(AUC = %0.2f)' %(i+1, roc_auc)])
            i += 1
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        mean_spe = np.mean(spes)
        mean_sen = np.mean(sens)
        mean_pre = np.mean(pres)
        mean_f1 = np.mean(f1s)
        mean_accu = np.mean(accus)
        mean_aupr = np.mean(auprs)
        mean_mcc = np.mean(mccs)

        return clf, mean_auc,mean_aupr,mean_mcc,mean_spe,mean_sen,mean_pre,mean_f1,mean_accu,(plot_data, mean_fpr, mean_tpr, tprs,aucs, np.std(aucs))

    def internal_eval(self,X_train,X_test,y_train,y_test,params):
        clf = self.Method[opt_clf].set_params(**params).fit(X_train,y_train)
        pred = clf.predict(X_test)
        probas = clf.predict_proba(X_test)
        fpr, tpr, thresholds = roc_curve(y_test, probas[:, 1])
        roc_auc = auc(fpr, tpr)
        sen = recall_score(y_test,pred)
        precision, recall, _ = precision_recall_curve(y_test,pred)
        aupr = auc(recall,precision)
        mcc = matthews_corrcoef(y_test,pred)
        TP = ((pred==1) & (y_test==1)).sum()
        FP = ((pred==1) & (y_test==0)).sum()
        TN = ((pred==0) & (y_test==0)).sum()
        FN = ((pred==0) & (y_test==1)).sum()
        spe = TN / float(FP + TN)
        pre = precision_score(y_test, pred)
        f1 = f1_score(y_test, pred)
        accu = accuracy_score(y_test, pred)
        

        return roc_auc,aupr,mcc,spe,sen,pre,f1,accu

    def cohort_validation(self,cohorts,metric):
        #cohorts = list(set(metadata[str(args.batch)]))
        evals = pd.DataFrame(columns = cohorts, index = cohorts)
        for i in cohorts:
            for j in cohorts:
                evals.loc[i,j] = self.internal_eval(eval("Cohort"+i),eval("Cohort"+j),eval("Group"+i),eval("Group"+j),eval("best_param"+i))[metric]
        for i in cohorts:
            evals.loc[i,i] = self.model_construction(eval("Cohort"+i), eval("Group"+i), eval("best_param"+i), k_fold=5)[metric+1]

        evals.loc['Average'] = evals.mean()
        evals = pd.concat([evals, pd.Series(name='LOCO').to_frame().T], ignore_index=False)
        #evals = evals.append(pd.Series(name='LOCO'))
        
        for i in cohorts:
            evals.loc["LOCO",i] = self.internal_eval(eval("LOCO"+i),eval("Cohort"+i),eval("LOCO_Group"+i),eval("Group"+i),eval("LOCO_best_param"+i))[metric]
            evals['Average'] = evals.mean(axis=1)
        return evals    

ML = machine_learning()

#Prepare data
cohorts = list(set(metadata[str(args.batch)]))
#cohorts = set(metadata[str(args.batch)])
datasets = locals()
for i in cohorts:
    datasets["Cohort"+str(i)] = opt_biomarker.loc[metadata[metadata[str(args.batch)]==i].index,]
    datasets["Group"+str(i)] = np.array([1 if i== str(args.exposure) else 0 for i in metadata[metadata[str(args.batch)]==i][str(args.group)]])
    datasets["LOCO"+str(i)] = opt_biomarker.loc[metadata[metadata[str(args.batch)]!=i].index,]
    datasets['LOCO_Group'+str(i)] = np.array([1 if i== str(args.exposure) else 0 for i in metadata[metadata[str(args.batch)]!=i][str(args.group)]])

#Prepare hyperparameters
dict = {'LRl1':ML.bayesian_optimise_l1,
'LRl2':ML.bayesian_optimise_l2,
'DT':ML.bayesian_optimise_dt,
'RF':ML.bayesian_optimise_rf,
'GB':ML.bayesian_optimise_gb,
'KNN':ML.bayesian_optimise_kn,
'SVC':ML.bayesian_optimise_svc}

params = locals()
for i in cohorts:
    tune_result = dict[opt_clf](eval("Cohort"+i),eval("Group"+i),ML.tune_parameter_cv,k_fold=5)
    params["best_param"+str(i)]= tune_result['params']
    
    tune_result = dict[opt_clf](eval("LOCO"+i),eval("LOCO_Group"+i),ML.tune_parameter_cv,k_fold=5)
    params["LOCO_best_param"+str(i)]= tune_result['params']

#Cohort-to-cohort & LOCO

def validation_plot(plot_heatmap, metric):
    max_plot = np.array(plot_heatmap.max()).max()
    min_plot = np.array(plot_heatmap.min()).min()

    sns.set(font_scale=1.5)
    grid_kws = {"height_ratios":(lc,1,1),"width_ratios":(lc,1)}
    f, axs= plt.subplots(3,2, gridspec_kw=grid_kws)
    f.suptitle("Internal validation "+metric, fontsize=20, fontweight='bold', y=1.05)

    sns.heatmap(plot_heatmap.iloc[0:lc,0:lc],cmap="YlGnBu", annot = True, ax = axs[0,0],vmin=min_plot,vmax=max_plot,cbar = False,fmt='.2')
    axs[0,0].xaxis.set_ticks_position('top')
    sns.heatmap(pd.DataFrame(plot_heatmap.loc['LOCO',cohorts]).T,cmap="YlGnBu", annot = True, ax = axs[2,0],vmin=min_plot,vmax=max_plot, cbar = False,xticklabels=False,fmt='.2')
    axs[2,0].set_yticklabels(["LOCO"],rotation='horizontal')
    sns.heatmap(pd.DataFrame(plot_heatmap.iloc[lc,0:lc]).T,cmap="YlGnBu", annot = True, ax = axs[1,0],vmin=min_plot,vmax=max_plot,cbar = False,xticklabels=False,fmt='.2')
    axs[1,0].set_yticklabels(["Average"],rotation='horizontal')
    sns.heatmap(pd.DataFrame(plot_heatmap.iloc[0:lc,lc]),cmap="YlGnBu", annot = True, ax = axs[0,1],vmin=min_plot,vmax=max_plot,cbar = False,yticklabels=False,fmt='.2')
    axs[0,1].xaxis.set_ticks_position('top')
    sns.heatmap(pd.DataFrame(plot_heatmap.loc['Average','Average'],index = ["Average"],columns = ["Average"]),cmap="YlGnBu", annot = True, ax = axs[1,1],vmin=min_plot,vmax=max_plot,cbar = False,xticklabels=False,yticklabels=False,fmt='.2')
    sns.heatmap(pd.DataFrame(plot_heatmap.loc['LOCO','Average'],index = ["LOCO"],columns = ["Average"]),cmap="YlGnBu", annot = True, ax = axs[2,1],vmin=min_plot,vmax=max_plot,cbar = False,xticklabels = False,yticklabels=False,fmt='.2')
    plt.subplots_adjust(wspace =.05, hspace =.1)

    return plt

metrics = ['AUC','AUPR','MCC','Specificity','Sensitivity','Precision','F1','Accuracy']
lc = len(cohorts)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
for i in range(8):
    evals = ML.cohort_validation(cohorts,i)
    evals.to_csv(args.Workplace+args.output+"_"+args.classifier+'_validation_'+metrics[i]+'.txt',sep = '\t')
    
    plot_heatmap = pd.DataFrame(evals, dtype='float')
    fig = validation_plot(plot_heatmap,metrics[i])
    #fig.title("Internal validation", fontsize=20, fontweight='bold', pad=20)
    fig.savefig(args.Workplace+args.output+"_"+args.classifier+'_validation_'+metrics[i]+'.pdf',bbox_inches = 'tight')
    fig.savefig(args.Workplace+args.output+"_"+args.classifier+'_validation_'+metrics[i]+'.svg',bbox_inches = 'tight',format = 'svg')
print("FINISH")

