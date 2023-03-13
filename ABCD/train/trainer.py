from ABCD.train.utils_tsai import *
from ABCD.train.cam import *
import torch.nn.functional as F
import torch.nn as nn
import torch
import os
from tsai.all import *
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

my_setup()

class ABCDTrainer:
    def __init__(self, config):
        self.pretrain = config['pretrain']
        self.freeze_epochs = config['freeze_epoches']
        self.metrics = ["Precision", "Recall", "accuracy", "RocAucBinary", "F1Score", "BalancedAccuracy"]
        # self.metrics = ["RocAucBinary"]
        self.subject_path = config['subject_path']
        self._result = _result = [[] for i in range(len(self.metrics))]
        self.model = config['model']
        self.lr = config['lr']
        self.loss_fn = LabelSmoothingCrossEntropyFlat()
        self.save_path = config['save_path']
        
        if config['add_genomics'] and config['add_c4']:
            self.feature_categories = ['ts', 'covs', 'cog_task', 'family_his', 'score', 'genomic', 'c4_feats', 'ts_ind']
        elif config['add_c4']:
            self.feature_categories = ['ts', 'covs', 'cog_task', 'family_his', 'score', 'c4_feats', 'ts_ind']
        else:
            self.feature_categories = ['ts', 'covs', 'cog_task', 'family_his', 'score', 'ts_ind']
        self.feature_names = get_feature_names(mapped_from=self.feature_categories)
        self.class_balance = config['class_balance']
        self.config = config
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)


    def get_model(self, dls):
        model = self.config['model']
        if model in ['InceptionTimePlus', 'XceptionTimePlus']:
            return eval(model)
        elif model == 'TST':
            return  nn.Sequential(nn.Conv1d(dls.vars, 128, kernel_size=(7,), stride=(1,), padding=(3,), bias=False), 
                                nn.BatchNorm1d(128, eps=1e-05, momentum=0.1, affine=True, track_running_stats=True),
                                nn.ReLU(), 
                                nn.AdaptiveAvgPool1d(output_size= 512), 
                                TSTPlus(128, dls.c, seq_len = 512, max_seq_len = 512, dropout = .1, store_attn = True))
        elif model == 'RNN':
            return RNN(dls.vars, dls.c)
        elif model == 'GRU_FCN':
            return GRU_FCN(dls.vars, dls.c, seq_len = 2880)
        elif model == 'MiniRocket':
            return  MiniRocket(dls.vars, dls.c, 2880)
        else:
            raise NotImplementedError('Model Not Implemented in this package yet')
    
    
    def fit(self, X_org, Y_org, n_tests, n_epochs,
            save_subject = False, ablation = True, ablation_target = ['W', 'G', 'cog', 'C1', 'C2','C3',],
            test_size = 0.2, random_state = 9999, 
            feature_importance = False, step_importance = False, grad_cam = False, record_best=False):
        assert len(self.feature_names) == X_org.shape[1]
        seed =  [i for i in range(random_state, random_state - n_tests, -1)]
        if ablation:
            index, query = ablation_ts(X_org, ablation_target, self.feature_categories)
            print(index, query)
            X_raw = X_org[:, index]
            self.feature_names = query
        Y_raw = 1 - Y_org
        for i in range(n_tests):
            if i > 0:
                for j in range(len(self.metrics)):
                    print(f'{i}/{n_tests} {self.metrics[j]}: {np.mean(self._result[j]):.3f} +/- {np.std(self._result[j]):.3f}')    
            else: print(f'{i}/{n_tests}')
            if self.class_balance:
                X, Y = class_balance(X_raw, Y_raw, seed = seed[i])
            else:
                X, Y = X_raw, Y_raw
            
            if save_subject:
                if not os.path.exists(subjectpath):
                    os.makedirs(subjectpath)
                np.save(os.path.join(self.subject_path,  f'{i}_subjectlist.npy'), subjectList)
                np.save(os.path.join(self.subject_path,  f'{i}_labels.npy'), Y)
            

            
            X_train , X_test, y_train, y_test = train_test_split(X, Y, stratify = Y, test_size = test_size, random_state = random_state)
            X, y, splits = combine_split_data([X_train, X_test], [y_train, y_test])
            self.X_test = X_test
            self.y_test = y_test
            tfms = [None, TSClassification()]
            dls = get_ts_dls(X, y, splits=splits, tfms=tfms) 
            model = self.get_model(dls)
            if self.pretrain:
                learn = ts_learner(dls, model, loss_func= self.loss_fn,
                                pretrained=True, weights_path=f'data/MVP/xception.pth', 
                                metrics=[Precision(), Recall(), accuracy, RocAucBinary(), F1Score(), accuracy])
                self._learn = learn
                learn.fine_tune(n_epochs, base_lr=2e-2, freeze_epochs=freeze_epochs)
            else:
                learn = ts_learner(dls, model, 
                        loss_func = self.loss_fn,
                        # metrics = accuracy,)
                        metrics=[Precision(), Recall(), accuracy, RocAucBinary(), F1Score(), BalancedAccuracy()])
                self._learn = learn
                learn.fit_one_cycle(n_epochs, self.lr)

            if feature_importance:
                self.feature_importance(prefix = str(i))
            if step_importance:
                self.step_importance(prefix = str(i), n_steps = self.config['step_importance_steps'])
            if grad_cam:
                self.grad_cam(X, prefix = str(i), size = self.config['Grad_CAM_size'], layer = self.config['Grad_CAM_layer'])
            
            if record_best:
                aurocs = np.array([learn.recorder.values[i][2+3] for i in range(n_epochs)])
                index = np.argmax(aurocs)
            for j in range(len(self.metrics)):
                if not record_best:
                    self._result[j].append(learn.recorder.values[-1][2+j])
                else:
                    self._result[j].append(learn.recorder.values[index][2+j])
        
        print('>>>>>>>>>>>>>>> Final Result <<<<<<<<<<<<<<<<<<')
        for j in range(len(self.metrics)):
            print(f'\n{self.metrics[j]}: {np.mean(self._result[j]):.3f} +/- {np.std(self._result[j]):.3f} in {n_tests} tests')
            print(self._result[j])


    def fit_sklearn(self, X_org, Y_org, n_tests, classifier = 'XGBoost',
            save_subject = False, ablation = True, ablation_target = ['W', 'G', 'cog', 'C1', 'C2','C3',],
            test_size = 0.2, random_state = 9999, ):
        
        assert len(self.feature_names) == X_org.shape[1]
        seed =  [i for i in range(random_state, random_state - n_tests, -1)]
        if ablation:
            index, query = ablation_ts(X_org, ablation_target, self.feature_categories)
            print(index, query)
            X_raw = X_org[:, index]
            self.feature_names = query
        Y_raw = 1 - Y_org
        for i in range(n_tests):
            if i > 0:
                for j in range(len(self.metrics)):
                    print(f'{i}/{n_tests} {self.metrics[j]}: {np.mean(self._result[j]):.3f} +/- {np.std(self._result[j]):.3f}')    
            else: print(f'{i}/{n_tests}')
            if self.class_balance:
                X, Y = class_balance(X_raw, Y_raw, seed = seed[i])
            else:
                X, Y = X_raw, Y_raw
            
            if save_subject:
                if not os.path.exists(subjectpath):
                    os.makedirs(subjectpath)
                np.save(os.path.join(self.subject_path,  f'{i}_subjectlist.npy'), subjectList)
                np.save(os.path.join(self.subject_path,  f'{i}_labels.npy'), Y)
            

            
            X_train , X_test, y_train, y_test = train_test_split(X, Y, stratify = Y, test_size = test_size, random_state = random_state)
            X_train = np.mean(X_train, axis = 2)
            X_test = np.mean(X_test, axis = 2)
            # Import the necessary modules
            # Import the necessary modules
            from sklearn.ensemble import RandomForestClassifier
            from sklearn.metrics import precision_score, recall_score, accuracy_score, roc_auc_score, f1_score, balanced_accuracy_score

            if classifier =='XGBoost':
                from xgboost import XGBClassifier
                clf = XGBClassifier()
            elif classifier == 'RF':
            # Create a random forest classifier with default parameters
                clf = RandomForestClassifier()
            else:
                raise NotImplementedError

            # Fit the classifier on the training data
            clf.fit(X_train, y_train)

            # Predict on the test data
            y_pred = clf.predict(X_test)
            y_prob = clf.predict_proba(X_test)[:, 1]

            # Calculate and print the metrics
            precision = precision_score(y_test, y_pred)
            recall = recall_score(y_test, y_pred)
            accuracy = accuracy_score(y_test, y_pred)
            auroc = roc_auc_score(y_test, y_prob)
            f1 = f1_score(y_test, y_pred)
            balanced_accuracy = balanced_accuracy_score(y_test, y_pred)
            result = [precision, recall, accuracy, auroc, f1, balanced_accuracy]
            for j in range(len(self.metrics)):
                self._result[j].append(result[j])
            print(f"Precision: {precision:.2f}")
            print(f"Recall: {recall:.2f}")
            print(f"Accuracy: {accuracy:.2f}")
            print(f"AUROC: {auroc:.2f}")
            print(f"F1: {f1:.2f}")
            print(f"Balanced Accuracy: {balanced_accuracy:.2f}")        
        print('>>>>>>>>>>>>>>> Final Result <<<<<<<<<<<<<<<<<<')
        for j in range(len(self.metrics)):
            print(f'\n{self.metrics[j]}: {np.mean(self._result[j]):.3f} +/- {np.std(self._result[j]):.3f} in {n_tests} tests')
            print(self._result[j])



    def grad_cam(self, X, size = 32, layer = 0, save = True, prefix = ''):
        '''
        X of size shape [N, C, L]
        Average of gram cap values with
        '''
        print(f'Performing Gard CAM of {prefix}...')
        data = torch.Tensor(X[0:0+size]).cuda()
        flag = False
        out = torch.zeros(X.shape[layer], X.shape[2])
        with Hook(self._learn.model[layer]) as hook:
            for i in range(data.shape[0]):
                with HookBwd(self._learn.model[layer]) as hookg:
                    output = self._learn.model.eval()(data) # last layer output [N, 1, num_classes]
                    act = hook.stored #act [N, C, L]
                    cls = torch.argmax(output, axis = 1) # [N, 1]
                    output[i,cls[i]].backward(retain_graph=True) # backpropagate each one
                    if not flag:
                        grad = hookg.stored # gradient [N, C, L]
                        flag = True
                    else:
                        grad += hookg.stored
        with torch.no_grad():
            w = grad.mean(dim=[2], keepdim=True) # average gradient by length [N, C, L(summed out)]
            
            cam_map = F.relu(w * act).sum((0,1)).detach().cpu().numpy() 
        if save:
            plt.plot(cam_map)
            plt.title('GradCAM Plot')
            plt.savefig(os.path.join(self.save_path, f"{prefix}_GradCAM.pdf"))
        print('Done') 

    def feature_importance(self, random_state = 44, method = 'permutation', save = True, prefix = ''):
        print(f'Performing Feature Importance of {prefix}...')
        self._learn.metrics = RocAucBinary()
        df = self._learn.feature_importance(self.X_test, self.y_test, random_state = 44, feature_names = self.feature_names, method = method)
        if save:
            df.to_csv(os.path.join(self.save_path,f'{prefix}_{self.model}_feature_importance_{method}.csv'))
        print('Done') 
    
    def step_importance(self, random_state = 44, n_steps=60, method = 'permutation', save = True, prefix = ''):
        print(f'Performing Step Importance of {prefix}...')
        self._learn.metrics = RocAucBinary()
        df = self._learn.step_importance(n_steps=n_steps, random_state = random_state, method = method)
        if save:
            df.to_csv(os.path.join(self.save_path,f'{prefix}_{self.model}_step_importance_{method}.csv'))
        print('Done') 

if __name__ == '__main__':
    X = np.random.randn(300, 2, 50)
    Y = np.zeros(300)
    Y[200:] = 1
    X_train , X_test, y_train, y_test = train_test_split(X, Y,  test_size = 0.2, random_state = 333)
    X, y, splits = combine_split_data([X_train, X_test], [y_train, y_test])
    tfms = [None, TSClassification()]
    dls = get_ts_dls(X, y, splits=splits, tfms=tfms) 
    learn = ts_learner(dls, InceptionTime, 
            loss_func = LabelSmoothingCrossEntropyFlat(),
            # metrics = accuracy,)
            metrics= RocAucBinary())
    learn.fit_one_cycle(2, 1e-4)
    learn.step_importance()
        
        