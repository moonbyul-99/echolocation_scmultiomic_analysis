import numpy as np 
import pandas as pd 
import scanpy as sc 
import lightgbm as lgb 
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score 
from sklearn.preprocessing import LabelEncoder, label_binarize
from sklearn.model_selection import train_test_split
from tqdm import tqdm 
import anndata as ad

'''
Basic functions of gbdt+lr model 
''' 

def load_data(path, use_raw = True):
    '''
    Args:
        path: path of h5ad file eg. /home/rsun@ZHANGroup.local/sly_data/data/final_hvg.h5ad
        use_raw: whether to use raw data or not 
    Returns:
        scdata: anndata object
    '''

    scdata = sc.read_h5ad(path)
    if use_raw:
        scdata = scdata.raw.to_adata()
    return scdata

def get_label_encoder(scdata, label_key = 'new_anno'):
    '''
    Args:
        scdata: anndata object
        label_key: key of label column in obs 
    Returns:
        scdata: anndata object with label encoded 
    '''

    # define label encoder
    labels = scdata.obs[label_key].values.astype(str)
    label_encoder = LabelEncoder()
    label_encoder.fit(labels)
    return label_encoder 

class GBDT_LR:

    def __init__(self, scdata, label_key):
        self.scdata = scdata 
        self.label_key = label_key
        self.label_encoder = get_label_encoder(scdata, label_key)
    
    def split_data(self, gbm_size = 0.5, random_state = 42):
        '''
        split data into two parts,

        one for gbdt classifier training, 
        one for lr classifier training
        '''

        gbm_id, lg_id = train_test_split(np.arange(self.scdata.shape[0]), 
                                         test_size=gbm_size, 
                                         random_state=random_state, 
                                         stratify= self.scdata.obs.loc[:,self.label_key].values.astype(str))

        gbm_data = self.scdata.X[gbm_id,:].toarray()
        lg_data = self.scdata.X[lg_id,:].toarray() 

        #gbm_label = self.scdata.obs.new_anno.values[gbm_id]
        #lg_label = self.scdata.obs.new_anno.values[lg_id]
        gbm_label = self.scdata.obs.loc[:,self.label_key].values[gbm_id]
        lg_label = self.scdata.obs.loc[:,self.label_key].values[lg_id]

        gbm_label = self.label_encoder.transform(gbm_label)
        lg_label = self.label_encoder.transform(lg_label)

        split_res = {'gbm_id': gbm_id, 
                    'lg_id': lg_id, 
                    'gbm_data': gbm_data, 
                    'lg_data': lg_data, 
                    'gbm_label': gbm_label, 
                    'lg_label': lg_label}
        self.split_res = split_res
        return None
    
    def gbm_fit(self, gbm_params = None):
        '''
        fit gbdt model
        '''

        if gbm_params is None:
            C = self.label_encoder.classes_.shape[0]

            # 设置LightGBM参数
            gbm_params = {
                'objective': 'multiclass',  # 多分类任务
                'num_class': C,  # 类别数量
                'metric': 'multi_logloss',  # 多分类的损失函数
                'boosting_type': 'gbdt',  # 使用GBDT算法
                'num_leaves': 15,  # 叶子节点数
                'learning_rate': 0.05,  # 学习率
                'feature_fraction': 0.9,  # 特征采样比例
                'bagging_fraction': 0.8,  # 数据采样比例
                'bagging_freq': 5,  # 每5次迭代进行一次bagging
                'verbose': 0  # 不输出详细信息
                }
        self.gbm_params = gbm_params 

        gbm_data = self.split_res['gbm_data']
        gbm_label = self.split_res['gbm_label']
        lg_data = self.split_res['lg_data']
        lg_label = self.split_res['lg_label']

        lgb_train = lgb.Dataset(gbm_data, label=gbm_label)
        lgb_eval = lgb.Dataset(lg_data, lg_label, reference=lgb_train)

        # 训练模型
        early_stopping_callback = lgb.early_stopping(stopping_rounds=10, verbose=True)
        model = lgb.train(gbm_params, lgb_train, num_boost_round=100, valid_sets=[lgb_eval], 
                        callbacks=[early_stopping_callback])

        # 在验证集上进行预测
        print(f'best iteration: {model.best_iteration}')
        y_pred_prob = model.predict(lg_data, num_iteration=model.best_iteration)  # 返回每个类别的概率

        # 计算多分类的ROC AUC Score
        roc_auc = roc_auc_score(lg_label, y_pred_prob, multi_class='ovr')  # 使用'ovr'策略
        print(f"ROC AUC Score on lg Data: {roc_auc}")
        
        self.gbm_model = model 
        return None 
        
    def gbm_transform(self):
        gbm_data = self.split_res['gbm_data']
        lg_data = self.split_res['lg_data']

        num_leaves = self.gbm_params['num_leaves']

        print('Start predicting...')
        # predict and get data on leaves, training data
        y_pred = self.gbm_model.predict(gbm_data, pred_leaf=True)

        print('Writing transformed training data')
        transformed_training_matrix = np.zeros([len(y_pred), len(y_pred[0]) * num_leaves],
                                            dtype=np.int64)  # N * num_tress * num_leafs
        for i in tqdm(range(0, len(y_pred))):
            temp = np.arange(len(y_pred[0])) * num_leaves + np.array(y_pred[i])
            transformed_training_matrix[i][temp] += 1


        y_pred = self.gbm_model.predict(lg_data, pred_leaf=True)
        print('Writing transformed testing data')
        transformed_testing_matrix = np.zeros([len(y_pred), len(y_pred[0]) * num_leaves], dtype=np.int64)
        for i in tqdm(range(0, len(y_pred))):
            temp = np.arange(len(y_pred[0])) * num_leaves + np.array(y_pred[i])
            transformed_testing_matrix[i][temp] += 1

        """
        set the transformed data into anndata object
        """

        new_data =  np.concatenate([transformed_training_matrix, transformed_testing_matrix])
        new_id = np.concatenate([self.split_res['gbm_id'], self.split_res['lg_id']])
        new_obs = self.scdata.obs.iloc[new_id,:].copy()

        new_obs.loc[:, 'split_label'] = ['gbm']*len(self.split_res['gbm_id']) + ['lg']*len(self.split_res['lg_id'])
        new_Umap = self.scdata.obsm['X_umap'][new_id,:]
        
        sc_new = sc_new = ad.AnnData(new_data, obs = new_obs)
        sc_new.obsm['X_umap'] = new_Umap 

        self.sc_new = sc_new 

        return None 
    
    def lg_fit(self):
        testing_id = np.arange(self.sc_new.shape[0])[self.sc_new.obs.split_label == 'lg']
        testing_label = self.sc_new.obs.loc[:,self.label_key][testing_id]
        transformed_testing_matrix = self.sc_new.X[testing_id, :]

        lg = LogisticRegression(penalty='l2',C=0.1) # logestic model construction
        lg.fit(transformed_testing_matrix,testing_label)  # fitting the data

        y_pred_proba = lg.predict_proba(transformed_testing_matrix)   # Give the probabilty on each label
        roc_auc = roc_auc_score(testing_label, y_pred_proba, multi_class='ovr')  # 使用'ovr'策略
        print(f"ROC AUC Score on lg split Data: {roc_auc}")


        training_id = np.arange(self.sc_new.shape[0])[self.sc_new.obs.split_label == 'gbm']
        training_label = self.sc_new.obs.loc[:,self.label_key][training_id]
        transformed_training_matrix = self.sc_new.X[training_id, :]

        y_pred_proba = lg.predict_proba(transformed_training_matrix)   # Give the probabilty on each label
        roc_auc = roc_auc_score(training_label, y_pred_proba, multi_class='ovr')  # 使用'ovr'策略
        print(f"ROC AUC Score on gbm split Data: {roc_auc}")

        self.lg_model = lg 
        return None



