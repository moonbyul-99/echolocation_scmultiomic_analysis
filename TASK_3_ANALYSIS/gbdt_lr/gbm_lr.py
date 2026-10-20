import numpy as np 
import pandas as pd 
import scanpy as sc 
import lightgbm as lgb 
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score 
from sklearn.preprocessing import LabelEncoder, label_binarize
from sklearn.model_selection import train_test_split
from tqdm import tqdm 

if __name__ == '__main__':
    scdata = sc.read_h5ad('/home/rsun@ZHANGroup.local/sly_data/data/final_hvg.h5ad')
    scdata = scdata.raw.to_adata()
    
    # define label encoder
    labels = scdata.obs.new_anno.values.astype(str)
    label_encoder = LabelEncoder()
    numeric_labels = label_encoder.fit_transform(labels) 

    gbm_id, lg_id = train_test_split(np.arange(scdata.shape[0]), test_size=0.5, random_state=42, stratify= scdata.obs.new_anno.values.astype(str))
    gbm_data = scdata.X[gbm_id,:].toarray()
    lg_data = scdata.X[lg_id,:].toarray() 

    gbm_label = scdata.obs.new_anno.values[gbm_id]
    lg_label = scdata.obs.new_anno.values[lg_id]

    gbm_label = label_encoder.transform(gbm_label)
    lg_label = label_encoder.transform(lg_label)
    print(gbm_data.shape, lg_data.shape)    

    lgb_train = lgb.Dataset(gbm_data, label=gbm_label)
    lgb_eval = lgb.Dataset(lg_data, lg_label, reference=lgb_train)

    C = label_encoder.classes_.shape[0]

    # 设置LightGBM参数
    params = {
        'objective': 'multiclass',  # 多分类任务
        'num_class': C,  # 类别数量
        'metric': 'multi_logloss',  # 多分类的损失函数
        'boosting_type': 'gbdt',  # 使用GBDT算法
        'num_leaves': 31,  # 叶子节点数
        'learning_rate': 0.05,  # 学习率
        'feature_fraction': 0.9,  # 特征采样比例
        'bagging_fraction': 0.8,  # 数据采样比例
        'bagging_freq': 5,  # 每5次迭代进行一次bagging
        'verbose': 0  # 不输出详细信息
    }

    # 训练模型
    early_stopping_callback = lgb.early_stopping(stopping_rounds=10, verbose=True)
    model = lgb.train(params, lgb_train, num_boost_round=100, valid_sets=[lgb_eval], 
                    callbacks=[early_stopping_callback])

    # 在验证集上进行预测
    print(f'best iteration: {model.best_iteration}')
    y_pred_prob = model.predict(lg_data, num_iteration=model.best_iteration)  # 返回每个类别的概率

    # 计算多分类的ROC AUC Score
    lg_label_bin = label_binarize(lg_label, classes=np.unique(lg_label))
    roc_auc = roc_auc_score(lg_label, y_pred_prob, multi_class='ovr')  # 使用'ovr'策略

    print(f"ROC AUC Score on Evaluation Data: {roc_auc}")

    print('Start predicting...')
    # predict and get data on leaves, training data
    y_pred = model.predict(gbm_data, pred_leaf=True)

    print(np.array(y_pred).shape)
    print(y_pred[:10])

    print('Writing transformed training data')
    transformed_training_matrix = np.zeros([len(y_pred), len(y_pred[0]) * params['num_leaves']],
                                        dtype=np.int64)  # N * num_tress * num_leafs
    for i in tqdm(range(0, len(y_pred))):
        temp = np.arange(len(y_pred[0])) * params['num_leaves'] + np.array(y_pred[i])
        transformed_training_matrix[i][temp] += 1

    y_pred = model.predict(lg_data, pred_leaf=True)
    print('Writing transformed testing data')
    transformed_testing_matrix = np.zeros([len(y_pred), len(y_pred[0]) * params['num_leaves']], dtype=np.int64)
    for i in tqdm(range(0, len(y_pred))):
        temp = np.arange(len(y_pred[0])) * params['num_leaves'] + np.array(y_pred[i])
        transformed_testing_matrix[i][temp] += 1

    lm = LogisticRegression(penalty='l2', C=0.1) # logestic model construction
    lm.fit(transformed_training_matrix, gbm_label)  # fitting the data
    y_pred_test = lm.predict_proba(transformed_testing_matrix)   # Give the probabilty on each label

    print(y_pred_test.shape)

    # 获取Logistic Regression的权重
    coef = lm.coef_
    # 计算每个特征的重要性（权重的绝对值之和）
    feature_importance = np.sum(np.abs(coef), axis=0)
    # 获取重要特征的索引并按重要性排序
    important_feature_indices = np.argsort(feature_importance)[::-1][:100]
    print(f"Number of important features: {len(important_feature_indices)}")

    # 将GBM生成的特征转换为原始特征名称组合的形式
    leaf_indices = np.arange(params['num_leaves'])
    feature_names = [f'leaf_{tree}_{leaf}' for tree in range(len(y_pred[0])) for leaf in leaf_indices]
    important_feature_names = [feature_names[i] for i in important_feature_indices]
    print("Important feature names:")
    print(important_feature_names)

    # 将GBM生成的特征映射回原始特征名称组合
    original_feature_names = scdata.var.index
    feature_map = {}
    for tree in range(len(y_pred[0])):
        for leaf in leaf_indices:
            leaf_index = tree * params['num_leaves'] + leaf
            feature_map[f'leaf_{tree}_{leaf}'] = []
            for i in range(len(gbm_data)):
                if y_pred[i][tree] == leaf:
                    feature_map[f'leaf_{tree}_{leaf}'].append(original_feature_names[i])

    # 打印重要特征的原始特征名称组合
    for feature_name in important_feature_names:
        print(f"{feature_name}: {feature_map[feature_name]}")


def parse_lightgbm_tree(tree, feature_names=None):
    """
    解析 LightGBM 树结构，生成叶节点路径条件表达式
    :param tree: 通过 model.dump_model()['tree_info'][0]['tree_structure'] 获取的树结构
    :param feature_names: 特征名称列表，若为 None 则使用 column_{index} 格式
    :return: 字典，键为叶节点索引，值为条件表达式字符串
    """
    def traverse(node, path, results):
        if 'leaf_index' in node:
            # 处理叶节点
            conditions = []
            for feat, thres, op in path:
                conditions.append(f"{feat} {op} {thres:.3f}")
            results[node['leaf_index']] = " and ".join(conditions)
        else:
            # 处理非叶节点
            feat_idx = node['split_feature']
            thres = node['threshold']
            feat_name = f"column_{feat_idx}" if feature_names is None else feature_names[feat_idx]
            
            # 左子树（<= 阈值）
            left_path = path + [(feat_name, thres, "<=")]
            traverse(node.get('left_child', {}), left_path, results)
            
            # 右子树（> 阈值）
            right_path = path + [(feat_name, thres, ">")]
            traverse(node.get('right_child', {}), right_path, results)
            
    results = {}
    traverse(tree, [], results)
    return results