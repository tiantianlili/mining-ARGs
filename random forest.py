import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend to avoid tkinter errors
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import LabelEncoder, StandardScaler, MinMaxScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import RFE
from scipy.stats import kruskal

# ============================================================================
# 1. Data Loading
# ============================================================================
arg_data = pd.read_csv('input_data.txt', sep='\t')
meta_data = pd.read_csv('metadata.txt', sep='\t')

# ============================================================================
# 2. Remove columns with over 80% zero values
# ============================================================================
zero_counts = (arg_data == 0).sum()
threshold_80pct = 0.8 * len(arg_data)
cols_to_drop = zero_counts[zero_counts > threshold_80pct].index
cols_to_drop = [c for c in cols_to_drop if c != 'Sample']  # Preserve 'Sample' column
arg_data.drop(columns=cols_to_drop, inplace=True)

# ============================================================================
# 3. Define target variable: all_arg_abundance (sum of remaining gene columns)
# ============================================================================
arg_data['all_arg_abundance'] = arg_data.drop(columns=['Sample']).sum(axis=1)

# ============================================================================
# 4. Data Merging
# ============================================================================
data = pd.merge(meta_data, arg_data[['Sample', 'all_arg_abundance']], on='Sample')

# Separate features (X) and target (y)
X = data.drop(columns=['Sample', 'all_arg_abundance'])
y = data['all_arg_abundance']

# ============================================================================
# 5. Stratified dataset split by group
# ============================================================================
if 'group' not in X.columns:
    raise ValueError("Column 'group' not found in metadata. Please verify column names.")

# Perform stratified split
X_train_full, X_test_full, y_train_full, y_test_full = train_test_split(
    X, 
    y, 
    test_size=0.3, 
    random_state=42, 
    stratify=X['group']
)

print("Training set group distribution:")
print(X_train_full['group'].value_counts())
print("Test set group distribution:")
print(X_test_full['group'].value_counts())

# ============================================================================
# 6. Numeric data normalization and outlier removal
# ============================================================================
# (1) Identify categorical/numeric columns
categorical_columns = ['MineT', 'LandM', 'group', 'Site', 'Geology', 'Soil type']
num_cols = [col for col in X_train_full.columns if col not in categorical_columns]

# (2) Remove outliers using IQR method (training set only)
X_train_num = X_train_full[num_cols].copy()
y_train = y_train_full.copy()
Q1 = X_train_num.quantile(0.25)
Q3 = X_train_num.quantile(0.75)
IQR = Q3 - Q1
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

mask = pd.DataFrame(True, index=X_train_num.index, columns=X_train_num.columns)
for col in X_train_num.columns:
    mask[col] = (X_train_num[col] >= lower_bound[col]) & (X_train_num[col] <= upper_bound[col])
valid_index_train = mask.all(axis=1)

X_train_full = X_train_full.loc[valid_index_train]
y_train_full = y_train.loc[valid_index_train]

# (3) Normalization/Standardization
scaler = StandardScaler()
X_train_num_scaled = scaler.fit_transform(X_train_full[num_cols])
X_train_full[num_cols] = X_train_num_scaled

# (4) Apply transformation to test set
X_test_full_num = X_test_full[num_cols].copy()
X_test_full_num_scaled = scaler.transform(X_test_full_num)
X_test_full[num_cols] = X_test_full_num_scaled

# ============================================================================
# 7. RFE Feature Selection
# ============================================================================
# Encode categorical variables
for col in categorical_columns:
    if col in X_train_full.columns:
        le = LabelEncoder()
        X_train_full[col] = le.fit_transform(X_train_full[col])
        X_test_full[col] = le.transform(X_test_full[col])

# Prepare final datasets
X_train_final = X_train_full.drop(columns=['group'])
X_test_final = X_test_full.drop(columns=['group'])
y_train_final = y_train_full
y_test_final = y_test_full

# Perform RFE
rfe_estimator = RandomForestRegressor(n_estimators=100, random_state=0, n_jobs=-1)
n_features_to_select = max(1, X_train_final.shape[1] // 2)

rfe = RFE(estimator=rfe_estimator, n_features_to_select=n_features_to_select)
rfe.fit(X_train_final, y_train_final)

selected_cols = X_train_final.columns[rfe.support_]
print(f"RFE selected features count: {len(selected_cols)}")
print(f"RFE selected features: {list(selected_cols)}")

X_train_sel = X_train_final[selected_cols]
X_test_sel = X_test_final[selected_cols]

# ============================================================================
# 8. Hyperparameter Tuning
# ============================================================================
param_candidates = {
    'n_estimators': [100, 500, 1000, 2000],
    'max_depth': [5, 10, 20],
    'min_samples_split': [2, 4],
    'min_samples_leaf': [1, 2, 4]
}

best_params = None
best_model = None
achieved_criteria = False

# Parameter search loop
for n_est in param_candidates['n_estimators']:
    for md in param_candidates['max_depth']:
        for mss in param_candidates['min_samples_split']:
            for msl in param_candidates['min_samples_leaf']:
                
                model = RandomForestRegressor(
                    n_estimators=n_est,
                    max_depth=md,
                    min_samples_split=mss,
                    min_samples_leaf=msl,
                    random_state=42,
                    n_jobs=-1
                )
                
                model.fit(X_train_sel, y_train_final)
                y_pred_train = model.predict(X_train_sel)
                mse_train = mean_squared_error(y_train_final, y_pred_train)
                r2_train = r2_score(y_train_final, y_pred_train)
                
                if mse_train < 0.001 and r2_train > 0.6:
                    best_params = {
                        'n_estimators': n_est,
                        'max_depth': md,
                        'min_samples_split': mss,
                        'min_samples_leaf': msl
                    }
                    best_model = model
                    achieved_criteria = True
                    print("Target metrics achieved:", best_params)
                    print("Training MSE:", mse_train)
                    print("Training R²:", r2_train)
                    break
            if achieved_criteria:
                break
        if achieved_criteria:
            break
    if achieved_criteria:
        break

# Fallback if criteria not met
if not achieved_criteria:
    print("Target metrics not achieved. Consider expanding search space or relaxing thresholds.")
    best_model = model  
    best_params = {
        'n_estimators': n_est,
        'max_depth': md,
        'min_samples_split': mss,
        'min_samples_leaf': msl
    }

# Final training metrics
y_pred_train_final = best_model.predict(X_train_sel)
train_mse_final = mean_squared_error(y_train_final, y_pred_train_final)
train_r2_final = r2_score(y_train_final, y_pred_train_final)

print("Final training parameters:", best_params)
print("Final training MSE:", train_mse_final)
print("Final training R²:", train_r2_final)

# ============================================================================
# 9. Test Set Evaluation
# ============================================================================
y_pred_test = best_model.predict(X_test_sel)
test_mse = mean_squared_error(y_test_final, y_pred_test)
test_r2 = r2_score(y_test_final, y_pred_test)
print("Test set MSE:", test_mse)
print("Test set R²:", test_r2)

# ============================================================================
# 10. Visualization: Train vs Test Metrics
# ============================================================================
plt.figure(figsize=(6, 4))
plt.bar(['Train', 'Test'], [train_r2_final, test_r2], color=['#1f77b4', '#ff7f0e'])
plt.title('R² Comparison: Train vs Test')
plt.ylabel('R²')
plt.tight_layout()
plt.savefig('r2_comparison.svg')
plt.close()

plt.figure(figsize=(6, 4))
plt.bar(['Train', 'Test'], [train_mse_final, test_mse], color=['#2ca02c', '#d62728'])
plt.title('MSE Comparison: Train vs Test')
plt.ylabel('MSE')
plt.tight_layout()
plt.savefig('mse_comparison.svg')
plt.close()

# ============================================================================
# 11. Group-wise Analysis
# ============================================================================
# a) Overall feature importance
feature_importances = pd.DataFrame({
    'Feature': selected_cols,
    'Importance': best_model.feature_importances_
}).sort_values('Importance', ascending=False)
feature_importances.to_csv('feature_importance.csv', index=False)

# b) Group-specific analysis
groups_to_check = ['P', 'CK', 'F']
all_group_importances = []

for grp in groups_to_check:
    grp_mask = X_train_full['group'] == grp
    X_grp = X_train_full.loc[grp_mask, selected_cols]
    y_grp = y_train_full.loc[grp_mask]
    
    if len(X_grp) < 5:
        print(f"Insufficient samples in group {grp}. Skipping analysis.")
        continue
    
    model_grp = RandomForestRegressor(**best_params, random_state=42, n_jobs=-1)
    model_grp.fit(X_grp, y_grp)

    grp_importances = pd.DataFrame({
        'Feature': selected_cols,
        'Importance': model_grp.feature_importances_,
        'Group': grp
    })
    all_group_importances.append(grp_importances)

# c) Significance testing
if len(all_group_importances) > 0:
    group_importance_df = pd.concat(all_group_importances, axis=0)
    group_importance_df.to_csv('group_feature_importance.csv', index=False)

    # Kruskal-Wallis test
    unique_features = group_importance_df['Feature'].unique()
    p_values = []
    for f_col in unique_features:
        data_p = group_importance_df.loc[group_importance_df['Feature'] == f_col]
        groups_data = [grp_data['Importance'] for _, grp_data in data_p.groupby('Group')]
        
        if len(groups_data) < 2:
            p_val = np.nan
        else:
            _, p_val = kruskal(*groups_data)
        p_values.append([f_col, p_val])
    
    sig_df = pd.DataFrame(p_values, columns=['Feature', 'p_value'])
    sig_df.to_csv('feature_significance.csv', index=False)

    # Visualization
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=group_importance_df, x='Feature', y='Importance', hue='Group')
    plt.xticks(rotation=90)
    plt.title('Feature Importance by Group')
    plt.tight_layout()
    plt.savefig('group_feature_importance.svg')
    plt.close()
else:
    print("Insufficient data for group analysis.")

print("Script execution completed. All outputs saved.")
