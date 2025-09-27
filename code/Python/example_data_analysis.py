"""
Example Data Analysis Script for Computational Economics
This script demonstrates basic data analysis techniques using Python
Course: Econ-81360, Fall 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Set style for plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def generate_economic_data(n=1000, seed=42):
    """Generate synthetic economic data for demonstration"""
    np.random.seed(seed)
    
    # Generate correlated variables
    income = np.random.normal(50000, 15000, n)
    education = np.random.normal(16, 3, n)
    experience = np.random.exponential(10, n)
    
    # Create consumption with some realistic relationships
    consumption = (0.7 * income + 
                  2000 * education + 
                  500 * experience + 
                  np.random.normal(0, 5000, n))
    
    # Ensure positive values
    income = np.maximum(income, 10000)
    education = np.maximum(education, 8)
    consumption = np.maximum(consumption, 5000)
    
    # Create DataFrame
    data = pd.DataFrame({
        'income': income,
        'education': education,
        'experience': experience,
        'consumption': consumption
    })
    
    return data

def analyze_data(data):
    """Perform basic data analysis"""
    print("=== Data Analysis Results ===")
    print(f"Dataset shape: {data.shape}")
    print("\nDescriptive Statistics:")
    print(data.describe())
    
    print("\nCorrelation Matrix:")
    correlation_matrix = data.corr()
    print(correlation_matrix)
    
    return correlation_matrix

def create_visualizations(data, correlation_matrix):
    """Create various plots for data visualization"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Economic Data Analysis Visualizations', fontsize=16)
    
    # 1. Scatter plot: Income vs Consumption
    axes[0, 0].scatter(data['income'], data['consumption'], alpha=0.6)
    axes[0, 0].set_xlabel('Income ($)')
    axes[0, 0].set_ylabel('Consumption ($)')
    axes[0, 0].set_title('Income vs Consumption')
    
    # Add regression line
    z = np.polyfit(data['income'], data['consumption'], 1)
    p = np.poly1d(z)
    axes[0, 0].plot(data['income'], p(data['income']), "r--", alpha=0.8)
    
    # 2. Histogram: Income distribution
    axes[0, 1].hist(data['income'], bins=30, alpha=0.7, edgecolor='black')
    axes[0, 1].set_xlabel('Income ($)')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Income Distribution')
    
    # 3. Box plot: Education levels
    axes[1, 0].boxplot(data['education'])
    axes[1, 0].set_ylabel('Years of Education')
    axes[1, 0].set_title('Education Distribution')
    
    # 4. Correlation heatmap
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0,
                ax=axes[1, 1])
    axes[1, 1].set_title('Correlation Heatmap')
    
    plt.tight_layout()
    plt.savefig('economic_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def regression_analysis(data):
    """Perform regression analysis"""
    print("\n=== Regression Analysis ===")
    
    # Prepare data
    X = data[['income', 'education', 'experience']]
    y = data['consumption']
    
    # Fit linear regression model
    model = LinearRegression()
    model.fit(X, y)
    
    # Predictions
    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    
    # Results
    print(f"R-squared: {r2:.4f}")
    print("Coefficients:")
    for i, col in enumerate(X.columns):
        print(f"  {col}: {model.coef_[i]:.4f}")
    print(f"Intercept: {model.intercept_:.4f}")
    
    # Statistical tests
    residuals = y - y_pred
    mse = np.mean(residuals**2)
    rmse = np.sqrt(mse)
    
    print(f"\nModel Performance:")
    print(f"  MSE: {mse:.2f}")
    print(f"  RMSE: {rmse:.2f}")
    
    return model, y_pred, residuals

def main():
    """Main function to run the analysis"""
    print("Economic Data Analysis Example")
    print("=" * 40)
    
    # 1. Generate data
    print("Generating synthetic economic data...")
    data = generate_economic_data()
    
    # 2. Basic analysis
    correlation_matrix = analyze_data(data)
    
    # 3. Visualizations
    print("\nCreating visualizations...")
    create_visualizations(data, correlation_matrix)
    
    # 4. Regression analysis
    model, predictions, residuals = regression_analysis(data)
    
    # 5. Save results
    data.to_csv('economic_data.csv', index=False)
    print("\nData saved to 'economic_data.csv'")
    print("Visualization saved as 'economic_analysis.png'")

if __name__ == "__main__":
    main()