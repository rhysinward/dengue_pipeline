import pandas as pd

# Load the metadata CSV file
df = pd.read_csv('results/subsampled_Dengue_1_cleaned_infoTbl.csv')

# Check for duplicates in the index column (assuming it's the first column)
if df.iloc[:,0].duplicated().any():
    print("Duplicates found. Removing...")
    # Remove duplicates, keeping the first occurrence
    df = df.drop_duplicates(subset=df.columns[0], keep='first')
    # Save the cleaned file
    df.to_csv('results/subsampled_Dengue_1_cleaned_infoTbl_cleaned.csv', index=False)
    print("Cleaned file saved.")
else:
    print("No duplicates found.")
