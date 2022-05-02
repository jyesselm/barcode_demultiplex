import pandas as pd

def main():
    df = pd.read_csv("test.csv")
    df = df.drop_duplicates(['barcode1', 'barcode2'])
    df.to_csv("test_uniqe.csv", index=False)

if __name__ == '__main__':
    main()