import pickle
import argparse
import pprint
import pandas as pd

# pd.DataFrame.il
def main():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('subject', type=str, help="subject identifier")
    args = parser.parse_args()

    file = open (f"data/{args.subject}/blueberry/0timing.pkl", "rb")
    data = pickle.load(file)
    pprint.pprint(data)

    print(data['visual_regular_end'] - data['visual_regular_start'])
    print(data['auditory_regular_end'] - data['auditory_regular_start'])
    print(data['mental_regular_end'] - data['mental_regular_start'])
    print('')
    print(data['visual_random_end'] - data['visual_random_start'])
    print(data['auditory_random_end'] - data['auditory_random_start'])
    print(data['mental_random_end'] - data['mental_random_start'])
    

if __name__ == '__main__':
    main()

    # 1681326450.993 - 1681325890.426

    # 1681330648.911 - 1681329808.079