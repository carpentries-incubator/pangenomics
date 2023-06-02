---
title: "TDA in Pangenomes"
teaching: 30
exercises: 15
questions:
- "How can I apply TDA to describe Pangenomes"
objectives:
- "Describe Pangenomes using Gudhi"
keypoints:
- "Pangenomes can be described using TDA"
---

~~~
import pandas as pd
import numpy as np
import random

import warnings

warnings.filterwarnings("ignore")  # Deactivate Python warnings

m =5 # The number of genomes (Cols)
n = 10  # Specify the number of genes (rows)
shellsize=int(m/2)+1

core = 0.3  # Probability of being core filling a row with ones
shell = 0.6  # Probability of being shell filling the remaining zeros with ones
cloud =1-core-shell

# Create an empty DataFrame
# Create an empty DataFrame
df = pd.DataFrame(columns=[f'Column{i+1}' for i in range(m)])


# Fill each row with ones based on probabilities P1 and P2
for _ in range(n):
    P1=random.random()  # Generate a random value between 0 and 1

    #print("P1",P1)
    if (core>1-P1): 
      row = [1] * m
      #print("core",row)
  
    elif(core+shell>1-P1): 
        #If the row is not already filled with ones, randomly distribute m/2 ones based on P2
        genome_row = np.random.choice([shellsize, m-1])  # Fill the row with ones based on P1
        #print("number of genomes in which this family is present",genome_row)
        # Create an empty list representing the row
        row = [0] * m
        # Generate k unique random indices
        indices = random.sample(range(m), genome_row)
        # Set the values at the random indices to 1
        for index in indices:
          row[index] = 1
        #print("shell",row)
    else:
        row = [0] * m
        num_ones = 1  # Calculate the number of ones to distribute
        # Choose a random index within the row
        random_index = random.randint(0, m - 1)
        #print("rNDOM INDEX", random_index)  
        # Set the value at the random index to 1
        row[random_index] = 1
        #print("cloud",row)

    # row van a ser las columnas del dataframe
    df = df.append(pd.Series(row, index=df.columns), ignore_index=True)  # Append the row to the DataFrame

print(df)
~~~
{: .python}
