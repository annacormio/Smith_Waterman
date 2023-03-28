import pandas as pd
import numpy as np

#INPUT SEQUENCES
A='AAATCGTACAT'
B='CCTGTACAA'

'''
A= input("insert the first sequence")
B= input("insert the second sequence")
'''

#SCORING SCHEME
#dafault scoring scheme
match= 1
mismatch= -1
gap= -2
#user can perzonalize scoring scheme
dec=input("Do u want to change the scoring scheme? (If not just press enter and I will use the default one)")
if dec:
    match=int(input("match score = "))
    mismatch = int(input("mismatch score = "))
    gap = int(input("gap score = "))

#INITIALIZING MATRIX
df=pd.DataFrame(columns=list(' '+A), index=list(' '+B)) #empty matrix leave a space at 1st column and row
tb= pd.DataFrame(columns=list(' '+A), index=list(' '+B)) #initializing new frame for traceback tracking

#initialize first row to 0
for c in range(len(A)+1): #since we added a space to build the matrix correctly we need to add 1 to the length of the string in order to reach the end
    df.iloc[0, c] =0
    tb.iloc[0, c]=0

#initialize first column to 0
for r in range(len(B)+1):
    df.iloc[r, 0]= 0
    tb.iloc[r, 0] = 0



#MATRIX FILLING AND TRACEBACK

#iterations
for r in range (1,len(B)+1): #iterate on each row with nt.
    for c in range(1,len(A)+1): #iterate on each column with nt.

        #MATRIX FILLING
        s_r_gap = df.iloc[r, c - 1] + gap  # right gap
        s_u_gap = df.iloc[r - 1, c] + gap  # up gap
        if df.columns[c]==df.index[r]: #match situation
            s_match = df.iloc[r-1,c-1] + match   #r-1,c-1 is the diagonal number found up left of the computed one
            score = max(s_match, s_r_gap, s_u_gap,0) #select the max score among the possible one and assign it to the cell
            if score == s_match:
                tb.iloc[r, c] = 'diag'#traceback matrix
        else: #mismatch situation
            s_mismatch =  df.iloc[r-1,c-1] + mismatch #mismatch in diagonal
            score = max(s_mismatch, s_r_gap, s_u_gap,0) #select the max score among the possible one and assign it to the cell
            if score == s_mismatch:
                tb.iloc[r, c] = 'diag' #traceback matrix
        df.iloc[r, c] = score

        if score == s_r_gap:
            tb.iloc[r, c] = 'left'
        elif score == s_u_gap:
            tb.iloc[r, c] = 'up'

        elif tb.iloc[r,c]!='diag':
            tb.iloc[r,c] = 'stop'


#SEQUENCE ALIGNMENT FROM TRACEBACK MATRIX

a1='' #first sequence
mid='' #string of alignments
a2='' #second sequence

#setting starting point for traceback
best_score=df.max(axis=1).max() # finding max value in dataframe
mask = df.applymap(lambda x: True if x == best_score else False).to_numpy() #creating a mask dataframe and convert it to numpy (useful to apply the following funtion)
indexes = np.argwhere(mask)[0]#return list of indexes of the max value (where mask ==True)
print(r)
print(c)
#traceback
while tb.iloc[r,c] != 'stop':
    #if there is a match or mismatch I align the 2 nt.
    if tb.iloc[r,c]=='diag':
        a1 = a1 + A[c - 1]
        a2 = a2 + B[r - 1]
        if A[c-1]== B[r-1]:
            mid = mid+('|')
        else:
            mid = mid + ('.')
        c = c - 1
        r = r - 1


    #gaps
    elif tb.iloc[r,c]=='up':
        a1 = a1 + '-' #gap in the first sequence
        mid = mid + (' ')
        a2=a2+B[r-1]
        r =r - 1

    elif tb.iloc[r,c]=='left':
        a1=a1+A[c-1]
        mid = mid + (' ')
        a2 = a2 + '-' #gap in the second sequence
        c =c - 1



#OUTPUT
print(f'score matrix \n {df}') #printing the scored matrix
print(f'\n the score of the alignment is {score} \n') #printing score of the alignment
print(tb) #printing the matrix with traceback indications
#since I traceBACK the first alignment i get is of the sequences displayed backward
align=f'{a1[::-1]}\n{mid[::-1]}\n{a2[::-1]}' #using [::-1] i am reversing the strings to have the right alignment, not backward
print(align) #printing the 2 aligned sequences


#ATGTACATAGA
#GTATCGTA


