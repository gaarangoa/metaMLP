
import sys

idx = int(sys.argv[2])
minp = float(sys.argv[3])

predicted = [];

for i in open(sys.argv[1]):
    i=i.split();
    try:
        if float(i[2]) < minp: continue
        predicted.append([i[0].split("|")[idx], i[1].replace('__label__',"").replace('__',"")])
    except:
        print i
        pass

# # Get the total counts for the predicted class:
# predicted_counts = {};
# for i in open(sys.argv[2]):
#     i = i.split();
#     try:
#         predicted_counts[i[1]]+=1;
#     except:
#         predicted_counts[i[1]]=1

# # Get the total counts for each class:
# counts = {}
# for i in open(sys.argv[1]):
#     i = i.split();
#     try:
#         counts[i[1]]+=1;
#     except:
#         counts[i[1]]=1


from sklearn.metrics import classification_report
print(classification_report([i[0] for i in predicted], [i[1] for i in predicted]))



