# Deep learning test
## import python module
import os
from PIL import Image
import numpy as np
import matplotlib.pylpot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
import tensorflow as tf
from tensorflw import keras
from sklearn.model_selection import train_test_split

## import images car
filedir='~/Desktop/Python/cars/'
os.listdir(filedir)

file_list1=[]
for root, dirs, files in os.walk(filedir):
    for file in files:
        if os.path.splitext(file)[1] == ".jpg":
            file_list1.append(os.path.join(root,file))

## change image resolution
for filename in file_list1:
    try:
        im = Image.open(filename)
        new_im = im.resize ((128,128))
        new_im.save('~/Desktop/Python/cars_128/')+filename[12:-4]+'.jpg')
        print('image'+filename[12:-4]+'.jpg'+'resoulution convert finished')
    except OSError as e:
        print (e.args)

## rebuild image list
filedir = '~/Desktop/Python/cars_128/'
os.listdir(filedir)

file_list_1=[]
for root, dirs, files in os.walk(filedir):
    for file in files:
        if os.path.splitext(file)[1] == ".jpg":
            file_list_1.append(os.path.join(root,file))

## import images beauty
filedir='~/Desktop/Python/beauty/'
file_list2=[]
for root, dirs, files in os.walk(filedir):
    for file in files:
        if os.path.splitext(file)[1] == ".jpg":
            file_list2.append(os.path.join(root,file))

## change image resolution
for filename in file_list2:
    try:
        im = Image.open(filename)
        new_im = im.resize ((128,128))
        new_im.save('~/Desktop/Python/beauty_128/')+filename[14:-4]+'.jpg')
        print('image'+filename[14:-4]+'.jpg'+'resoulution convert finished')
    except OSError as e:
        print (e.args)

## rebuild image list
filedir = '~/Desktop/Python/beauty_128/'
os.listdir(filedir)

file_list_2=[]
for root, dirs, files in os.walk(filedir):
    for file in files:
        if os.path.splitext(file)[1] == ".jpg":
            file_list_2.append(os.path.join(root,file))

## merge data in the list
len(file_list_1)
len(file_list_2)
file_list_all = file_list_1 + file_list_2
len(file_list_all)

## convert image to array
M = []
for filename in file_list_all:
    im = Image.open(filename)
    width,height = im.size
    im_L = im.convert("L")
    Core = im_L.getdata()
    arr1 = np.array(Core,dtype='float32') / 255.0
    arr1.shape
    list_img = arr1.tolist()
    M.extend(list_img)

X = np.array(M).reshape(len(file_list_all),width,height)
X.shape

class_names = ['cars', 'beauty']

## store image into dict
dict_label = {0:'cars', 1:'beauty'}
print(dict_label[0])
print(dict_label[1])

## use list imput label
label = [0]*len(file_list_1) + [1]*len(file_list_2)
y = np.array(label)

## set train/test 4:1
train_images, test_images, train_labels, test_labels = train_test_split(
    X, y, test_size=0.2, random_state=0)

plt.figure()
plt.imshow(train_images[10])
plt.colorbar()
plt.grid(False)

#show top 25 and label; validation
plt.figure(figsize=(10,10))
for i in range(25):
    plt.subplot(5,5,i+1)
    plt.xticks([])
    plt.yticks([])
    plt.grid(False)
    plt.imshow(train_images[i],cmap=plt.cm.binary)
    plt.xlabel(class_names[train_labels[i]])

#1st layer: 128 nodes
#2nd layer(last): 2 nodes softmax layer - return 2 probability array, sum = 1
#every node includes a percentage, represent probability
model = keras.Sequential([
    keras.layers.Flatten(input_shape=(128,128)),
    keras.layers.Dense(128, activation=tf.nn.relu),
    keras.layers.Dense(2, activation=tf.nn.softmax)
])

#################################################################################
model.compile(optimizer=tf.train.AdamOptimizer(),
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])

model.fit(train_images,train_labels,epochs=5)
test_loss, test_acc = model.evaluate(test_images, test_labels)
print('Test accuracy:', test_acc)

# prediction
predictions = model.predict(test_images)
predictions[0]
np.argmax(predictions[0])
dict_label[np.argmax(predictions[0])]
predictions = model.predict(test_images)
predictions[0]
np.argmax(predictions[0])
dict_label[np.argmax(predictions[0])]

#define plot function
def plot_image(i, predictions_array, true_label, img):
    predictions_array, true_label, img = predictions_array[i], true_label[i], img[i]
    plt.grid(False)
    plt.xticks([])
    plt.yticks([])

    plt.imshow(img, cmap=plt.cm.binary)

    predicted_label = np.argmax(predictions_array)
    if predicted_label == true_label:
        color = '#00bc57'
    else:
        color = 'red'

    plt.xlabel("{} {:2.0f}% ({})".format(class_names[predicted_label],
                                         100*np.max(predictions_array),
                                         class_names[true_label]),
                                         color=color)

#############################################################################################
def plot_value_array(i, predictions_array, true_label):
    predictions_array, true_label = predictions_array[i], true_label[i]
    plt.grid(False)
    plt.xticks([])
    plt.yticks([])
    thisplot = plt.bar(range(len(class_names)), predictions_array,
                        color='#FF7F0E', width=0.2)
    plt.ylim([0,1])
    predicted_label = np.argmax(predictions_array)

    thisplot[predicted_label].set_color('red')
    thisplot[true_label].set_color('#00bc57')

################################################################################################
# try first image
i = 0
plt.figure(figsize=(6,3))
plt.subplot(1,2,1)
plot_image(i, predictions,test_labels, test_images)
plt.subplot(1,2,2)
plot_value_array(i, predictions, test_labels)

# try 12th image
i = 12
plt.figure(figsize=(6,3))
plt.subplot(1,2,1)
plot_image(i, predictions,test_labels, test_images)
plt.subplot(1,2,2)
plot_value_array(i, predictions, test_labels)

# probability boxplot; green=true; red=false
num_rows = 5
num_cols = 3
num_images = num_rows*num_cols
plt.figure(figsize=(2*2*num_cols, 2*num_rows))
for i in range(num_images):
    plt.subplot(num_rows, 2*num_cols, 2*i+1)
    plot_image(i, predictions, test_labels, test_images)
    plt.subplot(num_rows, 2*num_cols, 2*i+2)
    plot_value_array(i, predictions, test_labels)

# 15th image
img = test_images[14]
plt.imshow(img, cmap=plt.cm.binary)
#############################################################################################
# add member
img = (np.expand_dims(img,0))
print(img.shape)
# predict image:
predictions_single = model.predict(img)
print(predictions_single)
###############################################################################################
plot_value_array(0, predictions_single, test_labels)
_ = plt.xticks(range(2), class_names, rotation=45)
###############################################################################################
np.argmax(predictions_single[0])
dict_label[np.argmax(predictions_single[0])]

# new file for prediction
filedir="new/pred"
file_list_pred=[]
for root, dirs, files in os.walk(filedir):
    for file in files:
        if os.path.splitext(file)[1] == '.jpg':
            file_list_pred.append(os.path.join(root, file))
#change resolution
for filename in file_list_pred:
    try:
        im = Image.open(filename)
        new_im = im.resize((128,128))
        new_im.save('new/pred_128'+filename[12:-4]+'.jpg')
        print('new'+filename[12:-4]+'.jpg'+'convert finished')
    except OSError as e:
        print(e.args)

#get image list
filedir='new/pred_128/'
file_list_pred = []
for root, dirs, files in os.walk(filedir):
    for file in files:
        if os.path.splitext(file)[1] == '.jpg':
            file_list_pred.append(os.path.join(root, file))

# for one image
im = Image.open(file_list_pred[3])
width,height = im.size
im_L = im.convert("L")
Core = im_L.getdata()
arr1 = np.array(Core, dtype='float32')/255.0
arr1.shape
list_img = arr1.tolist()
pred_images = np.array(list_img).reshape(width,height)
pred_labels = np.array([0])
img = pred_images
print(img.shape)
plt.imshow(img, cmap=plt.cm.binary)
# add member
img = (np.expand_dims(img,0))
pring(img.shape)
###################################################################################
img = (np.expand_dims(img,0))
print(img.shape)
###################################################################################
# prediction
predictions_single = model.predict(img)
print(predictions_single)
###################################################################################
plot_value_array(0,predictions_single,pred_labels)
_ = plt.xticks(range(2), class_names, rotation=45)
###################################################################################
np.argmax(predictions_single[0])
dict_label[np.argmax(predictions_single[0])]
###################################################################################
# multiple image
N = []
for filename in file_list_pred:
    im = Image.open(filename)
    im_L = im.convert("L")
    Core = im_L.getdata()
    arr1 = np.array(Core, dtype='float32')/255.0
    arr1.shape
    list_img = arr1.tolist()
    N.extend(list_img)

pred_images = np.array(N).reshape(len(file_list_pred),width,height)
pred_labels = [0,0,0,0,0,1,1,1,1,1]

#prediction
predictions = model.predict(pred_images)

for i in range(len(file_list_pred)):
    img = pred_images[i]
    plt.imshow(img, cmap=plt.cm.binary)
    # add member
    img = (np.expand_dims(img,0))
    # prediction
    predictions_single = model.predict(img)
    ###################################################################################
    plot_value_array(0, predictions_single, pred_labels)
    _ = plt.xticks(range(2), class_names, rotation=45)
    print (str(i)+'th image is: '+dict_label[np.argmax(predictions_single[0])])
