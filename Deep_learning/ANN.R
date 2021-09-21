## case 1
library(data.table)
library(ggplot2)
set.seed(69)
X<-fread("~/Desktop/CADISPAD/process/0920/X.txt")
Y<-fread("~/Desktop/CADISPAD/process/0920/Y.txt")
df<-cbind(X,Y)
df<-df[,-c(1,4)]
colnames(df)<-c("x1","x2","y")
df<-df[-1,]
ggplot(df)+geom_point(aes(x=x1,y=x2,color=y))

## trait test split
train_test_split_index <- 0.8 * nrow(df)
train <- df[1:train_test_split_index,]
test <- df[(train_test_split_index+1): nrow(df),]
X_train <- scale(train[, c(1:2)])

y_train <- train$y
dim(y_train) <- c(length(y_train), 1) # add extra dimension to vector

X_test <- scale(test[, c(1:2)])

y_test <- test$y
dim(y_test) <- c(length(y_test), 1) # add extra dimension to vector
X_train <- as.matrix(X_train, byrow=TRUE)
X_train <- t(X_train)
y_train <- as.matrix(y_train, byrow=TRUE)
y_train <- t(y_train)

X_test <- as.matrix(X_test, byrow=TRUE)
X_test <- t(X_test)
y_test <- as.matrix(y_test, byrow=TRUE)
y_test <- t(y_test)
getLayerSize <- function(X, y, hidden_neurons, train=TRUE) {
  n_x <- dim(X)[1]
  n_h <- hidden_neurons
  n_y <- dim(y)[1]   
  
  size <- list("n_x" = n_x,
               "n_h" = n_h,
               "n_y" = n_y)
  
  return(size)
}
layer_size <- getLayerSize(X_train, y_train, hidden_neurons = 4)
initializeParameters <- function(X, list_layer_size){
  
  m <- dim(data.matrix(X))[2]
  
  n_x <- list_layer_size$n_x
  n_h <- list_layer_size$n_h
  n_y <- list_layer_size$n_y
  
  W1 <- matrix(runif(n_h * n_x), nrow = n_h, ncol = n_x, byrow = TRUE) * 0.01
  b1 <- matrix(rep(0, n_h), nrow = n_h)
  W2 <- matrix(runif(n_y * n_h), nrow = n_y, ncol = n_h, byrow = TRUE) * 0.01
  b2 <- matrix(rep(0, n_y), nrow = n_y)
  
  params <- list("W1" = W1,
                 "b1" = b1, 
                 "W2" = W2,
                 "b2" = b2)
  
  return (params)
}
init_params <- initializeParameters(X_train, layer_size)
lapply(init_params, function(x) dim(x))
sigmoid <- function(x){
  return(1 / (1 + exp(-x)))
}
forwardPropagation <- function(X, params, list_layer_size){
  
  m <- dim(X)[2]
  n_h <- list_layer_size$n_h
  n_y <- list_layer_size$n_y
  
  W1 <- params$W1
  b1 <- params$b1
  W2 <- params$W2
  b2 <- params$b2
  
  b1_new <- matrix(rep(b1, m), nrow = n_h)
  b2_new <- matrix(rep(b2, m), nrow = n_y)
  
  Z1 <- W1 %*% X + b1_new
  A1 <- sigmoid(Z1)
  Z2 <- W2 %*% A1 + b2_new
  A2 <- sigmoid(Z2)
  
  cache <- list("Z1" = Z1,
                "A1" = A1, 
                "Z2" = Z2,
                "A2" = A2)
  
  return (cache)
}
fwd_prop <- forwardPropagation(X_train, init_params, layer_size)
lapply(fwd_prop, function(x) dim(x))
computeCost <- function(X, y, cache) {
  m <- dim(X)[2]
  A2 <- cache$A2
  logprobs <- (log(A2) * y) + (log(1-A2) * (1-y))
  cost <- -sum(logprobs/m)
  return (cost)
}
cost <- computeCost(X_train, y_train, fwd_prop)
cost
backwardPropagation <- function(X, y, cache, params, list_layer_size){
  
  m <- dim(X)[2]
  
  n_x <- list_layer_size$n_x
  n_h <- list_layer_size$n_h
  n_y <- list_layer_size$n_y
  
  A2 <- cache$A2
  A1 <- cache$A1
  W2 <- params$W2
  
  dZ2 <- A2 - y
  dW2 <- 1/m * (dZ2 %*% t(A1)) 
  db2 <- matrix(1/m * sum(dZ2), nrow = n_y)
  db2_new <- matrix(rep(db2, m), nrow = n_y)
  
  dZ1 <- (t(W2) %*% dZ2) * (1 - A1^2)
  dW1 <- 1/m * (dZ1 %*% t(X))
  db1 <- matrix(1/m * sum(dZ1), nrow = n_h)
  db1_new <- matrix(rep(db1, m), nrow = n_h)
  
  grads <- list("dW1" = dW1, 
                "db1" = db1,
                "dW2" = dW2,
                "db2" = db2)
  
  return(grads)
}
back_prop <- backwardPropagation(X_train, y_train, fwd_prop, init_params, layer_size)
lapply(back_prop, function(x) dim(x))
updateParameters <- function(grads, params, learning_rate){
  
  W1 <- params$W1
  b1 <- params$b1
  W2 <- params$W2
  b2 <- params$b2
  
  dW1 <- grads$dW1
  db1 <- grads$db1
  dW2 <- grads$dW2
  db2 <- grads$db2
  
  
  W1 <- W1 - learning_rate * dW1
  b1 <- b1 - learning_rate * db1
  W2 <- W2 - learning_rate * dW2
  b2 <- b2 - learning_rate * db2
  
  updated_params <- list("W1" = W1,
                         "b1" = b1,
                         "W2" = W2,
                         "b2" = b2)
  
  return (updated_params)
}
update_params <- updateParameters(back_prop, init_params, learning_rate = 0.01)
lapply(update_params, function(x) dim(x))
trainModel <- function(X, y, num_iteration, hidden_neurons, lr){
  
  layer_size <- getLayerSize(X, y, hidden_neurons)
  init_params <- initializeParameters(X, layer_size)
  cost_history <- c()
  for (i in 1:num_iteration) {
    fwd_prop <- forwardPropagation(X, init_params, layer_size)
    cost <- computeCost(X, y, fwd_prop)
    back_prop <- backwardPropagation(X, y, fwd_prop, init_params, layer_size)
    update_params <- updateParameters(back_prop, init_params, learning_rate = lr)
    init_params <- update_params
    cost_history <- c(cost_history, cost)
    
    if (i %% 10000 == 0) cat("Iteration", i, " | Cost: ", cost, "\n")
  }
  
  model_out <- list("updated_params" = update_params,
                    "cost_hist" = cost_history)
  return (model_out)
}
#We’re going to train our model, with 40 hidden neurons, for 60000 epochs with a learning rate of 0.9. We will print out the loss after every 10000 epochs.
EPOCHS = 60000
HIDDEN_NEURONS = 40
LEARNING_RATE = 0.9

train_model <- trainModel(X_train, y_train, hidden_neurons = HIDDEN_NEURONS, num_iteration = EPOCHS, lr = LEARNING_RATE)

lr_model <- glm(y ~ x1 + x2, data = train)
lr_pred <- round(as.vector(predict(lr_model, test[, 1:2])))
makePrediction <- function(X, y, hidden_neurons){
  layer_size <- getLayerSize(X, y, hidden_neurons)
  params <- train_model$updated_params
  fwd_prop <- forwardPropagation(X, params, layer_size)
  pred <- fwd_prop$A2
  
  return (pred)
}
y_pred <- makePrediction(X_test, y_test, HIDDEN_NEURONS)
y_pred <- round(y_pred)
tb_nn <- table(y_test, y_pred)
tb_lr <- table(y_test, lr_pred)
calculate_stats <- function(tb, model_name) {
  acc <- (tb[1] + tb[4])/(tb[1] + tb[2] + tb[3] + tb[4])
  recall <- tb[4]/(tb[4] + tb[3])
  precision <- tb[4]/(tb[4] + tb[2])
  f1 <- 2 * ((precision * recall) / (precision + recall))
  
  cat(model_name, ": \n")
  cat("\tAccuracy = ", acc*100, "%.")
  cat("\n\tPrecision = ", precision*100, "%.")
  cat("\n\tRecall = ", recall*100, "%.")
  cat("\n\tF1 Score = ", f1*100, "%.\n\n")
}

### case 2
library(keras)
library(mlbench)
library(dplyr)
library(magrittr)
library(neuralnet)

data("BostonHousing")
data <- BostonHousing
str(data)
#rank order
data %<>% mutate_if(is.factor, as.numeric)
## NN visualization
n <- neuralnet(medv ~ crim+zn+indus+chas+nox+rm+age+dis+rad+tax+ptratio+b+lstat,
               data = data,
               hidden = c(12,7),
               linear.output = F,
               lifesign = 'full',
               rep=1)
plot(n,col.hidden = 'darkgreen',     
     col.hidden.synapse = 'darkgreen',
     show.weights = F,
     information = F,
     fill = 'lightblue')
data <- as.matrix(data)
dimnames(data) <- NULL
## Data Partition
set.seed(123)
ind <- sample(2, nrow(data), replace = T, prob = c(.7, .3))
training <- data[ind==1,1:13]
test <- data[ind==2, 1:13]
trainingtarget <- data[ind==1, 14]
testtarget <- data[ind==2, 14]
str(trainingtarget)
str(testtarget)
## Scaling
m <- colMeans(training)
s <- apply(training, 2, sd)
training <- scale(training, center = m, scale = s)
test <- scale(test, center = m, scale = s)
## Model Creation
model <- keras_model_sequential()
model %>%
  layer_dense(units = 5, activation = 'relu', input_shape = c(13)) %>%
  layer_dense(units = 1)
## Model Compilation
model %>% compile(loss = 'mse',
                  optimizer = 'rmsprop', 
                  metrics = 'mae') 
## Model Fitting
mymodel <- model %>%          
  fit(training,trainingtarget,
      epochs = 100,
      batch_size = 32,
      validation_split = 0.2)
## Prediction
model %>% evaluate(test, testtarget)
pred <- model %>% predict(test)
mean((testtarget-pred)^2) 
## Scatter plot Original vs Predictied
plot(testtarget, pred) 
model %>%
  layer_dense(units = 100, activation = 'relu', input_shape = c(13)) %>%
  layer_dropout(rate=0.4)  %>%
  layer_dense(units = 50, activation = 'relu')  %>%
  layer_dropout(rate=0.2)  %>%
  layer_dense(units = 1)

## case 3
#install.packages('ISLR')
library(ISLR)

print(head(College,2))
## Data preprocessing
# Create Vector of Column Max and Min Values
maxs <- apply(College[,2:18], 2, max)
mins <- apply(College[,2:18], 2, min)

# Use scale() and convert the resulting matrix to a data frame
scaled.data <- as.data.frame(scale(College[,2:18],center = mins, scale = maxs - mins))

# Check out results
print(head(scaled.data,2))
## Train and Test split
# Convert Private column from Yes/No to 1/0
Private = as.numeric(College$Private)-1
data = cbind(Private,scaled.data)

library(caTools)
set.seed(101)

# Create Split (any column is fine)
split = sample.split(data$Private, SplitRatio = 0.70)

# Split based off of split Boolean Vector
train = subset(data, split == TRUE)
test = subset(data, split == FALSE)
## Neural Network Function
feats <- names(scaled.data)

# Concatenate strings
f <- paste(feats,collapse=' + ')
f <- paste('Private ~',f)

# Convert to formula
f <- as.formula(f)

f
# Private ~ Apps + Accept + Enroll + Top10perc + Top25perc + F.Undergrad + 
#P.Undergrad + Outstate + Room.Board + Books + Personal + 
#PhD + Terminal + S.F.Ratio + perc.alumni + Expend + Grad.Rate
#install.packages('neuralnet')
library(neuralnet)
nn <- neuralnet(f,train,hidden=c(10,5),linear.output=FALSE,threshold=0.01)
#The threshold is set to 0.01, meaning that if the change in error during an iteration is less than 1%, then no further optimization will be carried out by the model
## neuron in each hidden layer
## Predictions and Evaluations
# Compute Predictions off Test Set
predicted.nn.values <- neuralnet::compute(nn,test[2:18])
#predicted.nn.values = (predicted.nn.values$net.result * (max(train$Private) - min(train$Private))) + min(train$Private)
# Check out net.result
#print(head(predicted.nn.values$net.result))
#plot(test$Private, predicted.nn.values, col='blue', pch=16, ylab = "predicted rating NN", xlab = "real rating")
#abline(0,1)
# Calculate Root Mean Square Error (RMSE)
RMSE.NN = (sum((test$Private - predicted.nn.values)^2) / nrow(test)) ^ 0.5
##Confusion matrix
predicted.nn.values$net.result <- sapply(predicted.nn.values$net.result,round,digits=0)
table(test$Private,predicted.nn.values$net.result)
##Visulizing the Neural Net
plot(nn)
## Accuracy
results <- data.frame(actual = test$Private, prediction = predicted.nn.values$net.result)
predicted=results$prediction * abs(diff(range(Private))) + min(Private)
actual=results$actual * abs(diff(range(Private))) + min(Private)
comparison=data.frame(predicted,actual)
deviation=((actual-predicted)/actual)
comparison=data.frame(predicted,actual,deviation)
accuracy=1-abs(mean(deviation))
accuracy

##Cross vaildation
## Cross validation of neural network model

# install relevant libraries
install.packages("boot")
install.packages("plyr")

# Load libraries
library(boot)
library(plyr)

# Initialize variables
set.seed(50)
k = 100
RMSE.NN = NULL

List = list( )

# Fit neural network model within nested for loop
for(j in 10:65){
  for (i in 1:k) {
    index = sample(1:nrow(data),j )
    
    trainNN = scaled[index,]
    testNN = scaled[-index,]
    datatest = data[-index,]
    
    NN = neuralnet(rating ~ calories + protein + fat + sodium + fiber, trainNN, hidden = 3, linear.output= T)
    predict_testNN = compute(NN,testNN[,c(1:5)])
    predict_testNN = (predict_testNN$net.result*(max(data$rating)-min(data$rating)))+min(data$rating)
    
    RMSE.NN [i]<- (sum((datatest$rating - predict_testNN)^2)/nrow(datatest))^0.5
  }
  List[[j]] = RMSE.NN
}

Matrix.RMSE = do.call(cbind, List)
## Prepare boxplot
boxplot(Matrix.RMSE[,56], ylab = "RMSE", main = "RMSE BoxPlot (length of traning set = 65)")
## Variation of median RMSE 
install.packages("matrixStats")
library(matrixStats)

med = colMedians(Matrix.RMSE)

X = seq(10,65)

plot (med~X, type = "l", xlab = "length of training set", ylab = "median RMSE", main = "Variation of RMSE with length of training set")
